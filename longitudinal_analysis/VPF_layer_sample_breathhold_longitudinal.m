%% load data
clear
clc

% some paramters
subid = '7218';
sessionid = '01';
TE = [2 6 10 14 17 20];
numecho = numel(TE);
vols = 19;
run = 1:2;
Tvols = numel(run)*vols;
res = [0.75 0.75 0.75]; %resolution of the anatomical!
sz  = [228 224 40 vols numecho];

N = [20 10]; %number of layers

use_echo = 3; %echo on which the realignment/coregistration should be estimated with

task_key_val = 'task_breathhold';
part_key = 'mag';
contrast_key_val = 'part-mag_bold';

% this is the overall path
globalpath  = '/media/pfaffenrot/PostDoc_data/projects/hippocampus_breathhold/FLASH/derivates/';

% this is the path where the actual data are saved.
currentpath = ['/pipeline/' subid '/ses-' sessionid '/func_bart'];

% this is the path where the freesurfer output is stored
structpath = [globalpath '/pipeline/' subid '/ses-01/anat/T1/presurf_MPRAGEise/presurf_UNI/'];
% structpath = '/media/pfaffenrot/Elements/postdoc/projects/hippocampus_VASO/derivatives/pipeline/7513/ses-01/anat/T1/presurf_MPRAGEise/presurf_UNI/';

% this is the path of the anatomical reference
% fspath = [globalpath 'freesurfer/' subid];

% this is the path where the non-layer stats like b-maps, additional nuisance regressors etc. are saved
statspath = [globalpath '/pipeline/' subid '/ses-' sessionid '/stats'];

% this is the path where all the layer-related data (sampled functionals, sampbled b-map etc.) are saved
layerpath = [globalpath '/pipeline/' subid '/ses-' sessionid '/layerfication'];


data = cell(numecho,numel(run));
for kk = 1:numel(run)
    for ii = 1:numecho
        data {ii,kk} = [globalpath currentpath '/run' num2str(kk) '/' part_key '/sub-' subid '_' task_key_val '_run-' sprintf('%02.0f',kk) '_' contrast_key_val '_echo-' num2str(ii) '.nii'];
        if ii == 1 && kk == 1
            tmp = spm_vol(data{ii,kk});
            hdr = repmat(tmp, [1 numecho numel(run)]);
            hdr(:,ii,kk) = tmp;
            clear tmp;
        else
            hdr(:,ii,kk) = spm_vol(data{ii,kk});
        end
    end
end


pos = strfind(data{1,1},'/');
pos = pos(end);

%--> read structural data
mprage = [structpath 'sub-' subid '_UNI_MPRAGEised_biascorrected_denoised.nii'];
mprage_mask = [structpath 'sub-' subid '_UNI_MPRAGEised_brainmask.nii'];
% mprage = [structpath 'sub-7513_UNI_MPRAGEised_biascorrected_denoised.nii'];
% mprage_mask = [structpath 'sub-7513_UNI_MPRAGEised_brainmask.nii'];
hdr_mprage = spm_vol(mprage);
dims       = hdr_mprage.dim;            
%<--

surf_path = ['/media/pfaffenrot/PostDoc_data/projects/hippocampus_breathhold/FLASH/derivates/hippunfold/' subid];

TR = 40.992;
bs = 3;
hrf_flags = struct('TR',TR,'name','hrf','length',bs*TR,'order',1,'T',4,'vols',vols);
hpf_flags = struct('TR',TR,'row',1:vols,'HParam',2*bs*TR);

paradigm  = [repmat([0.75;0;1], [floor(vols/bs) 1])];
paradigm = [paradigm; 0];

% layer sampling
rule = '1 2 3 4 5 6';
subfields = numel(double(str2num(rule)));
if subfields == 0
    subfields = 1;
end

for jj = 1:numel(run)
    list = dir([globalpath currentpath '/run' num2str(jj) '/' part_key '/rsub-*_0*.nii']);
    sampled_img = cell(size(list));

    for ii = 1:length(list)
        sampled_img{ii} = [list(ii).folder '/' list(ii).name];
    end


    if jj == 1
        tmp = VPF_create_hippocampus_layers(sampled_img,str2double(subid),surf_path,mprage,N,rule);
        layers = cell([size(tmp),numel(run)]);
        layers(:,:,jj) = tmp;
    else
        layers(:,:,jj) = VPF_create_hippocampus_layers(sampled_img,str2double(subid),surf_path,mprage,N,rule);
    end

    %reshape into layers x volumes x echoes
    for hemi = 1:2
        for ss = 1:subfields
            if ss < 5
                layers{hemi,ss,jj} = reshape(layers{hemi,ss,jj},[size(layers{hemi,ss,jj},1) sum(N) vols numecho]); 
            else
                layers{hemi,ss,jj} = reshape(layers{hemi,ss,jj},[size(layers{hemi,ss,jj},1) 1 vols numecho]); 
            end
        end
    end
end

%runs = 5th dimension
layers = cellfun(@(x, y) cat(5, x, y), layers(:,:,1), layers(:,:,2), 'UniformOutput', false);

for uu = 1:size(layers,1)
    for vv = 1:size(layers,2)
    %highpass filter
    for ii = 1:numecho
        for jj = 1:numel(run)
            layers{uu,vv}(:,:,:,ii,jj) = VPF_HPF_SPM(layers{uu,vv}(:,:,:,ii,jj),hpf_flags);
        end
    end

    %average over runs and remove 1st and last volume due to highpass filter artifact
%     layers{uu,vv} =  mean(layers{uu,vv},5);
    layers{uu,vv} = layers{uu,vv}(:,:,:,:,2);
    layers{uu,vv} = squeeze(layers{uu,vv}(:,:,2:end-1,:));
    end
end

for hemi = 1:2
    for subfield = 5:6
        sz = size(layers{hemi,subfield});
        layers{hemi,subfield} = reshape(layers{hemi,subfield}, [sz(1) 1 sz(2:end)]);
    end
end


outpath = ['/media/pfaffenrot/PostDoc_data/projects/data/' subid '/breathhold'];

paradigm_use = paradigm(2:end-1);
dS_breathhold = layers;
for hemi = 1:2
    for subfield = 1:6

        OFF = mean(layers{hemi,subfield}(:,:,paradigm_use==0,:),3);
        ON = mean(layers{hemi,subfield}(:,:,paradigm_use==1,:),3);
        if subfield < 5
            OFF = squeeze(OFF);
            ON = squeeze(ON);
            sz = size(OFF);
            OFF = reshape(OFF, [sz(1:2) 1 sz(end)]);            
            ON = reshape(ON, [sz(1:2) 1 sz(end)]);
        end

        dS_breathhold{hemi,subfield} = cat(3,OFF,ON);
    end
end

save([outpath '/dS_breathhold.mat'],'dS_breathhold');