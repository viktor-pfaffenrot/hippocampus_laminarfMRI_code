function varargout = VPF_create_hippocampus_layers(sampled_img,subject_id,hippunfold_path,T1_path,N_layers,rule,SWI_file)
% This function transforms surfaces generated by hippunfold into the space of
% underlying T1 weighted image and equidistantly samples the input images
% using the midthickness surface generated by hippunfold as a landmark.

% The layers can be masked according to the subfield definitions provided
% by hippunfold. If an SWI image is provided as mask, the index of the
% vertices which should be summed to reduce venous bias are saved in the
% 2nd output variable (e.g. layers = squeeze(mean(layers(varargout(2),:),1)))

%Relies on some SPM12 functions. Make sure SPM12 is in your MATLAB path.
% Compressed .nii is supported

%INPUT:
%sampled_img  [cell, 1xN_volumes or 1x1 for 4D nifti]
%                       : cell of strings containing full path and filename
%                         of the nifti images which need to be sampled
%subject_id   [int]     : subject ID used in hippunfold. The surfaces are
%                         read in using that ID
%hippunfold_path [str]  : full path to the hippunfold output, e.g.
%                         ./hippunfold/sub-xxxx/
%T1_path      [str]     : full path and filename of the T1 weighted nifti

%OPTIONAL:
%N_layers     [int int] : number of sampling point, i.e. layers. Default
%                         [20 10]. The first number corresponds to layers
%                         sampled within the hippocampus. The 2nd numer are
%                         the layers sampled in SRLM.
%rule         [str]     : rule to select several subfields(1=subiculum, 2=CA1
%                       : 3=CA2,4=CA3,5=CA4, 6=DG),  e.g.rule = '1 2 3'.
%                         Make sure the string is sparated by spaces.
%                         DG is special as is is only one layer so it is
%                         treated seperately in the code. If you want the
%                         entire hippocampus, average over subfieds after
%                         sampling
%SWI_file     [str]     : full path and filename of the vessel masked SWI
%                         image in anatomy space (e.g. rSWI_vessel_masked.nii).


%OUTPUT:
%layers{2,N_subfields}[N_vertices,N_layers, N_volumes]
%                        : Cell containing the layers for each
%                          hemisphere (1st dim) and for each subfield
%                          (2nd dim).
%idx [N_vertices to use] : array of indices to average over to reduce venous
%                          bias provided that a SWI file was given



%TODO:
%
%

% Viktor Pfaffenrot <viktor.pfaffenrot@uni-due.de>; October 2023, ELH

if nargin < 7
    SWI_file = [];
    SWI_mask = [];
else
    SWI_thres = 3;
    tmp = load_nifti(SWI_file).vol;
    SWI_mask = zeros(size(tmp));
    SWI_mask(tmp>SWI_thres) = tmp(tmp>SWI_thres)./max(tmp(tmp>SWI_thres))*100*(-1);
end

if nargin < 6
    rule = '';
end

if nargin < 5 || isempty(N_layers)
    %the 2nd number is the number of layers outside of the inner surface
    %to also sample the SRLM, i.e. 20 layers between inner and outer and 10
    %layers outside the inner
    N_layers = [20 10];
end

DEBUG = false;

surf_path = [hippunfold_path '/surf'];

%load boundaries and subfield_labels
[layer_boundaries,subfield_labels] = VPF_load_hippocampus_layer_boundaries(surf_path,subject_id);

%load first image into memory. If it is a 4D dataset, its the only image in the cell
input_image = load_nifti(sampled_img{1}).vol;
if ndims(input_image) > 3
    N_volumes = size(input_image,4);
else
    N_volumes = numel(sampled_img);
end

mask = split(rule);

for volume = 1:N_volumes
    if ndims(input_image) == 3 && volume ~= 1 % we have a set of 3D images and the first image in the cell has already been loaded. No need to load again
        input_image = load_nifti(sampled_img{volume}).vol;
    end

    for jj = 1:2 %hemispheres
        %create corresponding boundaries and sample between them
        for subfield = 1:length(mask)
            single_rule = mask{subfield};
            layer_boundaries_sampled = VPF_transform_layers_to_matlab_space(layer_boundaries(jj,:),T1_path,...
                subfield_labels{jj},single_rule,DEBUG);

            if strcmp(single_rule,'6')
                layer_boundaries_sampled = layer_boundaries_sampled{2};
                CA4DG = 2;
            elseif strcmp(single_rule,'5')
                layer_boundaries_sampled = layer_boundaries_sampled{1};
                CA4DG = 1;
            else
                layer_boundaries_sampled = layer_boundaries_sampled{1};
                CA4DG = 0;
            end

            if ndims(input_image) > 3
                [tmp_layers,tmp_idx] = VPF_sample_layers(input_image(:,:,:,volume),layer_boundaries_sampled,N_layers,CA4DG,SWI_mask);
            else
                [tmp_layers,tmp_idx] = VPF_sample_layers(input_image,layer_boundaries_sampled,N_layers,CA4DG,SWI_mask);
            end


            if volume == 1 % create layer cell array. Can do it only now, because now I know how many vertices exist
                if subfield == 1 && jj == 1
                    layers = cell(2,length(mask));
                    idx = layers;
                end
                layers{jj,subfield} = zeros([size(tmp_layers),N_volumes]);
                idx{jj,subfield} = zeros(size(tmp_idx));
                if any(ismember(mask,'6'))
                    pos = ismember(mask,'6');
                    layers{jj,pos} = zeros([size(tmp_layers),N_volumes]);
                end
            end
            layers{jj,subfield}(:,:,volume) = tmp_layers;
            idx{jj,subfield} = tmp_idx;
        end
    end
end

varargout{1} = layers;
varargout{2} = idx;
end

function [layer_boundaries, subfield_labels] = VPF_load_hippocampus_layer_boundaries(surf_path,subject_id)
%loads the hippunfold layer boundaries for both hemispheres using the gifti tools provided by
%hippunfold.
%https://github.com/jordandekraker/hippunfold_toolbox

%layer_boundaries and subfield_labels are cells of size (2,1) corresponding to left/right hemisphere
%dimension of the layer_boundaries matrix is : [N_vertices x/y/z inner/middle/outer]
%dimension of the subfield_labels is : [N_vertices 1]

gifti_files = {[surf_path '/sub-' sprintf('%04.0f',subject_id) '_hemi-L_space-T2w_den-0p5mm_label-hipp_inner.surf.gii'],...
    [surf_path '/sub-' sprintf('%04.0f',subject_id) '_hemi-L_space-T2w_den-0p5mm_label-hipp_midthickness.surf.gii'],...
    [surf_path '/sub-' sprintf('%04.0f',subject_id) '_hemi-L_space-T2w_den-0p5mm_label-hipp_outer.surf.gii'],...
    [surf_path '/sub-' sprintf('%04.0f',subject_id) '_hemi-L_space-T2w_den-0p5mm_label-dentate_midthickness.surf.gii'];...
    [surf_path '/sub-' sprintf('%04.0f',subject_id) '_hemi-R_space-T2w_den-0p5mm_label-hipp_inner.surf.gii'],...
    [surf_path '/sub-' sprintf('%04.0f',subject_id) '_hemi-R_space-T2w_den-0p5mm_label-hipp_midthickness.surf.gii'],...
    [surf_path '/sub-' sprintf('%04.0f',subject_id) '_hemi-R_space-T2w_den-0p5mm_label-hipp_outer.surf.gii'],...
    [surf_path '/sub-' sprintf('%04.0f',subject_id) '_hemi-R_space-T2w_den-0p5mm_label-dentate_midthickness.surf.gii']};

layer_boundaries = cell(2,2);
for jj = 1:2
    for ii = 1:size(gifti_files,2)-1

        if ii == 1 && jj == 1
            tmp = gifti(gifti_files{jj,ii});
            N_vertices = size(tmp.vertices,1);
            layer_boundaries{1,1} = zeros(N_vertices,3,3);
            layer_boundaries{1,1}(:,:,ii) = tmp.vertices;
        elseif ii == 1 && jj == 2
            tmp = gifti(gifti_files{jj,ii});
            N_vertices = size(tmp.vertices,1);
            layer_boundaries{2,1} = zeros(N_vertices,3,3);
            layer_boundaries{2,1}(:,:,ii) = tmp.vertices;
        else
            tmp = gifti(gifti_files{jj,ii});
            layer_boundaries{jj,1}(:,:,ii) = tmp.vertices;
        end
    end

    tmp = gifti(gifti_files{jj,end});
    layer_boundaries{jj,2} =double(tmp.vertices);
end



%load the labels
label_files = {[surf_path '/sub-' sprintf('%04.0f',subject_id) '_hemi-L_space-T2w_den-0p5mm_label-hipp_atlas-bigbrain_subfields.label.gii']...
    [surf_path '/sub-' sprintf('%04.0f',subject_id) '_hemi-R_space-T2w_den-0p5mm_label-hipp_atlas-bigbrain_subfields.label.gii']};

subfield_labels = cell(2,1);
for jj = 1:2
    tmp = gifti(label_files{jj});
    subfield_labels{jj} = tmp.cdata;
end

end

function layer_boundaries = VPF_transform_layers_to_matlab_space(layer_boundaries,T1_path,subfield_labels,rule,DEBUG)

sz = size(layer_boundaries{1});
%--> read structural data
hdr_T1 = spm_vol(T1_path);
%<--

%mask layers based on subfields and the desired rule
rule_split = split(rule);
if ~strcmp(rule,'') && ~any(ismember(rule_split,'6'))
    cmd = '';
    for ii = 1:length(rule_split)
        if ~strcmp(rule_split{ii},'&')
            cmd = [cmd 'subfield_labels==' rule_split{ii}];
        else
            cmd = [cmd ' | '];
        end
    end
    layer_boundaries{1,1} = layer_boundaries{1,1}(eval(cmd),:,:);
end

%--> bring surfaces into matlab space
T1_mat = spm_get_space(T1_path);
boundaries_2_cat = ones(size(layer_boundaries{1,1}));
boundaries_2_cat = boundaries_2_cat(:,1,:);
layer_boundaries{1,1} = cat(2,layer_boundaries{1,1},boundaries_2_cat);

for ii = 1:sz(end)
    tmp = T1_mat\squeeze(layer_boundaries{1,1}(:,:,ii)).';
    layer_boundaries{1,1}(:,:,ii) = tmp.';
end

layer_boundaries{1,1} = layer_boundaries{1,1}(:,1:3,:);

%bring surface of DG into matlab space
if any(ismember(rule_split,'6'))
    boundaries_2_cat = ones(size(layer_boundaries{1,2}));
    boundaries_2_cat = boundaries_2_cat(:,1,:);
    layer_boundaries{1,2} = cat(2,layer_boundaries{1,2},boundaries_2_cat);

    tmp = T1_mat\squeeze(layer_boundaries{1,2}).';
    layer_boundaries{1,2} = tmp.';

    layer_boundaries{1,2} = layer_boundaries{1,2}(:,1:3);
end
%<--


if DEBUG==true
    %plot the surfaces onto the T1
    img_T1 = spm_read_vols(hdr_T1);
    for slice = 69
        idx = find(layer_boundaries{1,1}(:,3,1)>slice & layer_boundaries{1,1}(:,3,1)<slice+1);
        idx2 = find(layer_boundaries{1,1}(:,3,3)>slice & layer_boundaries{1,1}(:,3,3)<slice+1);
        idx3 = find(layer_boundaries{1,1}(:,3,2)>slice & layer_boundaries{1,1}(:,3,2)<slice+1);

        idx4 = find(layer_boundaries{1,2}(:,3)>slice & layer_boundaries{1,2}(:,3)<slice+1);
        figure;imagesc(squeeze(img_T1(:,:,slice)), [0 150]), colormap(gray(256)), title(['Slice ' num2str(slice)])
        hold on;

        plot(layer_boundaries{1,1}(idx,2,1),layer_boundaries{1,1}(idx,1,1),'color','r','Marker','.','MarkerSize',10,'LineStyle','none');
        plot(layer_boundaries{1,1}(idx3,2,2),layer_boundaries{1,1}(idx3,1,2),'color','y','Marker','.','MarkerSize',10,'LineStyle','none');
        plot(layer_boundaries{1,1}(idx2,2,3),layer_boundaries{1,1}(idx2,1,3),'color',[17,249,249]/255,'Marker','.','MarkerSize',10,'LineStyle','none');
        plot(layer_boundaries{1,2}(idx4,2),layer_boundaries{1,2}(idx4,1),'w.');

        legend('inner','midthickness','outer','DG');
    end
end
end

function [layers_out,idx] = VPF_sample_layers(sampled_img,layer_boundaries,N_layers,CA4DG,SWI_mask)

SAMPLING_ORDER = 3;
N_vertices = size(layer_boundaries,1);

switch CA4DG
    case 1 % CA4
        layers = zeros(N_vertices,3);
        SWI_sampled = layers;
        for k = 1:N_vertices
            X = squeeze(layer_boundaries(k,1,:));
            Y = squeeze(layer_boundaries(k,2,:));
            Z = squeeze(layer_boundaries(k,3,:));
            layers(k,:) = spm_sample_vol(sampled_img,X,Y,Z,SAMPLING_ORDER);

            if ~isempty(SWI_mask)
                SWI_sampled(k,:) = spm_sample_vol(SWI_mask,X,Y,Z,SAMPLING_ORDER);
            end
        end
        layers = mean(layers,2);
        SWI_sampled = mean(SWI_sampled,2);
    case 2 %DG
        layers= zeros(N_vertices,1);
        SWI_sampled = layers;
        for k = 1:N_vertices
            layers(k,:) = spm_sample_vol(sampled_img,layer_boundaries(k,1),layer_boundaries(k,2),layer_boundaries(k,3),SAMPLING_ORDER);
        end
        if ~isempty(SWI_mask)
            for k = 1:N_vertices
                SWI_sampled(k,:) = spm_sample_vol(SWI_mask,layer_boundaries(k,1),layer_boundaries(k,2),layer_boundaries(k,3),SAMPLING_ORDER);
            end
        end
    otherwise %all other subfields
        layers= zeros(N_vertices,sum(N_layers));
        SWI_sampled = layers;
        for k = 1:N_vertices

            %take the distance between the inner and the midthickness layer
            half_thick_X = layer_boundaries(k,1,2)-layer_boundaries(k,1,1);

            %create linear space between inner and midthickness of half the layer number
            X  = linspace(layer_boundaries(k,1,1),layer_boundaries(k,1,2),floor(N_layers(1)/2));

            %create linear space between inner layer and inner layer minus half the cortical thickness and midthickness
            %to sample in SRLM
            Xout = linspace(layer_boundaries(k,1,1)-half_thick_X,layer_boundaries(k,1,1),N_layers(2)+1);

            %concatenate inner to mid space with mid to outer space
            X = cat(2,X,linspace(layer_boundaries(k,1,2),layer_boundaries(k,1,3),floor(N_layers(1)/2)+1));
            %concatenate SRLM sampling space with sampling space within cortex
            X = cat(2,Xout,X);

            %remove doubled connection points between subspaces
            X = X([1:N_layers(2) N_layers(2)+2:floor(N_layers(1)/2)+N_layers(2) N_layers(2)+floor(N_layers(1)/2)+2:end]);


            half_thick_Y = layer_boundaries(k,2,2)-layer_boundaries(k,2,1);
            Y    = linspace(layer_boundaries(k,2,1),layer_boundaries(k,2,2),floor(N_layers(1)/2));
            Yout = linspace(layer_boundaries(k,2,1)-half_thick_Y,layer_boundaries(k,2,1),N_layers(2)+1);
            Y    = cat(2,Y,linspace(layer_boundaries(k,2,2),layer_boundaries(k,2,3),floor(N_layers(1)/2)+1));
            Y    = cat(2,Yout,Y);
            Y    = Y([1:N_layers(2) N_layers(2)+2:floor(N_layers(1)/2)+N_layers(2) N_layers(2)+floor(N_layers(1)/2)+2:end]);


            half_thick_Z = layer_boundaries(k,3,2)-layer_boundaries(k,3,1);
            Z    = linspace(layer_boundaries(k,3,1),layer_boundaries(k,3,2),floor(N_layers(1)/2));
            Zout = linspace(layer_boundaries(k,3,1)-half_thick_Z,layer_boundaries(k,3,1),N_layers(2)+1);
            Z    = cat(2,Z,linspace(layer_boundaries(k,3,2),layer_boundaries(k,3,3),floor(N_layers(1)/2)+1));
            Z    = cat(2,Zout,Z);
            Z    = Z([1:N_layers(2) N_layers(2)+2:floor(N_layers(1)/2)+N_layers(2) N_layers(2)+floor(N_layers(1)/2)+2:end]);

            layers(k,:) = spm_sample_vol(sampled_img,X,Y,Z,SAMPLING_ORDER);

            if ~isempty(SWI_mask)
                SWI_sampled(k,:) = spm_sample_vol(SWI_mask,X,Y,Z,SAMPLING_ORDER);
            end
        end
end

if ~isempty(SWI_mask)
    [~,mpos] = min(mean(SWI_sampled,1),[],2);
    idx = (find(SWI_sampled(:,mpos)>0));
else
    idx = 1:N_vertices;
end
layers_out = layers;
end
