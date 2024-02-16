function VPF_create_surface_Tmap_results_memory(layers,SPM,outpath,con_idx,hippunfold_path,FWHM,idx,ZTRANS)
% This function calculates statistical t-maps (FDR-corrected, alpha = 0.05)
% and projects them onto the inner and outer surface of the unfolded 
% hippocampus (averaged over hemispheres) created by hippunfold. The 
% surface is saved as .gii.

%INPUT:
%layers   [cell]     : cell of size (hemis,subfields,runs) containing the 
%                      sampled data (size = (N_vertices,N_layers,N_volumes).
%                      ASSUMES DG TO BE IN THERE, I.E. SUBFIELDS = 6!
%                      Output of "VPF_layer_pipeline_hippocampus_EPI.m"
%SPM      [struct]   : Full SPM.mat file containing the entire model. I ob-
%                      tained it by performing a full SPM analysis on non-
%                      smoothed volumes.
%outpath  [str]      : path to which the final .gii should be saved
%con_idx  [int]      : Index representing the contrast. 1=before button
%                      press vs. math, 2=after button press vs. math,
%                      3=before vs. after press, 4=average of before and   
%                      after press vs. math
% %hippunfold_path [str]     
%                    : path to the hippunfold output

%OPTIONAL:
%FWHM    [int]       : Smoothness Factor when using a Guassian average. 
%                      Default 0.9.
%idx     [cell]      : cell of size (hemis,subfields) containing the 
%                      indices of NON-vessel vertices. Output of 
%                      "VPF_layer_pipeline_hippocampus_EPI.m".The output 
%                      will contain only those vertices.

%OUTPUT:
%* .mat and .json files containing: T-map, contrast map, critial T-value at 
%  which p<0.05 FDR-corrected and maximum significant p-value (if not sig-
%  nificant: pmax = 1).
% 
%* sub-{subid}_hemi-avg_{surf}_T.shape.gii

%TODO:
%
%
% Viktor Pfaffenrot <viktor.pfaffenrot@uni-due.de>; January 2024, ELH

if nargin < 8 || isempty(ZTRANS)
    ZTRANS = false;
end

if nargin < 7 || isempty(idx)
    MASK_VEINS = false;
else
    MASK_VEINS = true;
end

if nargin < 6 || isempty(FWHM)
    SMOOTH = false;
else
    SMOOTH = true;
end


labels_dir = dir([hippunfold_path '/surf/sub-*_hemi-*_space-T2w_den-0p5mm_label-hipp_atlas-bigbrain_subfields.label.gii']);
subfields = size(layers,2);
vols = size(layers{1,1,1},3);
runs = size(layers,3);
Tvols = vols*runs;
N_layers = 3;

N_vertices = 7262;

W = SPM.xX.W;
GLM_use = SPM.xX.pKX;
con = SPM.xCon(con_idx).c;


%create layers with vertices: take inner and outer layer, average over
% hemispheres
layers_new = zeros([N_vertices,N_layers,vols, runs, 2]);
for hemi = 1:2
    labels = gifti([labels_dir(hemi).folder '/' labels_dir(hemi).name]).cdata;
    for run = 1:runs
        if ~MASK_VEINS
            for ii = 1:subfields-1
                inp = layers{hemi,ii,run}(:,[10 20 30],:);
                if SMOOTH
                    inp = smoothdata(inp,1,'gaussian','SmoothingFactor',FWHM);
                end
                layers_new(labels==ii,:,:,run,hemi) = inp;
            end
        else
            for ii = 1:subfields-1
                sz = size(layers{hemi,ii,run});
                sfield = zeros([sz(1) N_layers sz(end)]);
                inp = layers{hemi,ii,run}(idx{hemi,ii},[10 20 30],:);
                if SMOOTH
                    inp = smoothdata(inp,1,'gaussian','SmoothingFactor',FWHM);
                end
                sfield(idx{hemi,ii},:,:) = inp;
                layers_new(labels==ii,:,:,run,hemi) = sfield;
            end
        end
    end
end

layers = mean(layers_new,5);

clear layers_new

nonzeroidx = find(layers(:,1,1,1));
Y = reshape(layers(nonzeroidx,:,:,:),[numel(nonzeroidx)*N_layers vols runs]);

if ZTRANS
    %basline z-transform input. I take the math condition as baseline
    for run = 1:runs
        idx = find(SPM.xX.X(SPM.Sess(run).row,SPM.Sess(run).col(3))>0);
        m = mean(Y(:,idx,run),2);
        s = std(Y(:,idx,run),[],2);
        Y(:,:,run) = (Y(:,:,run)-m)./s;
    end
end

Y = reshape(Y,[size(Y,1) Tvols]);

KWY = spm_filter(SPM.xX.K,W*Y.');
b   = GLM_use*KWY;

res      = spm_sp('r',SPM.xX.xKXs,KWY);        %-Residuals
ResSS    = sum(res.^2);                    %-Residual SSQ
ResMS = ResSS / SPM.xX.trRV;
[tmpT,Tcrit,tmpcon,pmax] = VPF_Tmap_from_SPM(SPM,b,ResMS,con,0.05,'FDR');

tmpT = reshape(tmpT,[numel(nonzeroidx) N_layers]);
tmpcon= reshape(tmpcon,[numel(nonzeroidx) N_layers]);

T = zeros([N_vertices N_layers]);
con_array = T;
T(nonzeroidx,:) = tmpT;
con_array(nonzeroidx,:) = tmpcon;


results = struct('T',T,'con_array',con_array,'Tcrit',Tcrit,'pmax',pmax);
results_json = jsonencode(results,PrettyPrint=true);
if con_idx == 1
    outname = 'pre';
    results_pre = results;

elseif con_idx == 2
    outname = 'post';

    results_post = results;
elseif con_idx == 3
    outname = 'pre_vs_post';

    results_pre_vs_post = results;
else
    outname = 'memory_vs_math';
    results_memory_vs_math = results;
end

if ZTRANS
    outname_2 = [outname '_z'];
else
    outname_2 = outname;
end

if MASK_VEINS
    fulloutname = [outpath '/results_' outname_2 '_unfolded_vessel_masked.mat'];
else
    fulloutname = [outpath '/results_' outname_2 '_unfolded.mat'];
end
save(fulloutname,['results_' outname])
if MASK_VEINS
    fid = fopen([outpath '/results_' outname_2 '_unfolded_vessel_masked.json'],'w');
else
    fid = fopen([outpath '/results_' outname_2 '_unfolded.json'],'w');
end
fprintf(fid,'%s',results_json);
fclose(fid);


SURFS = {'inner','midthickness','outer'};

subfs_name = dir([hippunfold_path '/surf/sub-*_hemi-L_space-T2w_den-0p5mm_label-hipp_atlas-bigbrain_subfields.label.gii']);
mthick = dir([hippunfold_path '/surf/*_hemi-L_space-unfolded_den-0p5mm_label-hipp_midthickness.surf.gii']);

subid =  char(regexp(subfs_name.name,'sub-(\d+)_hemi','tokens','once'));
for SURF = 1:length(SURFS)
    inp = T(:,SURF);

    VPF_plot_hippocampus_unfolded(inp,[mthick.folder '/' mthick.name],[subfs_name.folder '/' subfs_name.name],'hot',[2 4],0);
    gaxis = gca;
    VPF_rot_hippocampus_flatmap(gaxis);
    pause(0.2)
    if MASK_VEINS
        title([SURFS{SURF} ' vessel masked']);
    else
        title(SURFS{SURF});
    end

    g = gifti(inp);
    if MASK_VEINS
        gii_outname = [outpath '/sub-' subid '_hemi-avg_' SURFS{SURF} '_T_' outname_2 '_vessel_masked.shape.gii'];
    else
        gii_outname = [outpath '/sub-' subid '_hemi-avg_' SURFS{SURF} '_T_' outname_2 '.shape.gii'];
    end
    save(g,gii_outname,'Base64Binary');    
end