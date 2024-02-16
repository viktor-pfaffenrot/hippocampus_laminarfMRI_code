function VPF_create_surface_tSNR_results_memory(layers,outpath,hippunfold_path,idx)
% This function calculates tSNR and projects it onto the inner and outer
% surface of the unfolded hippocampus (averaged over hemispheres) created 
% by hippunfold. The surface is saved as .gii.

%INPUT:
%layers   [cell]     : cell of size (hemis,subfields,runs) containing the 
%                      sampled data (size = (N_vertices,N_layers,N_volumes).
%                      ASSUMES DG TO BE IN THERE, I.E. SUBFIELDS = 6!
%                      Output of "VPF_layer_pipeline_hippocampus_EPI.m"
%outpath   [str]     : path to which the final .gii should be saved
%hippunfold_path [str]     
%                    : path to the hippunfold output

%OPTIONAL:
%idx     [cell]      : cell of size (hemis,subfields) containing the 
%                      indices of NON-vessel vertices. Output of 
%                      "VPF_layer_pipeline_hippocampus_EPI.m".The output 
%                      will contain only those vertices.

%OUTPUT:
% sub-{subid}_hemi-avg_{surf}_tSNR.shape.gii


%TODO:
%
%

% Viktor Pfaffenrot <viktor.pfaffenrot@uni-due.de>; January 2024, ELH


if nargin < 4 || isempty(idx)
    MASK_VEINS = false;
else
    MASK_VEINS = true;
end


labels_dir = dir([hippunfold_path '/surf/sub-*_hemi-*_space-T2w_den-0p5mm_label-hipp_atlas-bigbrain_subfields.label.gii']);
subfields = size(layers,2);
vols = size(layers{1,1,1},3);
runs = size(layers,3);
N_layers = 2;

N_vertices = 7262;

%create layers with vertices: take inner and outer layer, average over
% hemispheres
layers_new = zeros([N_vertices,N_layers,vols, runs, 2]);
for hemi = 1:2
    labels = gifti([labels_dir(hemi).folder '/' labels_dir(hemi).name]).cdata;
    for run = 1:runs
        if ~MASK_VEINS
            for ii = 1:subfields-1
                inp = layers{hemi,ii,run}(:,[10 30],:);

                layers_new(labels==ii,:,:,run,hemi) = inp;
            end
        else
            for ii = 1:subfields-1
                sz = size(layers{hemi,ii,run});
                sfield = zeros([sz(1) 2 sz(end)]);
                inp = layers{hemi,ii,run}(idx{hemi,ii},[10 30],:);

                sfield(idx{hemi,ii},:,:) = inp;
                layers_new(labels==ii,:,:,run,hemi) = sfield;
            end
        end
    end
end

layers = mean(layers_new,5);

clear layers_new

%average over runs and take the 1st 100 volumes
layers = mean(layers(:,:,1:100,:),4);

tSNR = mean(layers,3)./std(layers,[],3);
tSNR(isnan(tSNR)) = 0;



SURFS = {'inner','outer'};

subfs_name = dir([hippunfold_path '/surf/sub-*_hemi-L_space-T2w_den-0p5mm_label-hipp_atlas-bigbrain_subfields.label.gii']);
mthick = dir([hippunfold_path '/surf/*_hemi-L_space-unfolded_den-0p5mm_label-hipp_midthickness.surf.gii']);

subid =  char(regexp(subfs_name.name,'sub-(\d+)_hemi','tokens','once'));

for SURF = 1:length(SURFS)
    inp = tSNR(:,SURF);

    VPF_plot_hippocampus_unfolded(inp,[mthick.folder '/' mthick.name],[subfs_name.folder '/' subfs_name.name],'hot',[0 60],0);
    gaxis = gca;
    VPF_rot_hippocampus_flatmap(gaxis);
    pause(0.2)
    if MASK_VEINS
        title(['tSNR ' SURFS{SURF} ' vessel masked']);
    else
        title(['tSNR ' SURFS{SURF}]);
    end

    g = gifti(inp);
    
    if MASK_VEINS
        gii_outname = [outpath '/sub-' subid '_hemi-avg_' SURFS{SURF} '_tSNR_vessel_masked.shape.gii'];
    else
        gii_outname = [outpath '/sub-' subid '_hemi-avg_' SURFS{SURF} '_tSNR.shape.gii'];
    end

    save(g,gii_outname,'Base64Binary');
end
end