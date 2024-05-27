%% VPF_layer_pipeline_hippocampus_EPI
%% load data
clear;clc;close all

% some paramters
subid = '7560';
subid_BH = '7554';
sessionid = '1';

vols = 486;
runs = 3;
Tvols = runs*vols;

N = [20 10]; %number of layers

% this is the overall path
globalpath = ['/media/pfaffenrot/Elements/postdoc/projects/hippocampus_VASO/derivatives/pipeline/' subid '/ses-0' sessionid '/'];

% this is the path where the actual data are saved.
currentpath = 'func';

% this is the path where the freesurfer output is stored
structpath = ['/media/pfaffenrot/Elements/postdoc/projects/hippocampus_breathhold/FLASH/derivatives/pipeline/' subid_BH '/ses-01/anat/T1/presurf_MPRAGEise/presurf_UNI/'];



% this is the path where all the layer-related data (sampled functionals etc.) are saved
layerpath = [globalpath 'layerfication'];

%--> read structural data
mprage = [structpath 'sub-' subid_BH '_UNI_MPRAGEised_biascorrected_denoised.nii'];
mprage_mask = [structpath 'sub-' subid_BH '_UNI_MPRAGEised_brainmask.nii'];
hdr_mprage = spm_vol(mprage);
dims       = hdr_mprage.dim;
%<--

hippunfold_path = ['/media/pfaffenrot/Elements/postdoc/projects/hippocampus_breathhold/FLASH/derivatives/hippunfold/' subid_BH];
surf_path = [hippunfold_path '/surf'];

SPM_path = [globalpath 'func/rwls_stats_compcor_motion_confounds_w_WMmask/SPM.mat'];
load(SPM_path);

SWI_file = ['/media/pfaffenrot/Elements/postdoc/projects/hippocampus_breathhold/FLASH/derivatives/pipeline/' subid_BH '/ses-01/anat/SWI/rSWI_vessel_masked.nii'];



%% layer sampling
rule = '1 2 3 4 5 6';
subfields = numel(double(str2num(rule)));
if subfields == 0
    subfields = 1;
end

for run = 1:runs
    list = dir([globalpath '/func/run' num2str(run) '/func/mag_*Warped-to-Anat.nii.gz']);
    sampled_img = cell(size(list));

    for ii = 1:length(list)
        sampled_img{ii} = [list(ii).folder '/' list(ii).name];
    end

    if run == 1
        [tmp,idx] = VPF_create_hippocampus_layers(sampled_img,str2double(subid_BH),hippunfold_path,mprage,N,rule,SWI_file);
        layers = cell([size(tmp),runs]);
        layers(:,:,run) = tmp;
    else
        [layers(:,:,run),idx] = VPF_create_hippocampus_layers(sampled_img,str2double(subid_BH),hippunfold_path,mprage,N,rule,SWI_file);
    end
end


for run = 1:runs
    for subfield = 1:subfields
        for hemi = 1:2
            sz = size(layers{hemi,subfield,run});
            tmp = reshape(layers{hemi,subfield,run},[prod(sz(1:2)) sz(3)]);
            tmp =  tmp -(SPM.xX.K(run).X0*(SPM.xX.K(run).X0.'*tmp.')).';
            layers{hemi,subfield,run} = reshape(tmp,sz);
        end
    end
end

save([layerpathbase '/layers.mat'],'layers','idx','-v7.3')

%%
colorcode = VPF_create_hippocampus_colorcode();
for ZTRANS = [true,false]

    if ZTRANS
        layerpath = [layerpathbase '/z_transformed'];
    else
        layerpath = [layerpathbase '/native'];
    end

    for con_idx = 1:4 %1=pre_vs_math,2=post_vs_math,3=pre_vs_post,4_memory_vs_math
        for uu = 1:2 %1 == no vessels masked, 2 == vessels masked

            if uu == 1
                tmp = VPF_create_layer_results_memory(layers,SPM,con_idx,[],ZTRANS);
                results = repmat(tmp,[2 1]);
            else
                results(uu) = VPF_create_layer_results_memory(layers,SPM,con_idx,idx,ZTRANS);
            end

            titles = [cellstr('Subiculum'), cellstr('CA1'), cellstr('CA2'),...
                cellstr('CA3'),cellstr('CA4/DG')];

            plotspecs = struct('FontName','Arial','FontSize',22,'colormap',jet(256),...
                'xtick',[],'xlim',[1 30],'LineWidth',3);
            plotspecs.color = colorcode(:,1);

            if con_idx == 1
                figname = 'pre vs math';
            elseif con_idx == 2
                figname = 'post vs math';
            elseif con_idx == 3
                figname = 'pre vs post';
            else
                figname = 'memory vs math';
            end

            if uu == 2
                if ZTRANS
                    figname_out = [figname '_z_vessels_masked'];
                else
                    figname_out = [figname '_vessels_masked'];
                end
            else
                if ZTRANS
                    figname_out = [figname '_z'];
                else
                    figname_out = figname;
                end
            end


            f = figure('Name',figname_out,'NumberTitle','off');
            fields_to_plot=length(titles);
            if ZTRANS
                plotspecs.ytick = -0.4:0.2:2;
                plotspecs.ylim = [-0.4 2];
                mytitle = '\beta_{memory} - \beta_{math} [z]';
            else
                plotspecs.ytick = -10:5:30;
                plotspecs.ylim = [-10 30];
                mytitle = '\beta_{memory} - \beta_{math} [a. u.]';
            end
            h = VPF_show(@plot,1:30,results(uu).con_array(1:fields_to_plot,:).',[],[],[],mytitle,plotspecs);
            for kk = 1:length(h)
                set(h(kk), 'Color',plotspecs.color{kk});
            end
            l = line([10,10],plotspecs.ylim,'LineWidth',plotspecs.LineWidth,'Color','black');
            legend(titles);

            results_json = jsonencode(results(uu),PrettyPrint=true);
            if con_idx == 1
                outname = 'results_pre';
                results_pre = results(uu);

            elseif con_idx == 2
                outname = 'results_post';

                results_post = results(uu);
            elseif con_idx == 3
                outname = 'results_pre_vs_post';

                results_pre_vs_post = results(uu);
            else
                outname = 'results_memory_vs_math';
                results_memory_vs_math = results(uu);
            end

            if uu == 2
                if ZTRANS
                    fulloutname = [layerpath '/' outname '_full_tSNR_z_vessel_masked.mat'];
                else
                    fulloutname = [layerpath '/' outname '_full_tSNR_vessel_masked.mat'];
                end
            else
                if ZTRANS
                    fulloutname = [layerpath '/' outname '_full_tSNR_z.mat']; 
                else
                    fulloutname = [layerpath '/' outname '_full_tSNR.mat']; 
                end
            end
            save(fulloutname,outname)
            fid = fopen([fulloutname(1:end-4) '.json'],'w');
            fprintf(fid,'%s',results_json);
            fclose(fid);
        end

        for ii = 1:2
            if ii == 1
                idx_inp = [];
            else
                idx_inp = idx;
            end

            if con_idx == 3 || con_idx == 4
                VPF_create_surface_Tmap_results_memory(layers,...
                    SPM,layerpath,con_idx,hippunfold_path,0.6,idx_inp,ZTRANS);
            end

            if con_idx == 3 && ~ZTRANS
                VPF_create_surface_tSNR_results_memory(layers,...
                    layerpath,hippunfold_path,SPM,idx_inp);
                pause(0.1)
            end
        end
    end
end




