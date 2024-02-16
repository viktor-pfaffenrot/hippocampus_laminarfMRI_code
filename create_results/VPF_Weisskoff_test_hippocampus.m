clear;%clc;%close all
VOXELSPACE = 0;
PER_SUBFIELD = 1;
MASK_VEINS = 0;

data_hippunfold = dir('/media/pfaffenrot/Elements/postdoc/projects/data/7*');

subjects = {{'7512','7495'},{'7513','7566'},{'7551','7218'},{'7536','7474'},...
    {'7546','7498'},{'7568','7511'},{'7567','7553'},{'7560','7554'},{'7540','7491'}};

N_subjects = length(subjects);
hemis = {'L','R'};

rng(1,'twister');
w_test= cell(N_subjects,1);
for subject = 1:N_subjects
    globalpath = ['/media/pfaffenrot/Elements/postdoc/projects/hippocampus_VASO/derivatives/pipeline/'...
        subjects{subject}{1} '/ses-01/func'];
    data_base = [data_hippunfold(subject).folder '/' subjects{subject}{2}];

    %pick random run
%     run = randsample(3,1);
    run = 1;
    %calculate std only over the 'rest' volume, i.e. only over math condition
    load([globalpath '/rwls_stats_compcor_motion_confounds_w_WMmask/SPM.mat']);
    idx_vol = find(SPM.xX.X(SPM.Sess(run).row,SPM.Sess(run).col(3))>0);

    if VOXELSPACE == 1
        %create mask
        for hemi = 1:2
            hdr = load_nifti([data_base '/hippunfold/anat/sub-' subjects{subject}{2} '_hemi-' hemis{hemi} '_space-T2w_desc-subfields_atlas-bigbrain_dseg.nii.gz']);
            tmp = hdr.vol>0;
            if hemi == 1
                mask = tmp;
            end
            mask(tmp) = tmp(tmp);
        end

        currentpath = [globalpath '/run' num2str(run) '/func'];

        list = dir([currentpath '/*Warped-to-Anat.nii.gz']);
        for ii = 1:length(list)
            tmp = load_nifti([list(ii).folder '/' list(ii).name]).vol;
            tmp = tmp(mask);
            if ii == 1
                img = zeros([size(tmp,1) length(list)]);
            end
            img(:,ii) = tmp;
        end
    else %in vertex space
        if strcmp(subjects{subject}{2},'7495') || strcmp(subjects{subject}{2},'7566')
            load([data_base '/memory/ses-01/layers.mat']);
        else
            load([data_base '/memory/layers.mat']);
        end
        if ~PER_SUBFIELD
            vertices_2_use = 6000;
            N_vertices = 7262;
            N_layers = size(layers{1,1,1},2);
            N_vols = size(layers{1,1,1},3);
            labels_dir = dir([data_base '/hippunfold/surf/sub-*_hemi-*_space-T2w_den-0p5mm_label-hipp_atlas-bigbrain_subfields.label.gii']);
            img = zeros([N_vertices,N_layers,N_vols,2]);
            for hemi = 1:2
                labels = gifti([labels_dir(hemi).folder '/' labels_dir(hemi).name]).cdata;

                for ii = 1:5
                    sz = size(layers{hemi,ii,run});
                    img(labels==ii,:,:,hemi) = layers{hemi,ii,run};
                end
            end
            img = permute(img,[3 2 1 4]);
            img = permute(img(:,:,:),[3 1 2]);

            tmp = zeros(vertices_2_use,size(img,3));
            for layer = 1:N_layers
                parfor iter = 1:vertices_2_use
                    voxel = randsample(size(img,1),iter);
                    dat = img(voxel,idx_vol,layer);
                    dat(dat==0) = nan;
                    dat = mean(dat,1,'omitnan');
                    tmp(iter,layer) = squeeze(std(dat,[],2));
                end
            end
            w_test_result = tmp;
            w_test{subject} = tmp;
        else %PER_SUBFIELD
            N_subfields = size(layers,2);
            img = cell(N_subfields,1);
            if ~MASK_VEINS
                for subfield = 1:N_subfields
                    if subfield == N_subfields-1
                        img{subfield} = cat(1,mean(layers{1,subfield,run},2),mean(layers{2,subfield,run},2));
                    else
                        img{subfield} = cat(1,layers{1,subfield,run},layers{2,subfield,run});
                    end
                end
            else %MASK_VEINS
                for subfield = 1:N_subfields
                    if subfield == N_subfields-1
                        img{subfield} = cat(1,mean(layers{1,subfield,run}(idx{1,subfield},:,:),2),mean(layers{2,subfield,run}(idx{2,subfield},:,:),2));
                    else
                        img{subfield} = cat(1,layers{1,subfield,run}(idx{1,subfield},:,:),layers{2,subfield,run}(idx{2,subfield},:,:));
                    end
                end
            end

            clear layers

            img{5} = cat(1,img{5},img{6});
            img(6)= [];

            out_cell = cell(N_subfields-1,1);
            for subfield = 1:N_subfields-1
                fprintf(['subfield ' num2str(subfield) '\n']);
                sz = size(img{subfield});
                out = zeros(sz(1:2));
                if subfield == 2 || subfield == 3
                    layer_start = 1;
                else
                    layer_start = 10;
                end
                if subfield == 5
                    layer_start = 1;
                end
                for layer = layer_start:size(img{subfield},2)
                    tmp = squeeze(img{subfield}(:,layer,:));
                    parfor iter = 1:size(img{subfield},1)
                        voxel = randsample(size(tmp,1),iter);
                        dat = tmp(voxel,idx_vol);
%                         dat(dat==0) = nan;
                        dat = mean(dat,1);
                        out(iter,layer) = squeeze(std(dat,[],2));
                    end
                end
                out_cell{subfield} = out;
            end
            if subfield==N_subfields-1 && subject == 3
                keyboard
            end
                
            w_test_result = out_cell;
            w_test{subject} = out_cell;
            clear out_cell;
        end
    end

    out = 'Weisskoff_test_result';
    if PER_SUBFIELD==1
        out = [out '_per_subfield'];
    end

    if MASK_VEINS==1
        out = [out '_vessels_masked'];
    end


    if strcmp(subjects{subject}{2},'7495') || strcmp(subjects{subject}{2},'7566')
        save([data_base '/memory/ses-01/' out '.mat'],'w_test_result');
    else
        save([data_base '/memory/' out '.mat'],'w_test_result');
    end

end
save([data_hippunfold(1).folder '/Weisskoff_test_result_all_persubfield.mat'],'w_test')
