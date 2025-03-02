clear;clc;

datpath = '/media/pfaffenrot/PostDoc_data/projects/data';
subject_ids = [7218 7474 7491 7495 7498 7511 7553 7554 7566];

locations = {'anterior', 'posteriorbody', 'tail'};

out = zeros(2,7262,length(subject_ids));
idx = cell(2,5,length(locations));
for subject = 1:length(subject_ids)

    hippunfoldpath = [datpath '/' sprintf('%04.0f',subject_ids(subject)) '/hippunfold/surf'];
    
    % load the labels
    label_files = {[hippunfoldpath '/sub-' sprintf('%04.0f',subject_ids(subject)) '_hemi-L_space-T2w_den-0p5mm_label-hipp_atlas-bigbrain_subfields.label.gii']...
        [hippunfoldpath '/sub-' sprintf('%04.0f',subject_ids(subject)) '_hemi-R_space-T2w_den-0p5mm_label-hipp_atlas-bigbrain_subfields.label.gii']};

    subfield_labels = cell(2,1);
    for jj = 1:2 % [left/right]
        tmp = gifti(label_files{jj});
        subfield_labels{jj} = tmp.cdata;
    end


    if subject_ids(subject) ~= 7491
        load([datpath '/' sprintf('%04.0f',subject_ids(subject)) '/breathhold/dS_breathhold.mat']);
    end
    %load([datpath '/' sprintf('%04.0f',subject_ids(subject)) '/memory/layers.mat']);
    %load([datpath '/' sprintf('%04.0f',subject_ids(subject)) '/memory/SPM.mat']);
    

    % project the longitudinal mask onto a surface
    for location = 1:length(locations)
        inputpath = [datpath '/' sprintf('%04.0f',subject_ids(subject))...
                     '/hippunfold/anat/hippo_mask_' locations{location} '.nii.gz'];


        tmp = VPF_hippo2surf(inputpath,hippunfoldpath); % [surface hemisphere vertex]

        % fix interpolation errors
        tmp(tmp>=0.95) = 1;
        tmp(tmp<0.95) = 0;

        % take only inner surface
        tmp = squeeze(tmp(1,:,:));

        % get the index corresponding to the subfield label for each longitudinal label
        for subfield = 1:5
            for hemi = 1:2
                idx{hemi,subfield,location} = find(tmp(hemi,subfield_labels{hemi}==subfield)).';
            end
        end

        % average over hemispheres and vertices corresponding to the indices to get layer profiles
        % for each longitudinal direction

        % breath-hold experiment
        if subject_ids(subject) ~= 7491
            if location == 1
                tmp = VPF_create_layer_results_breathhold_longitudinal(dS_breathhold,idx(:,:,location));
                layers_longitudinal_breathhold = repmat(tmp,1,3);
            else
                layers_longitudinal_breathhold(:,location) = VPF_create_layer_results_breathhold_longitudinal(dS_breathhold,idx(:,:,location));
            end
        end

        % AM experiment
%         for con_idx = 3 % 3 = pre_vs_post, 4 = memory_vs_math
%             if location == 1
%                 tmp = VPF_create_layer_results_memory_longitudinal(layers,con_idx,idx(:,:,location),SPM,1);
%                 layers_longitudinal = repmat(tmp,1,3);
%             else
%                 layers_longitudinal(:,location) = VPF_create_layer_results_memory_longitudinal(layers,con_idx,idx(:,:,location),SPM,1);
%             end
%         end

        
    end

    %save([datpath '/' sprintf('%04.0f',subject_ids(subject)) '/memory/layers_longitudinal_pre_vs_post.mat'],'layers_longitudinal','idx','-v7.3')
    
    if subject_ids(subject) ~= 7491
      save([datpath '/' sprintf('%04.0f',subject_ids(subject)) '/breathhold/layers_longitudinal_breathhold.mat'],'layers_longitudinal_breathhold','-v7.3')
    end
end







