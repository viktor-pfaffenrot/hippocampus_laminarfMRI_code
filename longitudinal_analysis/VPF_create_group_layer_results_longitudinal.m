clear;clc;
data = dir('/media/pfaffenrot/PostDoc_data/projects/data/7*');
locations = {'anterior', 'posteriorbody', 'tail'};
N_locations = length(locations);
subjects = length(data);


%% breath-hold task
for subject = 1:subjects
    if ~strcmp(data(subject).name,'7491')
        load([data(subject).folder '/' data(subject).name...
            '/breathhold/layers_longitudinal_breathhold.mat']);
        if subject == 1
            X = zeros([size(layers_longitudinal_breathhold(1).con_array),N_locations,subjects]);
        end
        for location = 1:N_locations
            X(:,:,location,subject) = layers_longitudinal_breathhold(location).con_array;
        end
    end
end

X(:,:,:,3) = [];


outpath = '/media/pfaffenrot/PostDoc_data/projects/data/avg/breathhold';
outname = 'breathhold_longitudinal';
breathhold_longitudinal = reshape(permute(X,[1 4 2 3]),[40 30 3]);

save([outpath '/' outname '.mat'],'breathhold_longitudinal')


%% AM task

for subject = 1:subjects

    if strcmp(data(subject).name,'7495') ||  strcmp(data(subject).name,'7566')
        for session = 1:2
            load([data(subject).folder '/' data(subject).name '/memory/ses-0' num2str(session)...
                '/layers_longitudinal_pre_vs_post.mat'],'layers_longitudinal');
            if session == 1
                tmp = zeros([size(layers_longitudinal(1).con_array) N_locations]);
                if subject == 1
                    X = zeros([size(layers_longitudinal(1).con_array),N_locations,subjects]);
                end
            end
            for location = 1:N_locations
                tmp(:,:,location) = tmp(:,:,location)+layers_longitudinal(location).con_array;
                X(:,:,location,subject) = tmp(:,:,location)./2;
            end
        end
    else
        load([data(subject).folder '/' data(subject).name...
            '/memory/layers_longitudinal_pre_vs_post.mat'],'layers_longitudinal');
        if subject == 1
            X = zeros([size(layers_longitudinal(1).con_array),N_locations,subjects]);
        end
        for location = 1:N_locations
            X(:,:,location,subject) = layers_longitudinal(location).con_array;
        end
    end
end


outpath = '/media/pfaffenrot/PostDoc_data/projects/data/avg/memory';
outname = 'pre_vs_post_z_all_layers_no_masked_longitudinal';
pre_vs_post_longitudinal = reshape(permute(X,[1 4 2 3]),[45 30 3]);

save([outpath '/' outname '.mat'],'pre_vs_post_longitudinal')