clear;clc;%close all;

data = dir('/media/pfaffenrot/PostDoc_data/projects/data/7*');
data = data([1 2 4:end]);


subjects = length(data);

N_vertices = 7262;
subjects_select = [1 2 3 4 5 6 7 8];

hemis = {'L','R'};
SURFS = {'inner','outer'};
%%
%average surfaces
outname_base = '/media/pfaffenrot/Elements/postdoc/projects/data/avg/hippunfold';

for SURF = 1:2
    for hemi = 1:2
        outname = [outname_base '/canonical_hemi-' hemis{hemi} '_' SURFS{SURF} '.native.surf.gii'];
        outname_std = [outname_base '/canonical_hemi-' hemis{hemi} '_' SURFS{SURF} '_stdev.native.shape.gii'];
        cmd = ['wb_command -surface-average ' outname ' -stddev ' outname_std];
        for subject = 1:subjects
            cmd2cat = [' -surf ' data(subject).folder '/' data(subject).name '/hippunfold/surf'...
                '/sub-' data(subject).name '_hemi-' hemis{hemi}...
                '_space-T2w_den-0p5mm_label-hipp_' SURFS{SURF} '.surf.gii '];
            cmd = cat(2,cmd,cmd2cat);
        end
        system(cmd);
    end
    outname = [outname_base '/canonical_hemi-avg_' SURFS{SURF} '.native.surf.gii'];
    cmd = ['wb_command -surface-average ' outname];
    for hemi = 1:2
        cmd2cat = [' -surf ' outname_base '/canonical_hemi-' hemis{hemi} '_' SURFS{SURF} '.native.surf.gii '];
        cmd = cat(2,cmd,cmd2cat);
    end
    system(cmd);
end

for SURF = 1:2
    for hemi = 1:2
        outname = [outname_base '/canonical_hemi-' hemis{hemi} '_' SURFS{SURF} '.unfolded.surf.gii'];
        outname_std = [outname_base '/canonical_hemi-' hemis{hemi} '_' SURFS{SURF} '_stdev.unfolded.shape.gii'];
        cmd = ['wb_command -surface-average ' outname ' -stddev ' outname_std];
        for subject = 1:subjects
            cmd2cat = [' -surf ' data(subject).folder '/' data(subject).name '/hippunfold/surf'...
                '/sub-' data(subject).name '_hemi-' hemis{hemi}...
                '_space-unfolded_den-0p5mm_label-hipp_' SURFS{SURF} '.surf.gii '];
            cmd = cat(2,cmd,cmd2cat);
        end
        system(cmd);
    end
    outname = [outname_base '/canonical_hemi-avg_' SURFS{SURF} '.unfolded.surf.gii'];
    cmd = ['wb_command -surface-average ' outname];
    for hemi = 1:2
        cmd2cat = [' -surf ' outname_base '/canonical_hemi-' hemis{hemi} '_' SURFS{SURF} '.unfolded.surf.gii '];
        cmd = cat(2,cmd,cmd2cat);
    end
    system(cmd);
end

%% average angiograms, T2s and breathhold signal change

subfs_name = dir([data(1).folder '/' data(1).name '/hippunfold/surf/sub-*_hemi-L_space-T2w_den-0p5mm_label-hipp_atlas-bigbrain_subfields.label.gii']);
% subfs_name = dir('/media/pfaffenrot/Elements/postdoc/projects/data/avg/hippunfold/sub-7495_hemi-avg_space-T2w_den-0p5mm_label-hipp_atlas-bigbrain_subfields.label.gii');
mthick = dir('/media/pfaffenrot/PostDoc_data/projects/data/avg/hippunfold/canonical_hemi-avg_inner.unfolded.surf.gii');

outname_base = '/media/pfaffenrot/Elements/postdoc/projects/data/avg/breathhold';
% attribute = [{'angio'},{'T2s'},{'dSbreathhold'},{'dSbreathhold_weighted'}];
attribute = {'angio'};
titles = [{'vessels'}, {'T_2^*'},{'\DeltaS_{mean}'},{'\DeltaS_{mean,weighted}'}];
mycolor = [{'fireice'},{'inferno'},{'fireice'},{'fireice'}];
myscale = [[-10 10];[0 70];[-90 90];[-15 15]];


for ii = 1:length(attribute)
    for SURF = 1:2
        tmp_surface = zeros([N_vertices,subjects]);
        for subject = subjects_select
            g = gifti([data(subject).folder '/' data(subject).name '/breathhold/sub-'...
                data(subject).name '_hemi-avg_' SURFS{SURF} '_' attribute{ii} '.shape.gii']);
            tmp_surface(:,subject) = g.cdata;
        end
        avg_surface = mean(tmp_surface,2);

        if strcmp(attribute{ii},'angio')
            std_surface = zeros(size(tmp_surface,1),1);
            std_surface(avg_surface<0) = std(tmp_surface(avg_surface<0,:),[],2)*(-1);
            std_surface(avg_surface>0) = std(tmp_surface(avg_surface>0,:),[],2);
        else
            std_surface = std(tmp_surface,[],2);
        end

        VPF_plot_hippocampus_unfolded(avg_surface,[mthick.folder '/' mthick.name],[subfs_name.folder '/' subfs_name.name],mycolor{ii},myscale(ii,:),1);
        if ~strcmp(attribute{ii},'angio')
            colorbar;
        end
        gaxis = gca;
        VPF_rot_hippocampus_flatmap(gaxis);
        title([titles{ii} ' ' SURFS{SURF}])
        pause(0.2)

        if strcmp(attribute{ii},'angio')
            VPF_plot_hippocampus_unfolded(std_surface,[mthick.folder '/' mthick.name],[subfs_name.folder '/' subfs_name.name],mycolor{ii},[-10 10],1);
            if ~strcmp(attribute{ii},'angio')
                colorbar;
            end
            gaxis = gca;
            VPF_rot_hippocampus_flatmap(gaxis);
            title([titles{ii} ' ' SURFS{SURF} ' SD'])
            pause(0.2)       
        end

        avg_surface = gifti(avg_surface);
        save(avg_surface,[outname_base '/canonical_hemi-avg_' SURFS{SURF} '_' attribute{ii} '.shape.gii'],'Base64Binary');
        std_surface = gifti(std_surface);
        save(std_surface,[outname_base '/canonical_hemi-avg_' SURFS{SURF} '_stdev_' attribute{ii} '.shape.gii'],'Base64Binary');
    end
end
%% average memory beta-maps and tSNR
clear;clc;%close all;

dat_dir = '/media/pfaffenrot/Elements/postdoc/projects/data';
data = dir([dat_dir '/7*']);
subjects = length(data);

N_vertices = 7262;
subjects_select = [1 2 3 4 5 6 7 8];

hemis = {'L','R'};
SURFS = {'inner','outer'};


subfs_name = dir([data(1).folder '/' data(1).name '/hippunfold/surf/sub-*_hemi-L_space-T2w_den-0p5mm_label-hipp_atlas-bigbrain_subfields.label.gii']);
mthick = dir('/media/pfaffenrot/Elements/postdoc/projects/data/avg/hippunfold/canonical_hemi-avg_inner.unfolded.surf.gii');

outname_base = '/media/pfaffenrot/Elements/postdoc/projects/data/avg/memory';

attribute = {'tSNR','tSNR_vessel_masked','beta_map','beta_map_vessel_masked'};
titles = {'tSNR','tSNR vessel masked','beta map','beta map vessel masked'};

for ii = 1:length(attribute)
    for SURF = 1:2
        avg_surface = zeros([N_vertices,subjects]);
        for subject = subjects_select
            if strcmp(data(subject).name,'7495') || strcmp(data(subject).name,'7566') %subjects with 2 sessions
                g = zeros(N_vertices,1);
                for session = 1:2
                    tmp = gifti([data(subject).folder '/' data(subject).name '/memory/ses-0' num2str(session) '/native/sub-'...
                        data(subject).name '_hemi-avg_' SURFS{SURF} '_' attribute{ii} '.shape.gii']);
                    g = g + tmp.cdata;
                end
                g = g./2;
                avg_surface(:,subject) = g;
            else
                g = gifti([data(subject).folder '/' data(subject).name '/memory/native/sub-'...
                    data(subject).name '_hemi-avg_' SURFS{SURF} '_' attribute{ii} '.shape.gii']);
                avg_surface(:,subject) = g.cdata;
            end

        end
        avg_surface(avg_surface==0) = nan;
        std_surface = std(avg_surface,[],2,'omitnan');
        avg_surface = mean(avg_surface,2,'omitnan');

        VPF_plot_hippocampus_unfolded(avg_surface,[mthick.folder '/' mthick.name],[subfs_name.folder '/' subfs_name.name],'fireice',[0 40],1);
        colorbar;
        gaxis = gca;
        VPF_rot_hippocampus_flatmap(gaxis);
        title([titles{ii} ' ' SURFS{SURF}])
        pause(0.2)
        avg_surface = gifti(avg_surface);
        save(avg_surface,[outname_base '/canonical_hemi-avg_' SURFS{SURF} '_' attribute{ii} '.shape.gii'],'Base64Binary');
        std_surface = gifti(std_surface);
        save(std_surface,[outname_base '/canonical_hemi-avg_' SURFS{SURF} '_stdev_' attribute{ii} '.shape.gii'],'Base64Binary');
    end
end