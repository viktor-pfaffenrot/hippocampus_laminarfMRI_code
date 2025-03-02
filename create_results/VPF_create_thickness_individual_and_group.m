clear;clc;close all

data = dir('/media/pfaffenrot/PostDoc_data/projects/data/7*');
subjects = length(data);

subfields = 5;
n_vertices = 7262;

subids = cell(length(data),1);
thicknesses = zeros(subfields,subjects,2);
for subject = 1:subjects
    subids{subject} = data(subject).name;

    thickness = ft_read_cifti([data(subject).folder '/' data(subject).name '/' ...
        'hippunfold/surf/sub-' data(subject).name '_space-T2w_den-0p5mm_label-hipp_thickness.dscalar.nii']);

    labels_L = gifti([data(subject).folder '/' data(subject).name '/' ...
        'hippunfold/surf/sub-' data(subject).name '_hemi-L_space-T2w_den-0p5mm_label-hipp_atlas-bigbrain_subfields.label.gii']);

    labels_R = gifti([data(subject).folder '/' data(subject).name '/' ...
        'hippunfold/surf/sub-' data(subject).name '_hemi-R_space-T2w_den-0p5mm_label-hipp_atlas-bigbrain_subfields.label.gii']);

    thickness_per_subfield = zeros(subfields,1);
    thickness_per_subfield_with_SRLM = zeros(subfields,1);

    thickness_left = thickness.dscalar(1:n_vertices);
    thickness_right = thickness.dscalar(n_vertices+1:end);
    for subfield = 1:subfields
        
        thickness_per_subfield(subfield) = (median(thickness_left(labels_L.cdata==subfield)) + ...
                                            median(thickness_left(labels_R.cdata==subfield)))./2;

        thickness_per_subfield_with_SRLM(subfield) = thickness_per_subfield(subfield);

        if (subfield == 2 || subfield == 3)
            thickness_per_subfield_with_SRLM(subfield) = thickness_per_subfield_with_SRLM(subfield)+...
                thickness_per_subfield_with_SRLM(subfield)/2;
        end
    end

    outname = [data(subject).folder '/' data(subject).name '/' ...
        'hippunfold/thickness_per_sufield.mat'];
    save(outname,'thickness_per_subfield');
    
    thicknesses(:,subject,1) = thickness_per_subfield;
    thicknesses(:,subject,2) = thickness_per_subfield_with_SRLM;
end

avg_thickness = squeeze(mean(thicknesses,2));
std_thickness = squeeze(std(thicknesses,[],2));

save('/media/pfaffenrot/PostDoc_data/projects/data/avg/hippunfold/avg_thickness_per_subfield.mat','avg_thickness','std_thickness');