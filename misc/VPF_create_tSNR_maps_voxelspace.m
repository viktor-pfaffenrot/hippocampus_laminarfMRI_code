clear;clc;
main_path = '/media/pfaffenrot/Elements/postdoc/projects/hippocampus_VASO/derivatives/pipeline/7513/ses-02/func';
SPM_path = [main_path '/rwls_stats_compcor_motion_confounds_w_WMmask/SPM.mat'];
load(SPM_path);

math_vols = find(SPM.xX.X(SPM.Sess(1).row,SPM.Sess(1).col(3))>0);
math_vols = math_vols(1:100);


data = dir([main_path '/run1/func/*Warped-to-Anat.nii.gz']);

for vol = 1:length(math_vols)
    if vol == 1
        hdr = load_nifti([data(math_vols(vol)).folder '/' data(math_vols(vol)).name]);
        img = zeros([size(hdr.vol), 100]);
    end
    img(:,:,:,vol) = load_nifti([data(math_vols(vol)).folder '/' data(math_vols(vol)).name]).vol;
end



tSNR = mean(img,4)./std(img,[],4);


hdr.vol = tSNR;
save_nifti(hdr,'tSNR_ses02.nii.gz');
%%