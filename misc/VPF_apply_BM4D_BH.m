clear;clc;

runs = 2;
echoes = 6;
load('/media/pfaffenrot/Elements/postdoc/projects/library/misc/BM4D_settings.mat');

subid = '7566';
dat_path = ['/media/pfaffenrot/Elements/postdoc/projects/hippocampus_breathhold/FLASH/derivatives/pipeline/' subid '/ses-02/func_bart/'];
for run = 1:runs
    for echo = 1:echoes
        fname = [dat_path 'run' num2str(run) '/mag/sub-' subid '_task_breathhold_run-0' num2str(run)...
                 '_part-mag_bold_echo-' num2str(echo) '.nii'];

        fname_out = [dat_path 'run' num2str(run) '/mag/sub-' subid '_task_breathhold_run-0' num2str(run)...
            '_part-mag_bold_echo-' num2str(echo) '.nii'];

        hdr = load_nifti(fname);
        
        img_denoised = zeros(size(hdr.vol));
        img = hdr.vol;
        parfor vol = 1:size(img_denoised,4)
            img_denoised(:,:,:,vol) = bm4d(img(:,:,:,vol),settings);
        end
        hdr.vol = img_denoised;
        
        save_nifti(hdr,fname);
        
    end
end