clear;clc;

load('/media/pfaffenrot/Elements/postdoc/projects/hippocampus_breathhold/BM4D_settings.mat');
settings.parameters.denoise.sigma = 40;


types = [{'TOF'},{'SWI'}];

for type = types
    inp = ['/media/pfaffenrot/Elements/postdoc/projects/hippocampus_breathhold/FLASH'...
        '/derivatives/pipeline/7553/ses-01/anat/' type{1} '/' type{1} '.nii'];


    img = load_nifti(inp).vol;

    img_denoised = bm4d(img,settings);

    hdr = load_nifti(inp,1);
    hdr.vol = img_denoised;

    [filepath,filename,ext] = fileparts(inp);

    save_nifti(hdr,[filepath '/' filename '_denoised' ext]);
end