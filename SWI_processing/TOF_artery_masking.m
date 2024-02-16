clear;clc;

img = load_nifti('TOF_denoised.nii').vol;
m = load_nifti('mask.nii').vol;
sz = size(img);
img2 = zeros([2*sz(1:2) sz(3)]);
for ii = 1:sz(3)
    img2(:,:,ii) = imresize(img(:,:,ii),2*sz(1:2),'cubic');
end
[a,h,v,d] = haart2(img2,1);

b = (h+v+d);
b = b./max(b(:));

c = a.*b.*b;
c(c<0.1) = 0;
c(c>1) = 1;


c = logical(c);
for ii = 1:sz(3)
    c(:,:,ii) = imclose(c(:,:,ii), strel('disk', 2));
    m(:,:,ii) = imdilate(imfill(imclose(m(:,:,ii),strel('disk',12)),'holes'),strel('disk',12)); 
end

c = c.*m;
img2 = img.*c;

ma = max(img2(:)); 
mi = min(img2(:));
img2 = 100 * (img2 - mi) / (ma - mi);
img2(img2<1) = 0;




hdr = load_nifti('TOF_denoised.nii',1);
hdr.vol = img2;
save_nifti(hdr,'TOF_vessel_masked.nii');

%%

reslice_flags = struct('mask',0,'mean',0,'interp',4,'which',1);
mprage = '/media/pfaffenrot/Elements/postdoc/projects/hippocampus_breathhold/FLASH/derivatives/pipeline/7553/ses-01/anat/T1/presurf_MPRAGEise/presurf_UNI/sub-7553_UNI_MPRAGEised_biascorrected_denoised.nii';

tmp = table2array(readtable('/media/pfaffenrot/Elements/postdoc/projects/hippocampus_breathhold/FLASH/derivatives/pipeline/7553/ses-01/anat/SWI/M_itk.txt'));
M = ea_antsmat2mat(str2num(tmp{1,2})',str2num(tmp{2,2})');
tmp = spm_vol('TOF_vessel_masked.nii'); 
tmp_orient = spm_get_space(tmp.fname);
spm_get_space(tmp.fname,M*tmp_orient);  

spm_reslice(cellstr(char(mprage,['TOF_vessel_masked.nii'])),reslice_flags);

