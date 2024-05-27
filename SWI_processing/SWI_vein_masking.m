clear;clc;

img_hdr = load_nifti('mSWI_denoised.nii',1);
img = load_nifti('mSWI_denoised.nii').vol;
m = load_nifti('mask.nii').vol;
sz = size(img);
kern = [2 2 1];
for ii = 1:sz(3)
    m(:,:,ii) = imdilate(imfill(imclose(m(:,:,ii),strel('disk',12)),'holes'),strel('disk',12));
end



V  = VPF_Frangi_Vesselness_Filter(img,m,1);
V(isnan(V)) = 0;
V = imgaussfilt3(V,[2.5 2.5 1]);
%%
ma = max(V(:)); 
mi = min(V(:));
V = 100 * (V - mi) / (ma - mi);
V(V<1) = 0;


img_hdr.vol = V.*m;
save_nifti(img_hdr,'SWI_vessel_masked.nii');