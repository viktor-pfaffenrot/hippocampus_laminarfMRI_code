clear;clc;
%create roi of bad regions due to heartbeat and breathing based on the residuals of a first-pass GLM

sub_id = '7568';
sess_id ='01'; 
% mainpath = '/home/pfaffenrot/work/postdoc/projects/ANT_workdir/';
mainpath = '/media/pfaffenrot/Elements/postdoc/projects/hippocampus_VASO/derivatives/pipeline/';
statspath = [mainpath sub_id '/ses-' sess_id '/func/rwls_stats/'];
% statspath = '/home/pfaffenrot/work/postdoc/projects/ANT_workdir/rwls_stats/';
structpath =  '/media/pfaffenrot/Elements/postdoc/projects/hippocampus_breathhold/FLASH/derivatives/pipeline/7511/ses-01/anat/T1/presurf_MPRAGEise/presurf_UNI/sub-7511_UNI_MPRAGEised_biascorrected.nii';
%structpath =  '/media/pfaffenrot/Elements/postdoc/projects/hippocampus_VASO/derivatives/pipeline/7513/ses-01/anat/T1/presurf_MPRAGEise/presurf_UNI/sub-7513_UNI_MPRAGEised_biascorrected_denoised.nii';
mask_path = [mainpath sub_id '/ses-' sess_id '/func/run1/func/mag_POCS_r1_fixedMask.nii'];
unfoldpath = ['/media/pfaffenrot/Elements/postdoc/proje' ...
    'cts/hippocampus_breathhold/FLASH/derivatives/hippunfold/7511/anat/']; 
%unfoldpath = '/media/pfaffenrot/Elements/postdoc/projects/hippocampus_VASO/derivatives/hippunfold/7513/anat/'; 

%% 
%create hippocampus mask

hippo_mask_hdr = load_nifti([unfoldpath 'sub-7511_hemi-L_space-T2w_desc-subfields_atlas-bigbrain_dseg.nii.gz']);
hippo_mask_L = hippo_mask_hdr.vol>0;
hippo_mask_R = load_nifti([unfoldpath 'sub-7511_hemi-R_space-T2w_desc-subfields_atlas-bigbrain_dseg.nii.gz']).vol>0;

hippo_mask_L(hippo_mask_R>0) = hippo_mask_R(hippo_mask_R>0);
hippo_mask = imdilate(logical(hippo_mask_L),strel('disk',2));
hippo_mask_hdr.vol = hippo_mask;

save_nifti(hippo_mask_hdr,[unfoldpath 'hippo_mask.nii']);
%%
hippo_mask_path = [unfoldpath 'hippo_mask.nii'];

hdr = load_nifti([statspath 'ResMS.nii']);
ResMS = hdr.vol;

% ResMS(isnan(ResMS)) = 0;


roi = abs(ResMS)>3*std(ResMS(:),[],1,"omitnan");
hdr.vol = roi;

save_nifti(hdr,[statspath 'high_ResMS_mask.nii']);

roi_WM = load_nifti([unfoldpath 'WM_mask_midbrain.nii.gz']).vol;


%%
%bring the fixedmask of the ANTs layer pipeline into anatomical space and clean it up a little bit
transpath = '/home/pfaffenrot/work/postdoc/projects/ANT_workdir/custom_reg2anat_run1.txt';

reslice_flags = struct('mask',0,'mean',0,'interp',4,'which',[2 0]);

trans = table2array(readtable(transpath));
trans = ea_antsmat2mat(str2num(trans{1,2})',str2num(trans{2,2})');

spm_get_space(mask_path,trans*spm_get_space(mask_path));
spm_reslice(cellstr(char(structpath,mask_path)),reslice_flags);

%%
[p,m,ext] = fileparts(mask_path); 
mask = load_nifti([p '/r' m ext]).vol;
% mask = load_nifti(mask_path).vol;
mask(abs(mask)<1e-3) = 0;
mask(isnan(mask)) = 0;
mask = logical(mask);

for ii = 1:size(mask,3)
    mask(:,:,ii) = imdilate(imfill(imclose(mask(:,:,ii),strel('disk',12)),'holes'),strel('disk',12)); 
end
%%
% load the data into memory. To save space, mask and vectorize themmask = load_nifti([p '/r' m ext]).vol;

runs = 3;
N    = 486;

load([statspath 'SPM.mat']);

%some voxels with high residuals might be within the hippocampus. It is essential to not incorporate them into acompcor.
%Hence, use a hippocampus mask based on hippounfold segmentations, dilate it a bit and mask out those hippocampus voxels.
roi(hippo_mask) = 0;
roi = roi(mask);

roi_WM(hippo_mask) = 0;
roi_WM = roi_WM(mask);
%%
for run = 1:runs
%     mypath = ['/media/pfaffenrot/Elements/postdoc/projects/hippocampus_VASO/derivatives/pipeline/7512/ses-02/func/run' num2str(run) '/func/'];
mypath = [mainpath sub_id '/ses-' sess_id '/func/run' num2str(run) '/func/'];
% mypath = [mainpath 'mag_POCS_r' num2str(run) '_split/'];

    for vol = 1:N

        if vol == 1
            hdr = load_nifti([mypath 'mag_POCS_r' num2str(run) '_1' sprintf('%03d',vol-1) '_Warped-to-Anat.nii.gz']);
            tmp = hdr.vol(mask);
            img = zeros(length(tmp),N);
            img(:,vol) = tmp;
            clear tmp
        else
            tmp = load_nifti([mypath 'mag_POCS_r' num2str(run) '_1' sprintf('%03d',vol-1) '_Warped-to-Anat.nii.gz']).vol;
            img(:,vol) = tmp(mask);
            if vol == 10
                keyboard
            end
        end
    end
    

%     confounds = SPM.xX.pKX([(1:8)+8*(run-1) 25+(run-1)],1+(run-1)*N:run*N);
    
    confounds = SPM.xX.pKX((4:9)+9*(run-1),1+(run-1)*N:run*N);


    X = fmri_compcor(img,{roi, roi_WM},[0.5 5],'confounds',confounds');

    writematrix(X,[mypath 'compcor_motion_confounds_w_WM_mask.txt'],'Delimiter',' ')
end
