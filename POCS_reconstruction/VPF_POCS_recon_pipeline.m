function [] =  VPF_POCS_recon_pipeline(list,PF_factor)
%Pipeline to perform a post-hoc POCS reconstruction on already coil-combined
%data. Uses some SPM12 functions. Calls POCS reconstruction written by 
%Michael Völker.
%Michael Völker (2024). MRI Partial Fourier reconstruction with POCS 
%(https://www.mathworks.com/matlabcentral/fileexchange/39350-mri-partial-fourier-reconstruction-with-pocs), 
%MATLAB Central File Exchange. Retrieved March 21, 2024.

%INPUT:
%list [struct]     : Structure containing the folder and name of 3D nifti images.
%                    Easiest way to obtain this is by first converting all ma-
%                    gnitude and phase dicoms to nifti using VPF_dcm2nii.m.
%                    This results in 3D nifti images for both mag. and phase.
%                    Then use the dir command and fetch all .nii files in your
%                    folder, e.g. list = dir([mypath '/*.nii']);
%
%PF_factor [double]: Partial Fourier factor. Typically 6/8.

%OUTPUT
%                    4D-nifti images mag_POCS.nii and pha_POCS.nii. The 
%                    latter might also be useful for NORDIC. Consider 
%                    compressing the niftis using e.g. pigz.

%Example usage:


N = size(list,1)/2;
    parfor ii = 1:N
        M = size(list,1)/2;
        mag_file = [list(ii).folder '/' list(ii).name];
        pha_file = [list(ii+M).folder '/' list(ii+M).name];
        VPF_POCS3D_nifti(mag_file,pha_file,PF_factor);
        delete(mag_file);
        delete(pha_file);
    end
    
    %%
    
    for ii = 1:N
        if ii == 1
            tmp = spm_vol([list(ii).folder '/' list(ii).name(1:end-4) '_POCS.nii']);
            mag = repmat(tmp,[N, 1]);
            mag(ii) = tmp;

            tmp = spm_vol([list(ii+N).folder '/' list(ii+N).name(1:end-4) '_POCS.nii']);
            pha = repmat(tmp,[N, 1]);
            pha(ii) = tmp;
        else
            mag(ii) = spm_vol([list(ii).folder '/' list(ii).name(1:end-4) '_POCS.nii']);
            pha(ii) = spm_vol([list(ii+N).folder '/' list(ii+N).name(1:end-4) '_POCS.nii']);
        end
    end
    
    spm_file_merge(mag,'mag_POCS.nii');
    spm_file_merge(pha,'pha_POCS.nii');
    %%
    delete(mag.fname)
    delete(pha.fname)
end