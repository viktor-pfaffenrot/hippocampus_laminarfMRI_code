function out = VPF_hippo2surf(inputpath,hippunfold_path,space)

if nargin < 3
    space = 'T2w';
end
hemis = {'L','R'};
SURFS = {'inner','outer'};

[folder,name,ext] = fileparts(inputpath);

[~,name,ext2] = fileparts(name);
ext = [ext2 ext];

surfdex = 1;
for SURF = 1:length(SURFS)
    hemidex = 1;
    for hemi = 1:length(hemis)
        outname = [folder '/' name '_to_' SURFS{SURF} '_' hemis{hemi} '_' space '.shape.gii'];
        try
            dat = gifti(outname);
        catch

            cmd = ['wb_command -volume-to-surface-mapping ' inputpath ' '...
                hippunfold_path '/*hemi-' hemis{hemi} '_space-' space '_den-0p5mm_label-hipp_' SURFS{SURF} '.surf.gii '...
                outname ' -cubic'
                ];

            system(cmd);
            dat = gifti(outname);
        end

        if strcmp(hemis{hemi},hemis{1}) && strcmp(SURFS{SURF},SURFS{1})
            out = zeros([2,2,size(dat.cdata)]);
        end
        out(surfdex,hemidex,:)= dat.cdata;
        hemidex = hemidex + 1;
    end
    surfdex = surfdex + 1;
end
end