clear;clc;

subj = '7218';
thres_SWI = 3;
thres_TOF = 1;
dir = ['/media/pfaffenrot/Elements/postdoc/projects/hippocampus_breathhold/FLASH/derivatives/hippunfold/' subj '/surf'];
% dir = ['/media/pfaffenrot/Elements/postdoc/projects/hippocampus_VASO/derivatives//hippunfold/' subj '/surf'];
datdir = ['/media/pfaffenrot/Elements/postdoc/projects/hippocampus_breathhold/FLASH/derivatives/pipeline/' subj '/ses-01/anat'];
% datdir = ['/media/pfaffenrot/Elements/postdoc/projects/hippocampus_VASO/derivatives/pipeline/' subj '/ses-01/anat'];

places = [{'inner'},{'outer'}];

for place = places  

    hemis = ['L','R'];

    count = 1;
    for hemi = hemis
        dat = gifti([datdir '/SWI/SWI_to_' place{1} '_' hemi '.shape.gii']);
        dat_TOF = gifti([datdir '/TOF/TOF_to_' place{1} '_' hemi '.shape.gii']);
        subfs = [dir '/sub-' subj '_hemi-' hemi '_space-T2w_den-0p5mm_label-hipp_atlas-bigbrain_subfields.label.gii'];

        if strcmp(hemi,hemis(1))
            dat_to_plot = zeros([2,size(dat.cdata)]);
        end
        %threshold masked SWI image
        dat_to_plot(count,dat.cdata>thres_SWI)= dat.cdata(dat.cdata>thres_SWI);
        %convert scale of SWI to be between 0 and -100
        dat_to_plot(count,dat.cdata>thres_SWI) = dat_to_plot(count,dat.cdata>thres_SWI)./max(dat_to_plot(count,dat.cdata>thres_SWI))*100*(-1);

        %threshold masked TOF image
        dat_to_plot(count,dat_TOF.cdata>thres_TOF)= dat_TOF.cdata(dat_TOF.cdata>thres_TOF);
        %convert scale of TOF to be between 0 and 100
        dat_to_plot(count,dat_TOF.cdata>thres_TOF) = dat_to_plot(count,dat_TOF.cdata>thres_TOF)./max(dat_to_plot(count,dat_TOF.cdata>thres_TOF))*100;

        count = count +1;
    end

    dat_to_plot = squeeze(mean(dat_to_plot,1));
    g = gifti(dat_to_plot.');
    save(g,[datdir '/sub-' subj '_hemi-avg_' place{:} '_angio.shape.gii'],'Base64Binary');    

    mthick = [dir '/sub-' subj '_hemi-L_space-unfolded_den-0p5mm_label-hipp_midthickness.surf.gii'];
    
    VPF_plot_hippocampus_unfolded(dat_to_plot,mthick,subfs,'fireice',[-30 30],false);
    gaxis = gca;
    VPF_rot_hippocampus_flatmap(gaxis);



    %calculate the density of big veins and arteries per subfield. I do this by
    % counting vertices, i.e. I take the number of all negative vertices for veins
    % and of all positive vertices for arteries for each subfield and then divide
    % by the total number of vertices in each subfield
    subfs = gifti(subfs);
    subfields =  unique(subfs.cdata);

    if strcmp(place,places(1))
        vessel_density = zeros(2,length(subfields),2); %[veins/arteries, subfield, inner/outer]
        idx = 1;
    else
        idx = 2;
    end
    for subfield = 1:length(subfields)
        denom = length(subfs.cdata(subfs.cdata==subfield));

        vessels = dat_to_plot(subfs.cdata==subfield);
        vein_nominator = length(vessels(vessels<0));

        vessel_density(1,subfield,idx) = vein_nominator/denom*1000;
        vessel_density(2,subfield,idx) = length(vessels(vessels>0))/denom*1000;
    end
end

save([datdir '/vessel_density.mat'],'vessel_density')

