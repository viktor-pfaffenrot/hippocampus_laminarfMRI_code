function VPF_sample_india_ink_stain(experiment_mean,india_ink_path)

if nargin < 2
    india_ink_path = '/home/pfaffenrot/work/my_publics/hippocampus_layers/images/F_8/india_ink_sampling';
end

subfield_label = {'Subiculum','CA1','CA2','CA3','CA4/DG'};
N_main = 30;
profiles = nan(length(subfield_label),N_main);
for subfield = 1:length(subfield_label)
    if subfield == 2 || subfield == 3
        start = 1;
        N = N_main;
    else
        start = 10;
        N = N_main - start;
    end

    if subfield == 5
        subfield_folder = 'CA4';
    else
        subfield_folder = subfield_label{subfield};
    end

    vasc = spm_read_vols(spm_vol([india_ink_path '/' subfield_folder '/' subfield_folder '.nii']));

    rim = spm_read_vols(spm_vol([india_ink_path '/' subfield_folder '/rim_layers_equidist.nii']));
    rim(rim==0)= nan;




    profile = zeros(N,1);
    for ii = 1:N
        profile(ii) = mean(vasc(rim==ii),'omitnan');
    end
    profile = profile(end:-1:1);
    profile = smoothdata(profile,"gaussian",10);
    profile = 2*mean(profile,'omitnan')-profile;

    profile = profile./mean(profile)*experiment_mean(subfield);
    if ~(subfield == 2 || subfield == 3)
        profile = interp1(linspace(1,20,20).',profile,linspace(1,21,21).','linear','extrap');
    end
    
    profiles(subfield,start:end) = profile;
end
%%
colorcode = VPF_create_hippocampus_colorcode();

plotspecs = struct('FontName','Arial','FontSize',22,'colormap',jet(256),...
    'xtick',[],'xlim',[1 30],'LineWidth',2);
plotspecs.color = colorcode(:,1);

plotspecs.ylim = [0.9 1.4];
plotspecs.ytick = 0.9:0.1:1.4;

figure,
h = VPF_show(@plot,1:30,profiles,[],'India Ink Stain',[],[],plotspecs);
for kk = 1:length(h)
    set(h(kk), 'Color',plotspecs.color{kk});
end
hold on
plotspecs.Marker = '.';
plotspecs.MarkerSize = 30;
h = VPF_show(@plot,10,profiles(end,10),[],[],[],[],plotspecs);
set(h, 'Color',plotspecs.color{end});
h = VPF_show(@plot,30,profiles(end,end),[],'India Ink Stain',[],[],plotspecs);
set(h, 'Color',plotspecs.color{end});
hold off
l = line([10,10],plotspecs.ylim,'LineWidth',plotspecs.LineWidth,'Color','black');
legend(subfield_label,'Location','northwest');
end