clear;clc;

pre_vs_post_flag = true;

N_subfields = 5;
N_subjects = 9;
N_layers = 30;

if pre_vs_post_flag == true
    load('/media/pfaffenrot/Elements/postdoc/projects/data/avg/memory/pre_vs_post_all_layers_masked.mat');
    load('/media/pfaffenrot/Elements/postdoc/projects/data/avg/memory/pre_vs_post_all_layers_no_masked.mat');



    masked = permute(reshape(pre_vs_post_all_layers_masked,[N_subfields N_subjects N_layers]),[1 3 2]);
    no_masked = permute(reshape(pre_vs_post_all_layers_no_masked,[N_subfields N_subjects N_layers]),[1 3 2]);
else
    load('/media/pfaffenrot/Elements/postdoc/projects/data/avg/memory/memory_vs_math_all_layers_masked.mat');
    load('/media/pfaffenrot/Elements/postdoc/projects/data/avg/memory/memory_vs_math_all_layers_no_masked.mat');



    masked = permute(reshape(memory_vs_math_all_layers_masked,[N_subfields N_subjects N_layers]),[1 3 2]);
    no_masked = permute(reshape(memory_vs_math_all_layers_no_masked,[N_subfields N_subjects N_layers]),[1 3 2]);
end
d = no_masked-masked;
mean_d = mean(d,3,'omitnan');
SE_d = std(d,[],3,'omitnan')./sqrt(N_subjects);
%%

[colorcode,bla] = VPF_create_hippocampus_colorcode();
mylegend = [cellstr('Subiculum'), cellstr('CA1'), cellstr('CA2'),...
    cellstr('CA3'),cellstr('CA4/DG')];

plotspecs = struct('FontName','Arial','FontSize',22,'colormap',jet(256),...
    'xtick',[],'xlim',[1 N_layers],'LineWidth',3);


plotspecs.ytick = -0.2:0.05:0.3;
plotspecs.ylim = [-0.2 0.3];

figure(),
for subfield = 1:N_subfields-1
    l = line([10,10],plotspecs.ylim,'LineWidth',plotspecs.LineWidth,'Color','black');
    plotspecs.color = colorcode{subfield,1};
    h = shadedErrorBar(1:N_layers,mean_d(subfield,:),SE_d(subfield,:),plotspecs,1);
    hold on
    if subfield == 1
        legendax = h(subfield).mainLine;
        hold on
    else
        legendax = cat(2,legendax,h.mainLine);
    end
end

plotspecs.Marker = '.';
plotspecs.MarkerSize = 50;
h = VPF_show(@plot,10,mean_d(end,10),[],[],[],[],plotspecs);
set(h, 'Color',colorcode{end,1});
h = errorbar(10,mean_d(end,10),SE_d(end,10),'LineWidth',plotspecs.LineWidth);
set(h, 'Color',colorcode{end,1});

h = VPF_show(@plot,N_layers,mean_d(end,end),[],[],[],[],plotspecs);
set(h, 'Color',colorcode{end,1});
h2 = errorbar(N_layers,mean_d(end,end),SE_d(end,end),'LineWidth',plotspecs.LineWidth);
set(h2, 'Color',colorcode{end,1});
legendax = cat(2,legendax,h);
g = gca;
set(g,'FontName',plotspecs.FontName,'FontSize',plotspecs.FontSize)
legend(legendax,mylegend,'Location','Northwest');
ylabel('no masking - masking [\Deltaz]')

if pre_vs_post_flag == true
    title('pre > post')
else
    title('memory > math')
end
