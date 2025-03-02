clear;clc;
% breat-hold signal change
N_subjects = 8;
N_subfields = 5;
N_layers = 30;
N_locations = 3;

load('/media/pfaffenrot/PostDoc_data/projects/data/avg/breathhold/breathhold_longitudinal.mat')
breathhold_longitudinal = reshape(breathhold_longitudinal,[N_subfields N_subjects N_layers N_locations]);
Xm = permute(squeeze(mean(breathhold_longitudinal,2)), [2 1 3]);
Xs = permute(squeeze(std(breathhold_longitudinal,[],2)), [2 1 3])./sqrt(N_subjects);

colorcode = VPF_create_hippocampus_colorcode();
mytitles = [cellstr('Subiculum'), cellstr('CA1'), cellstr('CA2'),...
    cellstr('CA3'),cellstr('CA4/DG')];

mylegend = [cellstr('Anterior'), cellstr('Posterior body'), cellstr('Tail')];

plotspecs = struct('FontName','Arial','FontSize',28,'colormap',jet(256),...
    'xtick',[],'xlim',[1 30],'LineWidth',4,'LineStyle','-');
plotspecs.ytick = 0:0.4:6;
plotspecs.ylim = [0 6];

plotspecs.Marker = 'None';

my_LineStyle = [cellstr('-'),cellstr('--'),cellstr('-.')];


for ss = 1:N_subfields-1
    clear h
    figure,
    for kk = 1:N_locations
        plotspecs.color = colorcode{ss,1};
        plotspecs.LineStyle = my_LineStyle{kk};
        if kk ~= N_subfields
            h(kk) = shadedErrorBar(1:N_layers,Xm(:,ss,kk),...
                Xs(:,ss,kk),plotspecs,1);
            if kk == 1
                legendax = h(kk).mainLine;
                hold on
            else
                legendax = cat(2,legendax,h(kk).mainLine);
            end
        end
    end
    line([10,10],plotspecs.ylim,'LineWidth',plotspecs.LineWidth,'Color','black');

    ylabel('\DeltaS_{weighted} [a.u.]');
    title(mytitles{ss})
    legend(legendax,mylegend,'Location','Northwest');
    set(gca,'FontSize',plotspecs.FontSize);
end


my_Marker = [cellstr('.'),cellstr('square'),cellstr('diamond')];
my_MarkerSize = [42,18,18];
clear h
figure,
plotspecs.color = colorcode{end,1};
for kk = 1:N_locations
    plotspecs.Marker = my_Marker{kk};
    plotspecs.MarkerSize = my_MarkerSize(kk);
    h = VPF_show(@plot,10,Xm(10,end,kk),[],[],[],[],plotspecs);
    if kk == 1
        legendax = h;
        hold on
    else
        legendax = cat(2,legendax,h);
    end
    set(h, 'Color',plotspecs.color);
    h = errorbar(10,Xm(10,end,kk),Xs(10,end,kk),'LineWidth',plotspecs.LineWidth);
    set(h, 'Color',plotspecs.color);


    h = VPF_show(@plot,N_layers,Xm(end,end,kk),[],[],[],[],plotspecs);
    set(h, 'Color',plotspecs.color);
    h2 = errorbar(N_layers,Xm(end,end,kk),Xs(end,end,kk),'LineWidth',plotspecs.LineWidth);
    set(h2, 'Color',plotspecs.color);
end
    ylabel('\DeltaS_{weighted} [a.u.]');
title(mytitles{end})
legend(legendax,mylegend,'Location','Northwest');
ColorOrder = get(gca,'ColorOrder');
set(gca,'FontSize',plotspecs.FontSize);



%% layer profiles of contrast differences
N_subjects = 9;
N_subfields = 5;
N_layers = 30;
N_locations = 3;

load('/media/pfaffenrot/PostDoc_data/projects/data/avg/memory/pre_vs_post_z_all_layers_no_masked_longitudinal.mat')
pre_vs_post_longitudinal = reshape(pre_vs_post_longitudinal,[N_subfields N_subjects N_layers N_locations]);
Xm = permute(squeeze(mean(pre_vs_post_longitudinal,2)), [2 1 3]);
Xs = permute(squeeze(std(pre_vs_post_longitudinal,[],2)), [2 1 3])./sqrt(N_subjects);

colorcode = VPF_create_hippocampus_colorcode();
mytitles = [cellstr('Subiculum'), cellstr('CA1'), cellstr('CA2'),...
    cellstr('CA3'),cellstr('CA4/DG')];

mylegend = [cellstr('Anterior'), cellstr('Posterior body'), cellstr('Tail')];

plotspecs = struct('FontName','Arial','FontSize',28,'colormap',jet(256),...
    'xtick',[],'xlim',[1 30],'LineWidth',4,'LineStyle','-');
plotspecs.ytick = -0.5:0.5:2;
plotspecs.ylim = [-0.5 2];

plotspecs.Marker = 'None';

my_LineStyle = [cellstr('-'),cellstr('--'),cellstr('-.')];


for ss = 1:N_subfields-1
    clear h
    figure,
    for kk = 1:N_locations
        plotspecs.color = colorcode{ss,1};
        plotspecs.LineStyle = my_LineStyle{kk};
        if kk ~= N_subfields
            h(kk) = shadedErrorBar(1:N_layers,Xm(:,ss,kk),...
                Xs(:,ss,kk),plotspecs,1);
            if kk == 1
                legendax = h(kk).mainLine;
                hold on
            else
                legendax = cat(2,legendax,h(kk).mainLine);
            end
        end
    end
    line([10,10],plotspecs.ylim,'LineWidth',plotspecs.LineWidth,'Color','black');

    ylabel('\beta_{pre} - \beta_{post} [z]');
    title(mytitles{ss})
    legend(legendax,mylegend,'Location','Northwest');
    set(gca,'FontSize',plotspecs.FontSize);
end


my_Marker = [cellstr('.'),cellstr('square'),cellstr('diamond')];
my_MarkerSize = [42,18,18];
clear h
figure,
plotspecs.color = colorcode{end,1};
for kk = 1:N_locations
    plotspecs.Marker = my_Marker{kk};
    plotspecs.MarkerSize = my_MarkerSize(kk);
    h = VPF_show(@plot,10,Xm(10,end,kk),[],[],[],[],plotspecs);
    if kk == 1
        legendax = h;
        hold on
    else
        legendax = cat(2,legendax,h);
    end
    set(h, 'Color',plotspecs.color);
    h = errorbar(10,Xm(10,end,kk),Xs(10,end,kk),'LineWidth',plotspecs.LineWidth);
    set(h, 'Color',plotspecs.color);


    h = VPF_show(@plot,N_layers,Xm(end,end,kk),[],[],[],[],plotspecs);
    set(h, 'Color',plotspecs.color);
    h2 = errorbar(N_layers,Xm(end,end,kk),Xs(end,end,kk),'LineWidth',plotspecs.LineWidth);
    set(h2, 'Color',plotspecs.color);
end
ylabel('\beta_{pre} - \beta_{post} [z]');
title(mytitles{end})
legend(legendax,mylegend,'Location','Northwest');
ColorOrder = get(gca,'ColorOrder');
set(gca,'FontSize',plotspecs.FontSize);



