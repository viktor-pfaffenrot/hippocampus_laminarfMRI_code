clear;clc;
load('/media/pfaffenrot/Elements/postdoc/projects/data/avg/breathhold/results_avg.mat')
[colorcode,colorcode_points] = VPF_create_hippocampus_colorcode();

select = [1 2 3 4 5 6 7 8];
N = length(select);

mylegend = [cellstr('Subiculum'), cellstr('CA1'), cellstr('CA2'),...
    cellstr('CA3'),cellstr('CA4/DG')];

plotspecs = struct('FontName','Arial','FontSize',22,'colormap',jet(256),...
    'xtick',[],'xlim',[1 30],'LineWidth',2);

plotspecs.ytick = 0:2:20;
plotspecs.ylim = [0 20];

pos = 30;
figure
for kk = 1:5
    plotspecs.color = colorcode{kk,1};
    h(kk) = shadedErrorBar(1:pos,results_avg.dS_mean_echo_submean(:,kk),...
        results_avg.dS_mean_echo_substd(:,kk)./N,plotspecs,1);
    if kk == 1
        legendax = h(kk).mainLine;
        hold on
    else
        legendax = cat(2,legendax,h(kk).mainLine);
    end
    if kk == 5
        hold on
        plotspecs.Marker = '.';
        plotspecs.MarkerSize = 30;
        h = VPF_show(@plot,10,results_avg.dS_mean_echo_submean(10,kk),[],[],[],[],plotspecs);
        set(h, 'Color',plotspecs.color);
        h = VPF_show(@plot,pos,results_avg.dS_mean_echo_submean(end,kk),[],[],[],[],plotspecs);
        set(h, 'Color',plotspecs.color);
        hold off
    end
end
line([10,10],plotspecs.ylim,'LineWidth',plotspecs.LineWidth,'Color','black');
ylabel('\DeltaS_{mean} [a.u.]')
    legend(legendax,mylegend,'Location','Northwest');
set(gca,'FontSize',plotspecs.FontSize);

clear h
plotspecs.ytick = 0.4:0.2:2.6;
plotspecs.ylim = [0.4 2.6];

pos = 30;
figure
plotspecs.Marker = 'None';
for kk = 1:5
    plotspecs.color = colorcode{kk,1};
    h(kk) = shadedErrorBar(1:pos,results_avg.dS_mean_echo_weighted_submean(:,kk),...
        results_avg.dS_mean_echo_weighted_substd(:,kk)./N,plotspecs,1);
    if kk == 1
        legendax = h(kk).mainLine;
        hold on
    else
        legendax = cat(2,legendax,h(kk).mainLine);
    end
    if kk == 5
        hold on
        plotspecs.Marker = '.';
        plotspecs.MarkerSize = 30;
        h = VPF_show(@plot,10,results_avg.dS_mean_echo_weighted_submean(10,kk),[],[],[],[],plotspecs);
        set(h, 'Color',plotspecs.color);
        h = VPF_show(@plot,pos,results_avg.dS_mean_echo_weighted_submean(end,kk),[],[],[],[],plotspecs);
        set(h, 'Color',plotspecs.color);
        hold off
    end
end
line([10,10],plotspecs.ylim,'LineWidth',plotspecs.LineWidth,'Color','black');
ylabel('\DeltaS_{mean,weighted} [a.u.]');
    legend(legendax,mylegend,'Location','Northwest');
ColorOrder = get(gca,'ColorOrder');
set(gca,'FontSize',plotspecs.FontSize);

%%
clear h
plotspecs.ytick = 22:2:40;
plotspecs.ylim = [22 40];

pos = 30;
figure
plotspecs.Marker = 'None';
for kk = 1:5
    plotspecs.color = colorcode{kk,1};
    h(kk) = shadedErrorBar(1:pos,results_avg.T2s_submean(:,kk),...
        results_avg.T2s_substd(:,kk)./N,plotspecs,1);
    if kk == 1
        legendax = h(kk).mainLine;
        hold on
    else
        legendax = cat(2,legendax,h(kk).mainLine);
    end
    if kk == 5
        hold on
        plotspecs.Marker = '.';
        plotspecs.MarkerSize = 30;
        h = VPF_show(@plot,10,results_avg.T2s_submean(10,kk),[],[],[],[],plotspecs);
        set(h, 'Color',plotspecs.color);
        h = VPF_show(@plot,pos,results_avg.T2s_submean(end,kk),[],[],[],[],plotspecs);
        set(h, 'Color',plotspecs.color);
        hold off
    end
end
line([10,10],plotspecs.ylim,'LineWidth',plotspecs.LineWidth,'Color','black');
ylabel('T_2^* [ms]');
    legend(legendax,mylegend,'Location','Northwest');
ColorOrder = get(gca,'ColorOrder');
set(gca,'FontSize',plotspecs.FontSize);
