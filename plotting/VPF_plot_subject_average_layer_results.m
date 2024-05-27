clear;clc;
% breat-hold signal change
load('/media/pfaffenrot/Elements/postdoc/projects/data/avg/breathhold/results_avg.mat')
colorcode = VPF_create_hippocampus_colorcode();

N_subjects = 8;
N_subfields = 5;


results_avg.dS_mean_echo_substd = results_avg.dS_mean_echo_substd./sqrt(N_subjects);
results_avg.dS_mean_echo_weighted_substd = results_avg.dS_mean_echo_weighted_substd./sqrt(N_subjects);
mylegend = [cellstr('Subiculum'), cellstr('CA1'), cellstr('CA2'),...
    cellstr('CA3'),cellstr('CA4/DG')];

plotspecs = struct('FontName','Arial','FontSize',22,'colormap',jet(256),...
    'xtick',[],'xlim',[1 30],'LineWidth',4);

plotspecs.ytick = -4:2:24;
plotspecs.ylim = [-4 24];

clear h
N_layers = 30;
figure
for kk = 1:N_subfields-1
    plotspecs.color = colorcode{kk,1};
    h(kk) = shadedErrorBar(1:N_layers,results_avg.dS_mean_echo_submean(:,kk),...
        results_avg.dS_mean_echo_substd(:,kk),plotspecs,1);
    if kk == 1
        legendax = h(kk).mainLine;
        hold on
    else
        legendax = cat(2,legendax,h(kk).mainLine);
    end
end

line([10,10],plotspecs.ylim,'LineWidth',plotspecs.LineWidth-1,'Color','black');

hold on
plotspecs.color = colorcode{end,1};
plotspecs.Marker = '.';
plotspecs.MarkerSize = 42;
h = VPF_show(@plot,10,results_avg.dS_mean_echo_submean(10,end),[],[],[],[],plotspecs);
set(h, 'Color',plotspecs.color);
h = errorbar(10,results_avg.dS_mean_echo_submean(10,end),results_avg.dS_mean_echo_substd(10,end),'LineWidth',plotspecs.LineWidth);
set(h, 'Color',plotspecs.color);

h = VPF_show(@plot,N_layers,results_avg.dS_mean_echo_submean(end,end),[],[],[],[],plotspecs);
set(h, 'Color',plotspecs.color);
h2 = errorbar(N_layers,results_avg.dS_mean_echo_submean(end,end),results_avg.dS_mean_echo_substd(end,end),'LineWidth',plotspecs.LineWidth);
set(h2, 'Color',plotspecs.color);
hold off

legendax = cat(2,legendax,h);

legend(legendax,mylegend,'Location','Northwest');
set(gca,'FontSize',plotspecs.FontSize);
ylabel('\DeltaS_{mean} [a.u.]')


plotspecs = struct('FontName','Arial','FontSize',22,'colormap',jet(256),...
    'xtick',[],'xlim',[1 30],'LineWidth',4);
clear h
plotspecs.ytick = 0:0.2:3;
plotspecs.ylim = [0 3];

N_layers = 30;
figure
plotspecs.Marker = 'None';
for kk = 1:N_subfields-1
    plotspecs.color = colorcode{kk,1};
    h(kk) = shadedErrorBar(1:N_layers,results_avg.dS_mean_echo_weighted_submean(:,kk),...
        results_avg.dS_mean_echo_weighted_substd(:,kk),plotspecs,1);
    if kk == 1
        legendax = h(kk).mainLine;
        hold on
    else
        legendax = cat(2,legendax,h(kk).mainLine);
    end
end

line([10,10],plotspecs.ylim,'LineWidth',plotspecs.LineWidth,'Color','black');

hold on
plotspecs.color = colorcode{end,1};
plotspecs.Marker = '.';
plotspecs.MarkerSize = 42;
h = VPF_show(@plot,10,results_avg.dS_mean_echo_weighted_submean(10,end),[],[],[],[],plotspecs);
set(h, 'Color',plotspecs.color);
h = errorbar(10,results_avg.dS_mean_echo_weighted_submean(10,end),results_avg.dS_mean_echo_weighted_substd(10,end),'LineWidth',plotspecs.LineWidth);
set(h, 'Color',plotspecs.color);

h = VPF_show(@plot,N_layers,results_avg.dS_mean_echo_weighted_submean(end,end),[],[],[],[],plotspecs);
set(h, 'Color',plotspecs.color);
h2 = errorbar(N_layers,results_avg.dS_mean_echo_weighted_submean(end,end),results_avg.dS_mean_echo_weighted_substd(end,end),'LineWidth',plotspecs.LineWidth);
set(h2, 'Color',plotspecs.color);
hold off

legendax = cat(2,legendax,h);

legend(legendax,mylegend,'Location','Northwest');
ylabel('\DeltaS_{mean,weighted} [a.u.]');
set(gca,'FontSize',plotspecs.FontSize);

%% fitted T2*
load('/media/pfaffenrot/Elements/postdoc/projects/data/avg/breathhold/results_avg.mat')
colorcode = VPF_create_hippocampus_colorcode();
N_layers = 30;
N_subfields = 5;
N_subjects = 8;
results_avg.T2s_substd = results_avg.T2s_substd./sqrt(N_subjects);

mylegend = [cellstr('Subiculum'), cellstr('CA1'), cellstr('CA2'),...
    cellstr('CA3'),cellstr('CA4/DG')];

plotspecs = struct('FontName','Arial','FontSize',22,'colormap',jet(256),...
    'xtick',[],'xlim',[1 30],'LineWidth',4);
plotspecs.ytick = 22:2:42;
plotspecs.ylim = [22 42];

clear h
figure
for kk = 1:N_subfields-1
    plotspecs.color = colorcode{kk,1};
    h(kk) = shadedErrorBar(1:N_layers,results_avg.T2s_submean(:,kk),...
        results_avg.T2s_substd(:,kk),plotspecs,1);
    if kk == 1
        legendax = h(kk).mainLine;
        hold on
    else
        legendax = cat(2,legendax,h(kk).mainLine);
    end
end

line([10,10],plotspecs.ylim,'LineWidth',plotspecs.LineWidth-1,'Color','black');

hold on
plotspecs.color = colorcode{end,1};
plotspecs.Marker = '.';
plotspecs.MarkerSize = 42;
h = VPF_show(@plot,10,results_avg.T2s_submean(10,end),[],[],[],[],plotspecs);
set(h, 'Color',plotspecs.color);
h = errorbar(10,results_avg.T2s_submean(10,end),results_avg.T2s_substd(10,end),'LineWidth',plotspecs.LineWidth);
set(h, 'Color',plotspecs.color);

h = VPF_show(@plot,N_layers,results_avg.T2s_submean(end,end),[],[],[],[],plotspecs);
set(h, 'Color',plotspecs.color);
h2 = errorbar(N_layers,results_avg.T2s_submean(end,end),results_avg.T2s_substd(end,end),'LineWidth',plotspecs.LineWidth);
set(h2, 'Color',plotspecs.color);
hold off

legendax = cat(2,legendax,h);

legend(legendax,mylegend,'Location','Northwest');
set(gca,'FontSize',plotspecs.FontSize);
ylabel('T_2^* [ms]')







%% layer profiles of contrast differences
N_subjects = 9;
N_subfields = 5;
N_layers = 30;
load('/media/pfaffenrot/Elements/postdoc/projects/data/avg/memory/memory_vs_math_z_all_layers_masked.mat')
memory_vs_math = reshape(memory_vs_math,[N_subfields N_subjects N_layers]);
Xm = squeeze(mean(memory_vs_math,2)).';
Xs = squeeze(std(memory_vs_math,[],2)).'./sqrt(N_subjects);

colorcode = VPF_create_hippocampus_colorcode();
mylegend = [cellstr('Subiculum'), cellstr('CA1'), cellstr('CA2'),...
    cellstr('CA3'),cellstr('CA4/DG')];


plotspecs = struct('FontName','Arial','FontSize',22,'colormap',jet(256),...
    'xtick',[],'xlim',[1 30],'LineWidth',4);
plotspecs.ytick = -0.4:0.1:4;
plotspecs.ylim = [-0.4 1.4];

plotspecs.Marker = 'None';
clear h
figure,
for kk = 1:N_subfields-1
    plotspecs.color = colorcode{kk,1};
    if kk ~= N_subfields
        h(kk) = shadedErrorBar(1:N_layers,Xm(:,kk),...
            Xs(:,kk),plotspecs,1);
        if kk == 1
            legendax = h(kk).mainLine;
            hold on
        else
            legendax = cat(2,legendax,h(kk).mainLine);
        end
    end
end
line([10,10],plotspecs.ylim,'LineWidth',plotspecs.LineWidth,'Color','black');
plotspecs.color = colorcode{end,1};
plotspecs.Marker = '.';
plotspecs.MarkerSize = 42;
h = VPF_show(@plot,10,Xm(10,end),[],[],[],[],plotspecs);
set(h, 'Color',plotspecs.color);
h = errorbar(10,Xm(10,end),Xs(10,end),'LineWidth',plotspecs.LineWidth);
set(h, 'Color',plotspecs.color);

h = VPF_show(@plot,N_layers,Xm(end,end),[],[],[],[],plotspecs);
set(h, 'Color',plotspecs.color);
h2 = errorbar(N_layers,Xm(end,end),Xs(end,end),'LineWidth',plotspecs.LineWidth);
set(h2, 'Color',plotspecs.color);
hold off
legendax = cat(2,legendax,h);




ylabel('\beta_{memory} - \beta_{math} [z]');
legend(legendax,mylegend,'Location','Northwest');
ColorOrder = get(gca,'ColorOrder');
set(gca,'FontSize',plotspecs.FontSize);
%% tSNR
N_subjects = 9;
N_subfields = 5;
N_layers = 30;
load('/media/pfaffenrot/Elements/postdoc/projects/data/avg/memory/memory_vs_math_math_tSNR_all_layers_masked.mat')
colorcode = VPF_create_hippocampus_colorcode();
mylegend = [cellstr('Subiculum'), cellstr('CA1'), cellstr('CA2'),...
    cellstr('CA3'),cellstr('CA4/DG')];


plotspecs = struct('FontName','Arial','FontSize',22,'colormap',jet(256),...
    'xtick',[],'xlim',[1 30],'LineWidth',4);


plotspecs.ytick = 40:20:180;
plotspecs.ylim = [40 180];

plotspecs.Marker = 'None';
clear h
figure,
for kk = 1:N_subfields-1
    plotspecs.color = colorcode{kk,1};
    if kk ~= N_subfields
        h(kk) = shadedErrorBar(1:N_layers,m_tSNR(:,kk),...
            s_tSNR(:,kk),plotspecs,1);
        if kk == 1
            legendax = h(kk).mainLine;
            hold on
        else
            legendax = cat(2,legendax,h(kk).mainLine);
        end
    end
end
line([10,10],plotspecs.ylim,'LineWidth',plotspecs.LineWidth,'Color','black');
plotspecs.color = colorcode{end,1};
plotspecs.Marker = '.';
plotspecs.MarkerSize = 42;
h = VPF_show(@plot,10,m_tSNR(10,end),[],[],[],[],plotspecs);
set(h, 'Color',plotspecs.color);
h = errorbar(10,m_tSNR(10,end),s_tSNR(10,end),'LineWidth',plotspecs.LineWidth);
set(h, 'Color',plotspecs.color);

h = VPF_show(@plot,N_layers,m_tSNR(end,end),[],[],[],[],plotspecs);
set(h, 'Color',plotspecs.color);
h2 = errorbar(N_layers,m_tSNR(end,end),s_tSNR(end,end),'LineWidth',plotspecs.LineWidth);
set(h2, 'Color',plotspecs.color);
hold off
legendax = cat(2,legendax,h);




ylabel('tSNR [a. u.]');
legend(legendax,mylegend,'Location','Northwest');
ColorOrder = get(gca,'ColorOrder');
set(gca,'FontSize',plotspecs.FontSize);
