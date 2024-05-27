clear;clc;

%LME model similar to what is used for each individual contrast, but now 
%for the difference. 

globalpath = '/media/pfaffenrot/Elements/postdoc/projects/data/avg/memory';
load([globalpath '/pre_vs_post_all_layers_masked.mat']);
load([globalpath '/memory_vs_math_all_layers_masked.mat']);


N_subjects = 9;
N_subfields = 5;
N_layers = 30;
subfield_label = {'Subiculum','CA1','CA2','CA3','CA4/DG'};
[colorcode,colorcode2] = VPF_create_hippocampus_colorcode();

pre_vs_post_all_layers_masked = reshape(pre_vs_post_all_layers_masked,[N_subfields,N_subjects,N_layers]);
memory_vs_math_all_layers_masked = reshape(memory_vs_math_all_layers_masked,[N_subfields,N_subjects,N_layers]);


contrast_diff = memory_vs_math_all_layers_masked-pre_vs_post_all_layers_masked;
mean_contrast_diff = (squeeze(mean(contrast_diff,2))).';
SE_contrast_diff = (squeeze(std(contrast_diff,[],2))./(sqrt(N_subjects))).';
%%

system('Rscript ~/work/postdoc/projects/library/stats/VPF_compare_pre_vs_post_with_memory_vs_math_full_model_supplement.R');

load('/media/pfaffenrot/Elements/postdoc/projects/data/avg/memory/res_LMM_check_contrast_difference.mat');

Sub_struct = struct('p_value',Subiculum_p,'F_value',Subiculum_F,'coeff',Subiculum_coeff,'err',Subiculum_err);
CA1_struct = struct('p_value',CA1_p,'F_value',CA1_F,'coeff',CA1_coeff,'err',CA1_err);
CA2_struct = struct('p_value',CA2_p,'F_value',CA2_F,'coeff',CA2_coeff,'err',CA2_err);
CA3_struct = struct('p_value',CA3_p,'F_value',CA3_F,'coeff',CA3_coeff,'err',CA3_err);
CA4_struct = struct('p_value',CA4_p,'F_value',CA4_F,'coeff',CA4_coeff,'err',CA4_err);


results_LMM_check_contrast_effect = struct('Subiculum',Sub_struct,'CA1',CA1_struct,'CA2',CA2_struct,...
    'CA3',CA3_struct,'CA4DG',CA4_struct);
clear Sub_struct Subiculum_p Subiculum_F Subiculum_coeff Subiculum_err ...
    CA1_struct CA1_p CA1_F CA1_coeff CA1_err ...
    CA2_struct CA2_p CA2_F CA2_coeff CA2_err ...
    CA3_struct CA3_p CA3_F CA3_coeff CA3_err ...
    CA4_struct CA4_p CA4_F CA4_coeff CA4_err
%%
plotspecs = struct('FontName','Arial','FontSize',28,'colormap',jet(256),...
    'xtick',[],'xlim',[1 30],'LineWidth',3);
plotspecs.color = colorcode(:,1);

plotspecs.ytick = -1:0.25:1.5;
plotspecs.ylim = [-1 1.5];

plotspecs_eb = plotspecs;
plotspecs_eb.LineWidth = 4;

for subfield = 1:N_subfields
    figure()
    graph_to_show = squeeze(contrast_diff(subfield,:,:)).';
    h = VPF_show(@plot,1:30,graph_to_show,[],subfield_label{subfield},[],[],plotspecs);
    for kk = 1:length(h)
        set(h(kk), 'Color',[plotspecs.color{subfield},0.5]);
    end
    l = line([10,10],plotspecs.ylim,'LineWidth',plotspecs.LineWidth,'Color','black');
    hold on
    plotspecs_eb.color = 'k';
    shadedErrorBar(1:30,mean_contrast_diff(:,subfield),SE_contrast_diff(:,subfield),plotspecs_eb,1);

    if subfield == N_subfields
        plotspecs.Marker = '.';
        plotspecs.MarkerSize = 30;
        h = VPF_show(@plot,10,graph_to_show(10,:),[],[],[],[],plotspecs);
        for kk = 1:length(h)
            set(h(kk), 'Color',[plotspecs.color{subfield},0.5]);
        end
        h = VPF_show(@plot,N_layers,graph_to_show(end,:),[],subfield_label{subfield},[],[],plotspecs);
        for kk = 1:length(h)
            set(h(kk), 'Color',[plotspecs.color{subfield},0.5]);
        end
    end
    g = gca;
    set(g,'FontName',plotspecs.FontName,'FontSize',plotspecs.FontSize)
    my_ylabel = "(memory > math) - (pre > post) [z]";
    ylabel(my_ylabel);
end

%% plot LME fitting coefficients

beta = zeros(N_subfields-1,3);
err  = zeros(N_subfields-1,3);
for subfield = 1:N_subfields-1
    tmp = subfield_label{subfield};
    
    beta(subfield,1) = results_LMM_check_contrast_effect.(tmp).coeff.coeff_lin(2);
    beta(subfield,2) = results_LMM_check_contrast_effect.(tmp).coeff.coeff_quad(2);
    beta(subfield,3) = results_LMM_check_contrast_effect.(tmp).coeff.coeff_quad(3);

    err(subfield,1) = results_LMM_check_contrast_effect.(tmp).err.err_lin(2);
    err(subfield,2) = results_LMM_check_contrast_effect.(tmp).err.err_quad(2);
    err(subfield,3) = results_LMM_check_contrast_effect.(tmp).err.err_quad(3);
end


figure,
nil = zeros(4,1);
x = 1:4;
b1 = bar(x,[beta(:,1) nil nil],'grouped','LineWidth',2);
b1(1).FaceColor = 'flat';

g = gca;
set(g,'FontName','Arial','FontSize',28,'YLim',[-0.16 0.16])
ylabel('LME coefficient layer')
g.YAxis.Exponent = -2;


hold on
x_err = (1:4) -0.22;
er = errorbar(x_err,beta(:,1),err(:,1),'LineWidth', 2);
er.Color = [0 0 0];
er.LineStyle = 'none';

b2 = bar(x,[nil beta(:,2) nil],'grouped','LineWidth',2);
b2(2).FaceColor = 'flat';

x_err = (1:4);
er = errorbar(x_err,beta(:,2),err(:,2),'LineWidth', 2);
er.Color = [0 0 0];
er.LineStyle = 'none';



yyaxis right
b3 = bar(x,[nil nil beta(:,3)],'grouped','LineWidth',2);
b3(3).FaceColor = 'flat';

x_err = (1:4) +0.22;
er = errorbar(x_err,beta(:,3),err(:,3),'LineWidth', 2);
er(1).Color = [0 0 0];
er(1).LineStyle = 'none';

for ii = 1:4
    b1(1).CData(ii,:) = colorcode{ii,1};
    b2(2).CData(ii,:) = colorcode{ii,2};
    b3(3).CData(ii,:) = colorcode2{ii,2};
end

g = gca;
set(g,'xticklabel',{'Sub','CA1','CA2','CA3'})
set(g,'FontName','Arial','FontSize',24,'YLim',[-0.004 0.004])
ylabel('LME coefficient layer^2')
g.YAxis(2).Exponent = -3;
g.YAxis(2).Color = [128,128,128]/255;
title('LME fit')

h2 = bar(nan(3,3));
h2(1).FaceColor =  colorcode{1,1};
h2(2).FaceColor = colorcode{1,2};
h2(3).FaceColor = colorcode2{1,2};
legend(h2,{'layer linear model','layer quad. model','layer^2 quad. model'},'Location','northwest')
