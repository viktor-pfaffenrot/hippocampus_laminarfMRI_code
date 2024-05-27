clear;clc;%close all;

pre_vs_post_flag = true;
DO_MATLAB_LME = true;

data = dir('/media/pfaffenrot/Elements/postdoc/projects/data/7*');

[colorcode,colorcode2] = VPF_create_hippocampus_colorcode();

nScan = length(data);



if pre_vs_post_flag == true
    for Scan = 1:nScan
        if strcmp(data(Scan).name,'7495') ||  strcmp(data(Scan).name,'7566')
            for session = 1:2
                load([data(Scan).folder '/' data(Scan).name '/memory/ses-0' num2str(session)...
                    '/z_transformed/results_pre_vs_post_z_vessel_masked.mat']);
                if session == 1
                    tmp = zeros(size(results_pre_vs_post.con_array));
                    if Scan == 1
                        X = zeros([size(results_pre_vs_post.con_array),nScan]);
                    end
                end
                tmp = tmp+results_pre_vs_post.con_array;
                X(:,:,Scan) = tmp./2;
            end
        else
            load([data(Scan).folder '/' data(Scan).name...
                '/memory/z_transformed/results_pre_vs_post_z_vessel_masked.mat']);
            if Scan == 1
                X = zeros([size(results_pre_vs_post.con_array),nScan]);
            end
            X(:,:,Scan) = results_pre_vs_post.con_array;
        end
    end

else
    for Scan = 1:nScan
        if strcmp(data(Scan).name,'7495') ||  strcmp(data(Scan).name,'7566')
            for session = 1:2
                load([data(Scan).folder '/' data(Scan).name '/memory/ses-0' num2str(session)...
                    '/z_transformed/results_memory_vs_math_z_vessel_masked.mat']);
                if session == 1
                    tmp = zeros(size(results_memory_vs_math.con_array));
                    if Scan == 1
                        X = zeros([size(results_memory_vs_math.con_array),nScan]);
                    end
                end
                tmp = tmp+results_memory_vs_math.con_array;
                X(:,:,Scan) = tmp./2;
            end
        else
            load([data(Scan).folder '/' data(Scan).name...
                '/memory/z_transformed/results_memory_vs_math_z_vessel_masked.mat']);
            if Scan == 1
                X = zeros([size(results_memory_vs_math.con_array),nScan]);
            end
            X(:,:,Scan) = results_memory_vs_math.con_array;
        end
    end
end

X = X(:,:,[1 2 3 4 5 6 7 8 9]);

N_layers = size(X,2);
N_subjects = size(X,3);
N_subfields = size(X,1);
subfield_label = {'Subiculum','CA1','CA2','CA3','CA4/DG'};


%also load breath-hold data to plot on top of the AM results
load('/media/pfaffenrot/Elements/postdoc/projects/data/avg/breathhold/results_avg.mat')

results_BH = results_avg;
clear results_avg;
%%
if pre_vs_post_flag == true
    R_flag = 'TRUE';
else
    R_flag = 'FALSE';
end

system(['Rscript /home/pfaffenrot/work/postdoc/projects/library/stats/VPF_LMM_check_main_effect.R ' R_flag]);

if pre_vs_post_flag==true
    load('/media/pfaffenrot/Elements/postdoc/projects/data/avg/memory/res_LMM_check_main_effect_pre_vs_post.mat');
else
    load('/media/pfaffenrot/Elements/postdoc/projects/data/avg/memory/res_LMM_check_main_effect_memory_vs_math.mat');
end
Sub_struct = struct('p_value',Subiculum_p,'F_value',Subiculum_F,'coeff',Subiculum_coeff,'err',Subiculum_err);
CA1_struct = struct('p_value',CA1_p,'F_value',CA1_F,'coeff',CA1_coeff,'err',CA1_err);
CA2_struct = struct('p_value',CA2_p,'F_value',CA2_F,'coeff',CA2_coeff,'err',CA2_err);
CA3_struct = struct('p_value',CA3_p,'F_value',CA3_F,'coeff',CA3_coeff,'err',CA3_err);
CA4_struct = struct('p_value',CA4_p,'F_value',CA4_F,'coeff',CA4_coeff,'err',CA4_err);


results_LMM_check_main_effect = struct('Subiculum',Sub_struct,'CA1',CA1_struct,'CA2',CA2_struct,...
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

if pre_vs_post_flag == true
    plotspecs.ytick = -2:0.5:2.5;
    plotspecs.ylim = [-2 2.5];

else
    plotspecs.ytick = -0.5:0.5:2.5;
    plotspecs.ylim = [-0.5 2.5];
end


plotspecs_eb = plotspecs;
plotspecs_eb.LineWidth = 4;

coeff_out = zeros(5,1);
err_out = zeros(5,1);
if DO_MATLAB_LME == true
    % The code below shows individual subject profiles as a spaghetti plot. The
    % mean and SE over subjects is plotted on top. The rest of the code implements
    % an LME model similarly to that in R for comparison purposes. The main difference is
    % that matlab does not allow for a model comparison using the Kenward-Roger df
    % approximation for small sample sized, therefore the R implementation
    p_model = {};
    for subfield = 1:N_subfields
        subject = {};
        for ii = 1:N_subjects
            subject = cat(1,subject,repmat({['subject_' num2str(ii) ]},[size(X,2) 1]));
        end
        depths = repmat((1:N_layers).',[N_subjects,1]);

        beta_to_show = squeeze(X(subfield,:,:));
        beta = beta_to_show(:);
        depths(isnan(beta)) = [];
        subject(isnan(beta)) = [];
        beta(isnan(beta)) = [];

        tbl = table(subject,depths,beta);
        tbl.subject = categorical(tbl.subject);


        if subfield == 2 && ~pre_vs_post_flag
            formula = 'beta ~  depths^2 + depths + (1|subject)';
        else
            formula = 'beta ~ depths + (1|subject)';
        end
        lme = fitlme(tbl,formula,'FitMethod','REML');
        anv = anova(lme,'DFMethod','satterthwaite');
        %     p_model = cat(1,p_model,double(lme.Coefficients(:,end-2)));
        p_model = cat(1,p_model,anv.pValue);

        %         coeff = double(lme.Coefficients(:,2));
        if subfield == N_subfields
            ss = 'CA4DG';
        else
            ss = subfield_label{subfield};
        end

        coeff = results_LMM_check_main_effect.(ss).coeff;
        coeff_out(subfield) = coeff.coeff_lin(2);
        err_b = results_LMM_check_main_effect.(ss).err;
        err_out(subfield) =  err_b.err_lin(2);
        CI = [double(lme.Coefficients(:,end-1)),double(lme.Coefficients(:,end))];

        x2 = 1:30;
        if subfield == 2 && ~pre_vs_post_flag
            p_model{subfield} = p_model{subfield}(1:2:end);
            y_pred = coeff.coeff_quad(2)*x2 + coeff.coeff_quad(3)*x2.^2+coeff.coeff_quad(1);
            up = CI(3,2)*x2.^2 + CI(2,2)*x2+CI(1,2);
            err = double(err_b.err_quad(3))*x2.^2 +...
                double(err_b.err_quad(2))*x2 + ...
                double(err_b.err_quad(1));

        else
            y_pred = coeff.coeff_lin(2)*x2+coeff.coeff_lin(1);
            up = CI(2,2)*x2+CI(1,2);
            err = double(err_b.err_lin(2))*x2+double(err_b.err_lin(1));
        end
        y_pred(isnan(beta_to_show(:,1))) = nan;
        %     keyboard
        mean_beta = mean(beta_to_show,2,'omitnan');
        SE_beta  =std(beta_to_show,[],2,'omitnan')./sqrt(N_subjects);
        if subfield == N_subfields
            plotspecs.LineStyle = 'none';
        else
            plotspecs.LineStyle = '-';
        end
        figure()
        h = VPF_show(@plot,1:30,beta_to_show,[],subfield_label{subfield},[],[],plotspecs);
        for kk = 1:length(h)
            set(h(kk), 'Color',[plotspecs.color{subfield},0.5]);
        end


        hold on
        if subfield ~=N_subfields
            l = line([10,10],plotspecs.ylim,'LineWidth',plotspecs.LineWidth,'Color','black');
            plotspecs_eb.color = 'k';
            shadedErrorBar(x2,mean_beta,SE_beta,plotspecs_eb,1); % y_pred,err
        else
            %             l = line([10,10],plotspecs.ylim,'LineWidth',plotspecs.LineWidth,'Color','black');
            plotspecs.Marker = '.';
            plotspecs.MarkerSize = 50;
            h = VPF_show(@plot,10,beta_to_show(10,:),[],[],[],[],plotspecs);
            for kk = 1:length(h)
                set(h(kk), 'Color',[plotspecs.color{subfield},0.5]);
            end
            h2 = errorbar(10,mean_beta(10),SE_beta(10),'LineWidth',plotspecs_eb.LineWidth,'Marker','.','MarkerSize',plotspecs.MarkerSize);
            set(h2, 'Color',plotspecs_eb.color);


            h = VPF_show(@plot,N_layers,beta_to_show(end,:),[],subfield_label{subfield},[],[],plotspecs);
            for kk = 1:length(h)
                set(h(kk), 'Color',[plotspecs.color{subfield},0.5]);
            end
            h2 = errorbar(30,mean_beta(end),SE_beta(end),'LineWidth',plotspecs_eb.LineWidth,'Marker','.','MarkerSize',plotspecs.MarkerSize);
            set(h2, 'Color',plotspecs_eb.color);


        end
        %     plot(x2,y_pred,'LineWidth',plotspecs.LineWidth,'Color',[[200,33,55]./255,1]);
        g = gca;
        set(g,'FontName',plotspecs.FontName,'FontSize',plotspecs.FontSize)
        if pre_vs_post_flag == true
            my_ylabel = '\beta_{pre}-\beta_{post} [z]';
        else
            my_ylabel = '\beta_{memory}-\beta_{math} [z]';
        end
        ylabel(my_ylabel);

        if subfield ~=N_subfields
            yyaxis right
            h_bh = VPF_show(@plot,1:30,results_BH.dS_mean_echo_weighted_submean(:,subfield),[],subfield_label{subfield},[],[],plotspecs);
            set(h_bh,'Color',colorcode{subfield,2});
            set(h_bh,'LineStyle','--')
            g.YAxis(2).Color = colorcode{subfield,2};
            ylabel('\DeltaS_{weighted} [a.u.]')
        end
    end

    p_model = cell2mat(p_model);
    p_model = reshape(multicmp(p_model,'up',0.05),[2 N_subfields]);
end
%% plot LME fitting coefficients

beta = zeros(N_subfields-1,3);
err  = zeros(N_subfields-1,3);
for subfield = 1:N_subfields-1
    tmp = subfield_label{subfield};

    beta(subfield,1) = results_LMM_check_main_effect.(tmp).coeff.coeff_lin(2);
    beta(subfield,2) = results_LMM_check_main_effect.(tmp).coeff.coeff_quad(2);
    beta(subfield,3) = results_LMM_check_main_effect.(tmp).coeff.coeff_quad(3);

    err(subfield,1) = results_LMM_check_main_effect.(tmp).err.err_lin(2);
    err(subfield,2) = results_LMM_check_main_effect.(tmp).err.err_quad(2);
    err(subfield,3) = results_LMM_check_main_effect.(tmp).err.err_quad(3);
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


