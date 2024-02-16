clear;clc;close all;

data = dir('/media/pfaffenrot/Elements/postdoc/projects/data/7*');

colorcode = VPF_create_hippocampus_colorcode();
plotspecs = struct('FontName','Arial','FontSize',22,'colormap',jet(256),...
    'xtick',[],'xlim',[1 30],'LineWidth',2);
plotspecs.color = colorcode(:,1);
plotspecs.ytick = -2:0.5:2.5;
plotspecs.ylim = [-2 2.5];

plotspecs_eb = plotspecs;

nScan = length(data);
alph_FWE = 0.05;
bNeg = 0;
LAYERS = true;


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

X = X(:,:,[1 2 3 4 5 6 7 8 9]);

N_layers = size(X,2);
N_subjects = size(X,3);
N_subfields = size(X,1);
subfield_label = {'Subiculum','CA1','CA2','CA3','CA4/DG'};
%%



p_model = {};
for subfield = 1:5
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

    formula = 'beta ~ depths + (1|subject)';

    lme = fitlme(tbl,formula,'FitMethod','REML');
    p_model = cat(1,p_model,double(lme.Coefficients(:,end-2)));

    coeff = double(lme.Coefficients(:,2));
    CI = [double(lme.Coefficients(:,end-1)),double(lme.Coefficients(:,end))];

    x2 = 1:30;
    y_pred = coeff(2)*x2+coeff(1);
    up = CI(2,2)*x2+CI(1,2);

    y_pred(isnan(beta_to_show(:,1))) = nan;
    err = up-y_pred;

    figure()
    h = VPF_show(@plot,1:30,beta_to_show,[],subfield_label{subfield},[],[],plotspecs);
    for kk = 1:length(h)
        set(h(kk), 'Color',[plotspecs.color{subfield},0.5]);
    end
    l = line([10,10],plotspecs.ylim,'LineWidth',plotspecs.LineWidth,'Color','black');
    hold on
    plotspecs_eb.color = 'k';
    shadedErrorBar(x2,mean(beta_to_show,2),std(beta_to_show,[],2)./sqrt(N_subjects),plotspecs_eb,1);

    %     plot(x2,y_pred,'LineWidth',plotspecs.LineWidth,'Color',[[200,33,55]./255,1]);
    if subfield == N_subfields
        plotspecs.Marker = '.';
        plotspecs.MarkerSize = 30;
        h = VPF_show(@plot,10,beta_to_show(10,:),[],[],[],[],plotspecs);
        for kk = 1:length(h)
            set(h(kk), 'Color',[plotspecs.color{subfield},0.5]);
        end
        h = VPF_show(@plot,N_layers,beta_to_show(end,:),[],subfield_label{subfield},[],[],plotspecs);
        for kk = 1:length(h)
            set(h(kk), 'Color',[plotspecs.color{subfield},0.5]);
        end
    end    
    g = gca;
    set(g,'FontName',plotspecs.FontName,'FontSize',plotspecs.FontSize)
    my_ylabel = '\beta_{pre}-\beta_{post} [z]';
    ylabel(my_ylabel);
end
p_model = cell2mat(p_model);
p_model = reshape(multicmp(p_model,'up',0.05),[2 N_subfields]);
