clear;clc;close all

data = dir('/media/pfaffenrot/PostDoc_data/projects/data/7*');
data = data([1 2 4 5 6 7 8 9]);
subjects = length(data);

subfields = 5;

subids = cell(length(data),1);
for subject = 1:subjects
    subids{subject} = data(subject).name;
    if subject == 1
        load([data(subject).folder '/' data(subject).name '/' ...
            '/breathhold/vessel_density.mat']);
        sz = size(vessel_density);
        vessel_densities = zeros([sz subjects]);
    else
        load([data(subject).folder '/' data(subject).name '/' ...
            '/breathhold/vessel_density.mat']);
    end
    vessel_densities(:,:,:,subject) = vessel_density;

end
select = [1 2 3 4 5 6 7 8];

%%
[colorcode,colorcode_points] = VPF_create_hippocampus_colorcode();
p = zeros(5,1);
for subfield = 1:5
    if subfield == 1
        [p(subfield),~,tmp] = signrank(squeeze(vessel_densities(1,subfield,1,:)),squeeze(vessel_densities(1,subfield,2,:)),'alpha',0.05);
        stats = repmat(tmp,[5 1]);
    else
        [p(subfield),~,stats(subfield)] = signrank(squeeze(vessel_densities(1,subfield,1,:)),squeeze(vessel_densities(1,subfield,2,:)),'alpha',0.05);
    end
end
p = multicmp(p,'up',0.05);
%%
m = median(vessel_densities,4,'omitnan');
s = std(vessel_densities,[],4,'omitnan');
figure,
b = bar(1:5,squeeze(m(1,:,:)),'LineWidth',2);
ylim([0 550]);
b(1).FaceColor = 'flat';
b(2).FaceColor = 'flat';
b(1).CData(1,:) = colorcode{1,1};
b(2).CData(1,:) = colorcode{1,2};

b(1).CData(2,:) = colorcode{2,1};
b(2).CData(2,:) = colorcode{2,2};

b(1).CData(3,:) = colorcode{3,1};
b(2).CData(3,:) = colorcode{3,2};

b(1).CData(4,:) = colorcode{4,1};
b(2).CData(4,:) = colorcode{4,2};

b(1).CData(5,:) = colorcode{5,1};
b(2).CData(5,:) = colorcode{5,2};

hold on
for ii = 1:sz(3)
    if ii == 1
        x = (1:sz(2)) -0.15;
    else
        x = (1:sz(2)) + 0.15 * (ii- 1);
    end
    er = errorbar(x,squeeze(m(1,:,ii)),squeeze(s(1,:,ii)),'LineWidth', 2);
    er(1).Color = [0 0 0];
    er(1).LineStyle = 'none';

end

g = gca;
set(g,'xticklabel',{'Sub','CA1','CA2','CA3','DG/CA4'})
set(g,'FontName','Arial','FontSize',24)
ylabel('vessel density \rho [a.u.]')
hold on
for ii = 1:sz(3)
    if ii == 1
        x = (1:sz(2)) -0.15;
    else
        x = (1:sz(2)) + 0.15 * (ii- 1);
    end
    sc = scatter(x, squeeze(vessel_densities(1,:,ii,:)), 'filled', 'MarkerEdgeColor', 'k');
    for jj = 1:length(sc)
        newcdata = zeros([prod(sz(2)) 3]);

        for uu = 1:sz(2)
            newcdata(uu,:) = colorcode_points{uu,ii};
        end
        sc(jj).CData = newcdata;
    end
end
%%
subfield_labels = {'Subiculum','CA1','CA2','CA3','DG/CA4'};
vessel_densities_selected = vessel_densities(:,:,:,select);

for subject = 1:subjects
    load([data(subject).folder '/' data(subject).name '/' ...
        '/breathhold/results_breathhold.mat']);
    if subject == 1
        results_breathhold_all = repmat(results_breathhold,[subject 1]);
        results_breathhold_all(1) = results_breathhold;
    else
        results_breathhold_all(subject) = results_breathhold;
    end
end



d_vessel_densities = squeeze(vessel_densities_selected(:,:,1,:))-...
    squeeze(vessel_densities_selected(:,:,2,:));

y = zeros([subjects,subfields]);

for subject = 1:subjects
    y(subject,:) = results_breathhold_all(subject).dS_mean_echo_weighted_slope;
end

y = y(select,:);
x = squeeze(d_vessel_densities(1,:,:)).';


figure,
for jj = 1:subfields
    sc(jj) = scatter(x(:,jj),y(:,jj), 80, 'filled', 'MarkerEdgeColor', 'k');
    sc(jj).CData = colorcode{jj,1};
    hold on
    %k = convhull(x(:,jj),y(:,jj));

    %plot(x(k,jj),y(k,jj),'Color',colorcode{jj,1},'LineWidth',1,'LineStyle','-')
end

system('Rscript /home/pfaffenrot/work/postdoc/projects/library/stats/VPF_LMM_dSlope_vs_dDensity.R ');
load('/media/pfaffenrot/Elements/postdoc/projects/data/avg/breathhold/results_stat_density_change_vs_slope_change.mat')
results_LMM = struct('KR',KR,'Estimates',Estimates,'SE',SE);
clear KR Estimates SE

x2 = linspace(min(x(:)),max(x(:)),numel(x));
y_pred = results_LMM.Estimates(2)*x2+results_LMM.Estimates(1);
err = results_LMM.SE(2)*x2+results_LMM.SE(1);

plotspecs = struct('FontName','Arial','FontSize',22,'colormap',jet(256),...
    'xtick',[-500:100:300],'xlim',[-500 300],'LineWidth',2);

plotspecs.ytick = -0.08:0.02:0.1;
plotspecs.ylim = [-0.08 0.1];

plotspecs.color = 'k';
shadedErrorBar(x2,y_pred,err,plotspecs,1);



plot(x2,y_pred,'LineWidth',2,'LineStyle','-','Color','black','Marker','none')
ylim([-0.1 0.1])

g = gca;
set(g,'FontName','Arial','FontSize',24)
ylabel('signal change slope [a.u.]')
xlabel('\rho_{inner}-\rho_{outer} [a. u.]')
legend(sc(1:subfields),subfield_labels,'Location','NorthWest')


%%
dS_mean_echo = zeros([size(results_breathhold_all(1).dS_mean_echo),subjects]);
dS_mean_echo_weighted = dS_mean_echo;
T2s = dS_mean_echo;
dR2s = T2s;
for subject = 1:subjects
    dS_mean_echo(:,:,subject) = results_breathhold_all(subject).dS_mean_echo;
    dS_mean_echo_weighted(:,:,subject) = results_breathhold_all(subject).dS_mean_echo_weighted;
    T2s(:,:,subject) = 1./results_breathhold_all(subject).R2s(:,:,1)*1000;
    dR2s(:,:,subject) = results_breathhold_all(subject).dR2s;

end

%%
vessel_density_submean = squeeze(m(1,:,:));
vessel_density_substd = squeeze(s(1,:,:));

dS_mean_echo_submean = mean(dS_mean_echo(:,:,select),3);
dS_mean_echo_substd = std(dS_mean_echo(:,:,select),[],3);

dS_mean_echo_weighted_submean = mean(dS_mean_echo_weighted(:,:,select),3);
dS_mean_echo_weighted_substd = std(dS_mean_echo_weighted(:,:,select),[],3);



T2s_submean = mean(T2s(:,:,select),3);
T2s_substd = std(T2s(:,:,select),[],3);

dR2s_submean = mean(dR2s(:,:,select),3);
dR2s_substd = std(dR2s(:,:,select),[],3);



results_avg = struct('dS_mean_echo_submean',dS_mean_echo_submean,'dS_mean_echo_substd',dS_mean_echo_substd,...
    'dS_mean_echo_weighted_submean',dS_mean_echo_weighted_submean,'dS_mean_echo_weighted_substd',dS_mean_echo_weighted_substd,...
    'T2s_submean',T2s_submean,'T2s_substd',T2s_substd,...
    'dR2s_submean',dR2s_submean,'dR2s_substd',dR2s_substd,...
    'vessel_density_submean',vessel_density_submean,'vessel_density_substd',vessel_density_substd);
%%
outpath = '/media/pfaffenrot/Elements/postdoc/projects/data/avg/breathhold';

fulloutname = [outpath '/results_avg.mat'];

save(fulloutname,'results_avg');

results_json = jsonencode(results_avg,PrettyPrint=true);
fid = fopen([outpath '/results_avg.json'],'w');

fprintf(fid,'%s',results_json);
fclose(fid);
