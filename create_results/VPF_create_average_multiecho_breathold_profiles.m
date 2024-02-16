clear;clc;

data = dir('/media/pfaffenrot/Elements/postdoc/projects/hippocampus_breathhold/FLASH/derivatives/pipeline/7*');
data = data([1 3 5 6 8:end]);
sessionid = [1 2 3 2 1 1 1 2];
TE = [2 6 10 14 17 20]; %
N = [20 10]; %number of layers

N_subjects = length(data);

for subject = 1:N_subjects
    load([data(subject).folder '/' data(subject).name '/ses-0' num2str(sessionid(subject)) '/layerfication/OFF_layers.mat']);
    load([data(subject).folder '/' data(subject).name '/ses-0' num2str(sessionid(subject)) '/layerfication/ON_layers.mat']);
    if subject == 1
        OFF_all = zeros([size(OFF) N_subjects]);
        ON_all = OFF_all;
    end
    OFF_all(:,:,:,subject) = OFF;
    ON_all(:,:,:,subject) = ON;
end
% sOFF = std(OFF_all,[],4,'omitnan');
% sON = std(ON_all,[],4,'omitnan');
% OFF = mean(OFF_all,4,'omitnan');
% ON = mean(ON_all,4,'omitnan');

d = ON_all-OFF_all;
sd = std(d,[],4,'omitnan');
md = mean(d,4,'omitnan');
figtitles = [cellstr('Subiculum'), cellstr('CA1'),cellstr('CA2'),...
    cellstr('CA3'),cellstr('CA4/DG')]; % cellstr('all')


% Desired range for the mean
target_mean_range = [-1, 1];
for subfield = 1:length(figtitles)
    % Calculate the scale factor for mean rescaling
    tmp = md(:,:,subfield);
    mean_scale_factor = diff(target_mean_range) / (2 * std(tmp(:),'omitnan'));

    % Rescale the mean

    md(:,:,subfield) = (tmp - mean(tmp(:),'omitnan')) * mean_scale_factor;

    % Calculate the scale factor for standard deviation rescaling
    sd(:,:,subfield) = sd(:,:,subfield) * mean_scale_factor;

    tmp = md(:,:,subfield);
    fac = max(tmp(:));
    md(:,:,subfield) = md(:,:,subfield)./fac;
    sd(:,:,subfield) = sd(:,:,subfield)./fac;
end

sd = sd./N_subjects;


mylegend = cell(1,numel(TE));
count = 1;
for ii = 1:numel(TE)
    mylegend{count} = [num2str(TE(ii)) ' ms'];
    count = count + 1;
end



plotspecs = struct('FontName','Arial','FontSize',22,'colormap',jet(256),...
    'xtick',[],'xlim',[1 30],'LineWidth',2);
plotspecs.ytick = -1.8:0.2:1.8;
plotspecs.ylim = [-1.8 1.8];
% plotspecs.ytick = 4:2:20;
% plotspecs.ylim = [4 20];




pos = 30;
for subfield = 1:numel(figtitles)
    figure(),
    for kk = 1:numel(TE)
        set(gca,'ColorOrderIndex',kk);
        h(kk) = shadedErrorBar(1:pos,md(:,kk,subfield),...
            sd(:,kk,subfield),plotspecs,1);
        if kk == 1
            legendax = h(kk).mainLine;
            hold on
        else
            legendax = cat(2,legendax,h(kk).mainLine);
        end
    end

    if subfield ==numel(figtitles)
        plotspecs.Marker = '.';
        plotspecs.MarkerSize = 30;
        set(gca,'ColorOrderIndex',1);
        VPF_show(@plot,10,md(10,:,subfield),[],[],[],'\DeltaS_{norm} [a.u.]',plotspecs);
        set(gca,'ColorOrderIndex',1);
        VPF_show(@plot,sum(N),md(end,:,subfield),[],[],[],'\DeltaS_{norm} [a.u.]',plotspecs);
    end

    l = line([10,10],plotspecs.ylim,'LineWidth',plotspecs.LineWidth,'Color','black');
    lgd = legend(legendax,mylegend,'Location','NorthEast');
    lgd.NumColumns = 3;

    title([figtitles{subfield}])
    ylabel('\DeltaS_{norm,mean} [a.u.]')
    legend(legendax,mylegend,'Location','NorthEast');
    ColorOrder = get(gca,'ColorOrder');
    set(gca,'FontSize',plotspecs.FontSize)
end
%%
w = 1./TE;
% w = 1;
w = w./sum(w);
dw = bsxfun(@times,d,reshape(w,[1 numel(TE) 1 1]));
md = squeeze(mean(mean(dw,2),4));
sd = squeeze(std(mean(dw,2),[],4));
sd = sd./N_subjects;
colorcode = VPF_create_hippocampus_colorcode();
% plotspecs.ytick = 4:2:20;
% plotspecs.ylim = [4 20];
plotspecs.ytick = 0.8:0.2:6;
plotspecs.ylim = [0.8 2.6];
figure(),
for subfield = 1:length(figtitles)
    plotspecs.color = colorcode{subfield,1};
    h(subfield) = shadedErrorBar(1:pos,md(:,subfield),...
        sd(:,subfield),plotspecs,1);
    if subfield == 1
        legendax = h(subfield).mainLine;
        hold on
    else
        legendax = cat(2,legendax,h(subfield).mainLine);
    end
end
lgd = legend(legendax,figtitles,'Location','NorthEast');
l = line([10,10],plotspecs.ylim,'LineWidth',plotspecs.LineWidth,'Color','black');

ylabel('\DeltaS_{echomean,weighted} [a.u.]')
legend(legendax,figtitles,'Location','NorthEast');
ColorOrder = get(gca,'ColorOrder');
set(gca,'FontSize',plotspecs.FontSize)



md = squeeze(mean(mean(d,2),4));
sd = squeeze(std(mean(d,2),[],4));
sd = sd./N_subjects;
colorcode = VPF_create_hippocampus_colorcode();
plotspecs.ytick = 4:2:20;
plotspecs.ylim = [4 20];
figure(),
for subfield = 1:length(figtitles)
    plotspecs.color = colorcode{subfield,1};
    h(subfield) = shadedErrorBar(1:pos,md(:,subfield),...
        sd(:,subfield),plotspecs,1);
    if subfield == 1
        legendax = h(subfield).mainLine;
        hold on
    else
        legendax = cat(2,legendax,h(subfield).mainLine);
    end
end
    lgd = legend(legendax,figtitles,'Location','NorthEast');
    l = line([10,10],plotspecs.ylim,'LineWidth',plotspecs.LineWidth,'Color','black');

    ylabel('\DeltaS_{echomean} [a.u.]')
    legend(legendax,figtitles,'Location','NorthEast');
    ColorOrder = get(gca,'ColorOrder');
    set(gca,'FontSize',plotspecs.FontSize)

