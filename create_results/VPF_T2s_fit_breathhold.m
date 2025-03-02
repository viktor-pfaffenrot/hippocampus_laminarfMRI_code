clear;clc;close all;
subid = '7218';
sessionid = '01';
layerpath = ['/media/pfaffenrot/Elements/postdoc/projects/hippocampus_breathhold/FLASH/derivatives/pipeline/' subid '/ses-' sessionid '/layerfication'];
func_path = ['/media/pfaffenrot/Elements/postdoc/projects/hippocampus_breathhold/FLASH/derivatives/pipeline/' subid '/ses-' sessionid '/func_bart'];
load([layerpath '/OFF_layers.mat'])
load([layerpath '/ON_layers.mat'])
TE = [2 6 10 14 17 20]; %17 21
TEs = numel(TE);
subfields = 5;
layers = 30;
conditions = 2;

dS = ON-OFF;
dS_mean_echo = squeeze(mean(dS,2));


w = 1./TE;
w = w./sum(w);
dS_mean_echo_weighted = squeeze(mean(bsxfun(@times,dS,w),2));
%%
opts = optimset('lsqcurvefit');
opts.Display = 'none';
opts.MaxIter = 1e5;
opts.TolX = 1e-10;
opts.TolFun = 1e-10;

x0 = [1000 30];
lb = [0 0];
ub = [Inf Inf];
NLinModelFun = @(x,t)(x(1)*exp(-t/x(2)));

sz = size(OFF);
T2s = zeros([sz(1) sz(3) 2]);
R2_fit = T2s;

OFF(isnan(OFF)) = 0;
ON(isnan(ON)) = 0;

for layer = 1:layers
    for subfield = 1:subfields
        out = lsqcurvefit(NLinModelFun,x0,TE,OFF(layer,:,subfield),lb,ub,opts);
        T2s(layer,subfield,1) = out(2);

        imgfit = out(1).*exp(-TE./out(2));
        SSE = sum((squeeze(OFF(layer,:,subfield))-imgfit).^2);
        SST = sum((squeeze(OFF(layer,:,subfield))-mean(squeeze(OFF(layer,:,subfield)))).^2);
        R2_fit(layer,subfield,1) = 1-SSE/SST;

        out = lsqcurvefit(NLinModelFun,x0,TE,squeeze(ON(layer,:,subfield)),lb,ub,opts);
        T2s(layer,subfield,2) = out(2);
        imgfit = out(1).*exp(-TE./out(2));
        SSE = sum((squeeze(ON(layer,:,subfield))-imgfit).^2);
        SST = sum((squeeze(ON(layer,:,subfield))-mean(squeeze(ON(layer,:,subfield)))).^2);
        R2_fit(layer,subfield,2) = 1-SSE/SST;
    end
end
T2s(1:9,[1,4,5],:) = nan;
R2_fit(1:9,[1,4,5],:) = nan;


%%
T2s = T2s./1000;
R2s = 1./T2s;
dR2s = R2s(:,:,1)-R2s(:,:,2);

%%
colorcode = VPF_create_hippocampus_colorcode();
plotspecs = struct('FontName','Arial','FontSize',22,'colormap',jet(256),...
    'xtick',[],'xlim',[1 30],'LineWidth',2);
plotspecs.color = colorcode(:,1);
plotspecs.ytick = 20:2:60;
plotspecs.ylim = [20 60];


subfield_names = [cellstr('Subiculum'), cellstr('CA1'), cellstr('CA2'),...
    cellstr('CA3'),cellstr('CA4/DG')];
figure(),
h = VPF_show(@plot,1:30,squeeze(T2s(:,:,1)*1000),[],[],[],'T_2^* [ms]',plotspecs);
for kk = 1:subfield
    set(h(kk), 'Color',plotspecs.color{kk});
end
hold on
plotspecs.Marker = '.';
plotspecs.MarkerSize = 30;
h = VPF_show(@plot,10,squeeze(T2s(10,end,1)*1000),[],[],[],'T_2^* [ms]',plotspecs);
set(h, 'Color',plotspecs.color{end});
h = VPF_show(@plot,size(T2s,1),squeeze(T2s(end,end,1)*1000),[],[],[],'T_2^* [ms]',plotspecs);
set(h, 'Color',plotspecs.color{end});
hold off
l = line([10,10],plotspecs.ylim,'LineWidth',plotspecs.LineWidth,'Color','black');
lgd = legend(subfield_names);
%%
plotspecs.Marker = 'none';
plotspecs.ytick = -4.5:0.5:6;
plotspecs.ylim = [-4.5 6];

figure(),
h = VPF_show(@plot,1:30,dR2s,[],[],[],'\DeltaR_2^* [1/s]',plotspecs);
for kk = 1:subfield
    set(h(kk), 'Color',plotspecs.color{kk});
end
hold on
plotspecs.Marker = '.';
plotspecs.MarkerSize = 30;
set(gca,'ColorOrderIndex',numel(subfield_names));
h = VPF_show(@plot,10,dR2s(10,end),[],[],[],'\DeltaR_2^* [1/s]',plotspecs);
set(h, 'Color',plotspecs.color{end});
h = VPF_show(@plot,size(dR2s,1),dR2s(end,end),[],[],[],'\DeltaR_2^* [1/s]',plotspecs);
set(h, 'Color',plotspecs.color{end});
hold off
l = line([10,10],plotspecs.ylim,'LineWidth',plotspecs.LineWidth,'Color','black');
lgd = legend(subfield_names);
%%

dR2s_slope = (dR2s(10,:)-dR2s(end,:))./(20-1);
dS_mean_echo_slope = (dS_mean_echo(10,:)-dS_mean_echo(end,:))./(20-1);
dS_mean_echo_weighted_slope = (dS_mean_echo_weighted(10,:)-dS_mean_echo_weighted(end,:))./(20-1);
plotspecs.Marker = 'none';
plotspecs.ytick = -40:5:60;
plotspecs.ylim = [-40 60];
figure(),
h = VPF_show(@plot,1:30,dS_mean_echo,[],[],[],'\DeltaR_2^* [1/s]',plotspecs);
for kk = 1:subfield
    set(h(kk), 'Color',plotspecs.color{kk});
end
hold on
plotspecs.Marker = '.';
plotspecs.MarkerSize = 30;
h = VPF_show(@plot,10,dS_mean_echo(10,end),[],[],[],'\DeltaS_{mean} [a.u.]',plotspecs);
set(h, 'Color',plotspecs.color{end});
h = VPF_show(@plot,size(dS_mean_echo,1),dS_mean_echo(end,end),[],[],[],'\DeltaS_{mean} [a.u.]',plotspecs);
set(h, 'Color',plotspecs.color{end});
hold off
l = line([10,10],plotspecs.ylim,'LineWidth',plotspecs.LineWidth,'Color','black');
lgd = legend(subfield_names);
%%
plotspecs.Marker = 'none';
plotspecs.ytick = -6:1:6;
plotspecs.ylim = [-6 6];
figure(),
h = VPF_show(@plot,1:30,dS_mean_echo_weighted,[],[],[],'\DeltaR_2^* [1/s]',plotspecs);
for kk = 1:subfield
    set(h(kk), 'Color',plotspecs.color{kk});
end
hold on
plotspecs.Marker = '.';
plotspecs.MarkerSize = 30;
h = VPF_show(@plot,10,dS_mean_echo_weighted(10,end),[],[],[],'\DeltaS_{mean,weighted} [a.u.]',plotspecs);
set(h, 'Color',plotspecs.color{end});
h = VPF_show(@plot,size(dS_mean_echo_weighted,1),dS_mean_echo_weighted(end,end),[],[],[],'\DeltaS_{mean,weighted} [a.u.]',plotspecs);
set(h, 'Color',plotspecs.color{end});
hold off
l = line([10,10],plotspecs.ylim,'LineWidth',plotspecs.LineWidth,'Color','black');
lgd = legend(subfield_names);
%%

results_breathhold = {dS_mean_echo,dS_mean_echo_weighted,R2s,dR2s,R2_fit,dR2s_slope,dS_mean_echo_slope,dS_mean_echo_weighted_slope};
results_breathhold = cell2struct(results_breathhold,{'dS_mean_echo','dS_mean_echo_weighted','R2s','dR2s','R2_fit','dR2s_slope','dS_mean_echo_slope','dS_mean_echo_weighted_slope'},2);
save([layerpath '/results_breathhold.mat'],'results_breathhold','-v7.3')
results_breathhold = jsonencode(results_breathhold,PrettyPrint=true);
fid = fopen([layerpath '/results_breathhold.json'],'w');
fprintf(fid,'%s',results_breathhold);
fclose(fid);

%% with vertices

load([layerpath '/OFF_layers_w_vertices.mat'])
OFF = OFF_w_vertices;
load([layerpath '/ON_layers_w_vertices.mat'])
ON = ON_w_vertices;

opts = optimset('lsqcurvefit');
opts.Display = 'none';
opts.MaxIter = 1e5;
opts.TolX = 1e-10;
opts.TolFun = 1e-10;

x0 = [1000 30];
lb = [0 0];
ub = [Inf Inf];
NLinModelFun = @(x,t)(x(1)*exp(-t/x(2)));

sz = size(OFF);
T2s = zeros([sz(1) 2]);
R2_fit = T2s;

OFF(isnan(OFF)) = 0;
ON(isnan(ON)) = 0;

for layer = 1:sz(1)
        out = lsqcurvefit(NLinModelFun,x0,TE,squeeze(OFF(layer,1,:)).',lb,ub,opts);
        T2s(layer,1) = out(2);

        out = lsqcurvefit(NLinModelFun,x0,TE,squeeze(OFF(layer,2,:)).',lb,ub,opts);
        T2s(layer,2) = out(2);
end
save([func_path '/T2s_from_vertices.mat'],'T2s');
%%

dS = ON_w_vertices - OFF_w_vertices;
dS_mean_echo = squeeze(mean(dS,3));
dS_mean_echo_weighted = squeeze(mean(bsxfun(@times,dS,reshape(w,[1 1 6])),3));


SURFS = {'inner','outer'};
mypath =['/media/pfaffenrot/Elements/postdoc/projects/hippocampus_breathhold/FLASH/derivatives/hippunfold/' subid '/surf'];

subfs = [mypath '/sub-' subid '_hemi-L_space-T2w_den-0p5mm_label-hipp_atlas-bigbrain_subfields.label.gii'];
mthick = [mypath '/sub-' subid '_hemi-L_space-unfolded_den-0p5mm_label-hipp_midthickness.surf.gii'];

for SURF = 1:2
    VPF_plot_hippocampus_unfolded(dS_mean_echo(:,SURF),mthick,subfs,'fireice',[-100 100],0);
    gaxis = gca;
    VPF_rot_hippocampus_flatmap(gaxis);
    title(['\DeltaS_{mean} ' SURFS{SURF}])
    g = gifti(dS_mean_echo(:,SURF));
    save(g,[func_path '/sub-' subid '_hemi-avg_' SURFS{SURF} '_dSbreathhold.shape.gii'],'Base64Binary');    
end

for SURF = 1:2
    VPF_plot_hippocampus_unfolded(dS_mean_echo_weighted(:,SURF),mthick,subfs,'fireice',[-18 18],0);
    gaxis = gca;
    VPF_rot_hippocampus_flatmap(gaxis);
    title(['\DeltaS_{mean,weighted} ' SURFS{SURF}])
    g = gifti(dS_mean_echo_weighted(:,SURF));
    save(g,[func_path '/sub-' subid '_hemi-avg_' SURFS{SURF} '_dSbreathhold_weighted.shape.gii'],'Base64Binary');    
end

for SURF = 1:2
    VPF_plot_hippocampus_unfolded(T2s(:,SURF),mthick,subfs,'inferno',[0 80],0);
    gaxis = gca;
    VPF_rot_hippocampus_flatmap(gaxis);
    title(['T2s ' SURFS{SURF}])
    saveas(gcf,[func_path '/T2s_from_vertices_' SURFS{SURF} '.fig']);
    g = gifti(T2s(:,SURF));
    save(g,[func_path '/sub-' subid '_hemi-avg_' SURFS{SURF} '_T2s.shape.gii'],'Base64Binary');    
end





