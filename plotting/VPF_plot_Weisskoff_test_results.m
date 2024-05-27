clear;clc;close all;
load('/media/pfaffenrot/Elements/postdoc/projects/data/Weisskoff_test_result_all_persubfield.mat');

%%
colorcode = VPF_create_hippocampus_colorcode();
subfield_label = {'subiculum','CA1','CA2','CA3','CA4/DG'};
myclim = [0.8 2];
colVec = linspace(myclim(1), myclim(2), 128);
MM = cubehelix(256,1,-2,3,0.5);

N_subjects = size(w_test,1);
N_subfields = 5;
N_layers = 30;
mymax = zeros(N_subfields,1);
for subfield = 1:N_subfields
    mymax(subfield) = size(w_test{1}{subfield},1);
    for ii = 1:N_subjects
        tmp = w_test{ii}{subfield};
        tmp = size(tmp,1);
        if tmp>mymax(subfield)
            mymax(subfield) = tmp;
        end
    end
end
mymin = [3324,7139,1401,1979,4257];
subfield_cell = cell(N_subfields,1);

for subfield = 1:N_subfields
    subfield_cell{subfield} = zeros([mymax(subfield) N_layers]);

    tmp = 0;
    for ii = 1:N_subjects
        tmp = w_test{ii}{subfield};
        idx = size(tmp,1);
        if subfield == 1 || subfield == 4
            subfield_cell{subfield}(1:idx,10:end) = subfield_cell{subfield}(1:idx,10:end) + tmp(:,10:end);
        else
            subfield_cell{subfield}(1:idx,:) = subfield_cell{subfield}(1:idx,:) + tmp;
        end
    end
    subfield_cell{subfield} = subfield_cell{subfield}./N_subjects;
end

subfield_cell{1}(:,1:9) = nan;
subfield_cell{4}(:,1:9) = nan;



mdata_lineplot = zeros(N_subfields,N_layers);
for subfield = 1:N_subfields
    layers = 1:N_layers;
    N_samples = mymax(subfield);
    samples = sqrt(1:N_samples);
    
    mdata = log10(subfield_cell{subfield});

    s = 1:floor(mymin(subfield)/10)*10;
    mymin(subfield) = length(s);

    [X,Y] = meshgrid(layers,samples(s));
    mdata_show = mdata(1:mymin(subfield),:);

    [xk,yk] = meshgrid(-5:5);
    K = exp(-(xk.^2 + yk.^2)/10); K = K/sum(K(:));

    plotspecs = struct('FontName','Arial','FontSize',22,'colormap',jet(256));

    if subfield == 1 || subfield == 4 || subfield == 5
        mdata_show = mdata_show(:,10:end);
        mdata_show = conv2(mdata_show,K,'same')./conv2(ones(size(mdata_show)),K,'same');
        mdata_show2 = zeros([size(mdata_show,1) 30]);
        mdata_show2(:,10:end) = mdata_show;
        mdata_show = mdata_show2;
        clear mdata_show2;
    else
        mdata_show = conv2(mdata_show,K,'same')./conv2(ones(size(mdata_show)),K,'same');
    end


    mdata_lineplot(subfield,:) = mdata_show(end,:);

    figure,
    surf(X,Y,mdata_show,'edgecolor','k');colormap(MM);
    ylabel('$\sqrt{N}$','FontName',plotspecs.FontName,'FontSize',plotspecs.FontSize,Interpreter='latex'),
    zlabel('temporal noise','FontName',plotspecs.FontName,'FontSize',plotspecs.FontSize),
    title(subfield_label{subfield})
    shading interp
    hold on
    line(X,Y(end,1)*ones(size(Y)),mdata_show(end,:),'Color',colorcode{subfield,1},'LineWidth',7)
    line(10*ones(size(X)),Y(:,10),mdata_show(:,10),'Color','k','LineWidth',5)

    hidden off
    view(-156,21);
    set(gca,'ZScale','log')
    set(gca,'ZLim',[0 1.7])
    set(gca,'CLim',myclim)
    set(gca,'Xdir','reverse')
    set(gca,'XTick',[])
    set(gca,"FontSize",plotspecs.FontSize)

    [~,idx] = max(mdata_show(end,:));
    max_noise_line = mdata_show(:,idx);
    dmax_noise_line = diff(max_noise_line);
    smoothed_dmax = smoothdata(dmax_noise_line);
    pos = floor(mean(find(abs(smoothed_dmax)<1e-7)))
    figure, plot(smoothed_dmax),hold on, line(1:length(dmax_noise_line),zeros(length(dmax_noise_line),1))
    title(subfield_label{subfield})

    limit_reached = max_noise_line(pos);

    figure, plot(max_noise_line),hold on, line(1:length(dmax_noise_line),ones(length(dmax_noise_line))*limit_reached);
    title(subfield_label{subfield})

end
mdata_lineplot( mdata_lineplot==0) = nan;
mdata_lineplot(end,1:9) = nan;
%%
plotspecs = struct('FontName','Arial','FontSize',22,'colormap',jet(256),...
    'xtick',[],'xlim',[1 30],'LineWidth',2);
plotspecs.color = colorcode(:,1);

plotspecs.ylim = [0.9 1.4];
plotspecs.ytick = 0.9:0.1:1.4;

figure,
h = VPF_show(@semilogy,1:30,mdata_lineplot,[],'Weisskoff Test N = 9',[],'temporal noise',plotspecs);
for kk = 1:length(h)
    set(h(kk), 'Color',plotspecs.color{kk});
end
hold on
plotspecs.Marker = '.';
plotspecs.MarkerSize = 42;
h = VPF_show(@semilogy,10,mdata_lineplot(end,10),[],[],[],[],plotspecs);
set(h, 'Color',plotspecs.color{end});
h = VPF_show(@semilogy,30,mdata_lineplot(end,end),[],'Weisskoff Test N = 9',[],'temporal noise',plotspecs);
set(h, 'Color',plotspecs.color{end});
hold off
l = line([10,10],plotspecs.ylim,'LineWidth',plotspecs.LineWidth,'Color','black');
legend(subfield_label,'Location','northwest');

experiment_mean = mean(mdata_lineplot,2,'omitnan');

india_ink_path = '/home/pfaffenrot/work/my_publics/hippocampus_layers/images/F_8/india_ink_sampling';
VPF_sample_india_ink_stain(experiment_mean,india_ink_path);




