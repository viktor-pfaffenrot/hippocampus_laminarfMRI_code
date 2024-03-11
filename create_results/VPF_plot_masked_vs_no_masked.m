clear;clc; 


load('/media/pfaffenrot/Elements/postdoc/projects/data/avg/memory/pre_vs_post_aggregated_masked.mat');
load('/media/pfaffenrot/Elements/postdoc/projects/data/avg/memory/pre_vs_post_aggregated_no_masked.mat');

masked = permute(reshape(pre_vs_post_aggregated_masked,[5 9 4]),[1 3 2]);
no_masked = permute(reshape(pre_vs_post_aggregated_no_masked,[5 9 4]),[1 3 2]);

masked(end,end,:) = mean(masked(end,[end-1 end],:),2,'omitnan');
masked(end,3,:) = nan;

masked_std = std(masked,[],3,'omitnan')./sqrt(9);
masked_mean = mean(masked,3,'omitnan');

no_masked(end,end,:) = mean(no_masked(end,[end-1 end],:),2,'omitnan');
no_masked(end,3,:) = nan;

no_masked_std = std(no_masked,[],3,'omitnan')./sqrt(9);
no_masked_mean = mean(no_masked,3,'omitnan');

Xm_agg = no_masked_mean-masked_mean;
Xs_agg = no_masked_std-masked_std;



[colorcode,bla] = VPF_create_hippocampus_colorcode();
titles = [cellstr('Subiculum'), cellstr('CA1'), cellstr('CA2'),...
    cellstr('CA3'),cellstr('CA4/DG')];

plotspecs = struct('FontName','Arial','FontSize',22,'colormap',jet(256),...
    'xtick',[],'xlim',[1 30],'LineWidth',2);

plotspecs.ytick = -0.5:0.1:13;
plotspecs.ylim = [-0.5 1.3];

figure,
b = bar(1:5,Xm_agg,'LineWidth',2);
b(1).FaceColor = 'flat';
b(2).FaceColor = 'flat';
b(3).FaceColor = 'flat';
b(4).FaceColor = 'flat';

for ii = 1:5
    b(1).CData(ii,:) = bla{ii,1};
    b(2).CData(ii,:) = colorcode{ii,1};
    b(3).CData(ii,:) = bla{ii,2};
    b(4).CData(ii,:) = colorcode{ii,2};
end

hold on
for ii = 1:4
    if ii == 1
        x = (1:5) -0.27;
    else
        x = (1:5) + 0.18 * (ii-2.5);
    end
    er = errorbar(x,Xm_agg(:,ii),Xs_agg(:,ii),'LineWidth', 2);
    er(1).Color = [0 0 0];
    er(1).LineStyle = 'none';

end

g = gca;
set(g,'xticklabel',{'Sub','CA1','CA2','CA3','DG/CA4'})
set(g,'FontName','Arial','FontSize',24,'YLim',[-0.2 0.3])
ylabel('contrast_{no masking} - contrast_{masking} [\Deltaz]')
title('pre > post')

h2 = bar(nan(4,4));
h2(1).FaceColor = bla{2,1};
h2(2).FaceColor = colorcode{2,1};
h2(3).FaceColor = bla{2,2};
h2(4).FaceColor = colorcode{2,2};
legend(h2,{'SRLM','inner','midthickness','outer'})
