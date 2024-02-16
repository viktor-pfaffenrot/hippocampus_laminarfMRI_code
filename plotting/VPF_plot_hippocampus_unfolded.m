function VPF_plot_hippocampus_unfolded(data,surfpath,subfs_path,cmap,climinp,smooth)



figure;
sluh = gifti(surfpath);
sluh.vertices(:,[1 2 3]) = sluh.vertices(:,[2 1 3]);
sluh.vertices(:,1) = -sluh.vertices(:,1); % flip left
g = gca;
sluh.vertices(:,1) = sluh.vertices(:,1) + (g.XLim(1)-max(sluh.vertices(:,1))) -1; % translate to left edge
sluh.vertices(:,2) = sluh.vertices(:,2) -mean(sluh.vertices(:,2)); % align middle

plot_gifti(sluh,data,false,smooth);
cmap = eval([cmap '(256)']);
colormap(gca,cmap)

hold on

subfs = gifti(subfs_path);

plot_giftiborders(sluh,subfs.cdata);
clim(climinp);

axis image

end
