clear;clc;%close all;


%just a small helper to see MRI metrics on the folder surface. 'surf' is 
%either outer.native or inner.native. 'g' is the metric .shape.gii, e.g.
%angio, dSbreathhold, or T2s. Easiest way to see these plots is to go to 
%the public shiny app (data can also be downloaded from there:
%https://viktor-pfaffenrot.shinyapps.io/hippocampus_data_viewer/

N_vertices = 7262;

%%
colorcode = VPF_create_hippocampus_colorcode();
colorcode = cell2mat(colorcode(:,1));
g = gifti('//media/pfaffenrot/PostDoc_data/projects/data/avg/memory/canonical_hemi-avg_outer_tSNR_vessel_masked.shape.gii');

dat = g.cdata;

surf = '/media/pfaffenrot/PostDoc_data/projects/data/avg/hippunfold/canonical_hemi-L_outer.native.surf.gii';
sub_fs = '/media/pfaffenrot/PostDoc_data/projects/data/avg/hippunfold/canonical_hemi-L_atlas-bigbrain_subfields.label.gii';

VPF_plot_hippocampus_unfolded(dat,surf,sub_fs,'hot',[0 40],1);
% gaxis = gca;
% VPF_rot_hippocampus_flatmap(gaxis);
view(-12,45)