clear;clc;%close all;

%just a small helper to see MRI metrics on the folder surface. 'surf' is 
%either outer.native or inner.native. 'g' is the metric .shape.gii, e.g.
%angio, dSbreathhold, or T2s. Easiest way to see these plots is to go to 
%the public shiny app (data can also be downloaded from there:
%https://viktor-pfaffenrot.shinyapps.io/hippocampus_data_viewer/

data = dir('/media/pfaffenrot/Elements/postdoc/projects/data/7*');

subjects = length(data);

N_vertices = 7262;

%%
g = gifti('/media/pfaffenrot/Elements/postdoc/projects/data/avg/breathhold/canonical_hemi-avg_inner_angio.shape.gii');

dat = g.cdata;

surf = '/media/pfaffenrot/Elements/postdoc/projects/data/avg/hippunfold/canonical_hemi-L_outer.native.surf.gii';
sub_fs = '/media/pfaffenrot/Elements/postdoc/projects/data/7218/hippunfold/surf/sub-7218_hemi-L_space-T2w_den-0p5mm_label-hipp_atlas-bigbrain_subfields.label.gii';

VPF_plot_hippocampus_unfolded(dat,surf,sub_fs,'fireice',[-10 10],2);
view(-12,45)