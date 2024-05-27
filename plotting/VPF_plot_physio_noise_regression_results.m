clear;clc;close all;
%plots compcor vs no compcor results for testing
subjects = {'7495','7566'};

mainpath = '/media/pfaffenrot/Elements/postdoc/projects/data';

colorcode = VPF_create_hippocampus_colorcode();
plotspecs = struct('FontName','Arial','FontSize',22,'colormap',jet(256),...
    'xtick',[],'xlim',[1 30],'LineWidth',2);
plotspecs.color = colorcode(:,1);

plotspecs.ytick = -0.5:0.5:2.5;
plotspecs.ylim = [-0.5 2.5];

for subject = 1:length(subjects)
    for session = 1:2
        load([mainpath '/' subjects{subject} '/memory/ses-0' num2str(session) '/native/'...
            'results_memory_vs_math_vessel_masked_nocompcor.mat']);

        if subject ==1 && session == 1
            con_array = zeros([size(results_memory_vs_math.con_array),2,2,2]); %wo/w compcor,ses1/ses2,subject1/2
        end
        con_array(:,:,1,session,subject) = results_memory_vs_math.con_array;

        load([mainpath '/' subjects{subject} '/memory/ses-0' num2str(session) '/native/'...
            'results_memory_vs_math_vessel_masked.mat']);
        con_array(:,:,2,session,subject) = results_memory_vs_math.con_array;
    end
end


N_layers = size(con_array,2);
N_compcor = size(con_array,3);
N_subfields = size(con_array,1);
N_sessions = size(con_array,4);
N_subjects = size(con_array,5);
subfield_label = {'Subiculum','CA1','CA2','CA3','CA4/DG'};
%% left: no compcor, right: compcor
%wo/w compcor,ses1/ses2,subject1/2

%As I am just interest whether the direction of layers is 
%more similar with compcor, I calculate the rmse of the gradient
%across layers

RMSE_Grad = zeros(N_subfields,N_compcor,N_subjects);
for subfield = 1:N_subfields
    for docompcor = 1:N_compcor
        for subject = 1:N_subjects
            layers = squeeze(con_array(subfield,10:30,docompcor,:,subject));
            dlayers = [gradient(layers(:,1)) gradient(layers(:,2))];
         
            RMSE_Grad(subfield,docompcor,subject) = rmse(dlayers(:,1),dlayers(:,2));
        end
    end
end
rel_RMSE = squeeze(((RMSE_Grad(:,2,:)-RMSE_Grad(:,1,:))));


figure, 
b = bar(rel_RMSE,'LineWidth',2);

b(1).FaceColor = 'flat';
b(2).FaceColor = 'flat';

for ii = 1:5
    b(1).CData(ii,:) = colorcode{1,1};
    b(2).CData(ii,:) = colorcode{2,1};
end
g = gca;
set(g,'xticklabel',{'Sub','CA1','CA2','CA3','DG/CA4'})
set(g,'FontName','Arial','FontSize',24,'YLim',[-0.6 0.6])
ylabel('\DeltaRMSE_{gradient} [a. u.]')
title('Effect of CompCor')
legend({'suject 1','subject 2'})

%%
figure,
subplot(1,2,1)
plot(con_array(:,:,1,1,1).'),ylim([-10 30]),title('sub1 session1 no compcor')
legend(subfield_label)
subplot(1,2,2)
plot(con_array(:,:,1,2,1).'),ylim([-10 30]),title('sub1 session2 no compcor')
legend(subfield_label)

figure,
subplot(1,2,1)
plot(con_array(:,:,2,1,1).'),ylim([-10 30]),title('sub1 session1 compcor')
legend(subfield_label)
subplot(1,2,2)
plot(con_array(:,:,2,2,1).'),ylim([-10 30]),title('sub1 session2 compcor')
legend(subfield_label)
