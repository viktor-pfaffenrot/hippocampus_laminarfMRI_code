function lme = VPF_LME_slope_change_vs_density_change(x,y,w)
%x = change in vessel density
%y = slope of signal change curve
%w = weights
if nargin < 3
    w = ones(size(x));
end
N_subjects = size(x,1);
subfields = {'sub','CA1','CA2','CA3','CA4/DG'}.';
subjects = {'subject 1','subject 2','subject 3','subject 4','subject 5'...
    'subject 6','subject 7','subject 8'};
N_subfields = length(subfields);

% subject_label = {};
% for subject = 1:N_subjects
%     subject_label = cat(1,subject_label, repmat({['subject ' num2str(subject)]},[size(x,2) 1]));
% end
% 
% subfield_label = repmat(subfields,[N_subjects,1]);

subject_label = repmat(subjects,[1,N_subfields]).';

subfield_label = {};
for subfield = 1:N_subfields
    subfield_label = cat(1,subfield_label,repmat(subfields(subfield),[size(x,1) 1]));
end
x = reshape(x,[N_subjects*N_subfields 1]);
y = reshape(y,[N_subjects*N_subfields 1]);

tbl = table(subject_label,subfield_label,x,y,'VariableNames',...
    {'subject','subfield','density_change','slope_change'});

tbl.subject = categorical(tbl.subject);
tbl.subfield = categorical(tbl.subfield);
tbl(tbl.density_change==0,:) = [];

% tbl = tbl(1:5:end,:);
% tbl = tbl(:,3:4);

% formula = 'slope_change ~ 1 + density_change + (density_change|subject)';% + (density_change|subfield)';
formula = 'slope_change ~ density_change + (1|subject) + (1|subfield)';

lme = fitlme(tbl,formula,'FitMethod','REML');


keyboard
% lme = fitlm(tbl)
% keyboard
% ypred = predict(lme);
anova(lme,'DFMethod','satterthwaite')
R_lme = lme.Rsquared
% coeff = double(lme.Coefficients(:,2));
% p = double(lme.Coefficients(:,6));
% CI = [double(lme.Coefficients(:,7)),double(lme.Coefficients(:,8))];
end