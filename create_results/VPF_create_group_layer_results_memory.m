clear;clc;

data = dir('/media/pfaffenrot/Elements/postdoc/projects/data/7*');
subjects = length(data);
for subject = 1:subjects
    if strcmp(data(subject).name,'7495') ||  strcmp(data(subject).name,'7566')
        for session = 1:2
            load([data(subject).folder '/' data(subject).name '/memory/ses-0' num2str(session)...
                '/z_transformed/results_memory_vs_math_z.mat']);
            if session == 1
                tmp = zeros(size(results_memory_vs_math.con_array));
                if subject == 1
                    X = zeros([size(results_memory_vs_math.con_array),subjects]);
                end
            end
            tmp = tmp+results_memory_vs_math.con_array;
            X(:,:,subject) = tmp./2;
        end
    else
        load([data(subject).folder '/' data(subject).name...
            '/memory/z_transformed/results_memory_vs_math_z.mat']);
        if subject == 1
            X = zeros([size(results_memory_vs_math.con_array),subjects]);
        end
        X(:,:,subject) = results_memory_vs_math.con_array;
    end
end

X = X(:,:,[1 2 3 4 5 6 7 8 9]);

Xm = mean(X,3);
Xs = std(X,[],3)./sqrt(size(X,3));
outpath = '/media/pfaffenrot/Elements/postdoc/projects/data/avg/memory';
outname = 'memory_vs_math_z_all_layers_no_masked';
memory_vs_math = reshape(permute(X,[1 3 2]),[45 30]);

save([outpath '/' outname '.mat'],'memory_vs_math')
%% create tSNR 
clear h
subfields = 5;
pos = 30;
tSNR = zeros(subfields,pos,subjects);
for subject = 1:subjects
    if strcmp(data(subject).name,'7495') || strcmp(data(subject).name,'7566') %subjects with 2 sessions
        tmp = zeros(subfields,pos);
        for session = 1:2
            load([data(subject).folder '/' data(subject).name '/memory/ses-0' num2str(session) '/z_transformed/results_memory_vs_math_z.mat']);
            tmp = tmp + results_memory_vs_math.tSNR;
        end
        tmp = tmp./2;
        tSNR(:,:,subject) = tmp;
    else
        load([data(subject).folder '/' data(subject).name '/memory/z_transformed/results_memory_vs_math_z.mat']);
        tSNR(:,:,subject) = results_memory_vs_math.tSNR;
    end

end
m_tSNR = mean(tSNR,3).';
s_tSNR = (std(tSNR,[],3)./sqrt(subjects)).';


results = struct('con_array',Xm,'tSNR',m_tSNR.');

tSNR_outname = 'memory_vs_math_math_tSNR_all_layers_no_masked';
save([outpath '/' tSNR_outname '.mat'],'m_tSNR','s_tSNR');

results_json = jsonencode(results,PrettyPrint=true);
fid = fopen([outpath '/' outname '_math_tSNR.json'],'w');
fprintf(fid,'%s',results_json);
fclose(fid);
