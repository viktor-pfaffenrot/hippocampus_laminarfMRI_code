clear;clc;

data = dir('/media/pfaffenrot/Elements/postdoc/projects/data/7*');


nScan = length(data);
for Scan = 1:nScan
    if strcmp(data(Scan).name,'7495') ||  strcmp(data(Scan).name,'7566')
        for session = 1:2
            load([data(Scan).folder '/' data(Scan).name '/memory/ses-0' num2str(session)...
                '/z_transformed/results_pre_vs_post_z_vessel_masked.mat']);
            if session == 1
                tmp = zeros(size(results_pre_vs_post.con_array));
                if Scan == 1
                    X = zeros([size(results_pre_vs_post.con_array),nScan]);
                end
            end
            tmp = tmp+results_pre_vs_post.con_array;
            X(:,:,Scan) = tmp./2;
        end
    else
        load([data(Scan).folder '/' data(Scan).name...
            '/memory/z_transformed/results_pre_vs_post_z_vessel_masked.mat']);
        if Scan == 1
            X = zeros([size(results_pre_vs_post.con_array),nScan]);
        end
        X(:,:,Scan) = results_pre_vs_post.con_array;
    end
end

X = X(:,:,[1 2 3 4 5 6 7 8 9]);
sz = size(X);
%%
% X = mean(reshape(X,[prod(sz(1:2)),sz(3)]),1,'omitnan');
X = squeeze(mean(X(4,20:end,:),2,'omitnan'));
Xm = mean(X);
Xs = std(X);

save('avg_4_power_analysis_pre_vs_post_onlyCA3_outer.mat','Xm','Xs')