clear;clc;

clear;clc;

data = dir('/media/pfaffenrot/Elements/postdoc/projects/data/7*');


nScan = length(data);
alph_FWE = 0.05;
bNeg = 0;
for Scan = 1:nScan
    if strcmp(data(Scan).name,'7495') ||  strcmp(data(Scan).name,'7566')
        for session = 1:2
            load([data(Scan).folder '/' data(Scan).name '/memory/ses-0' num2str(session)...
                '/z_transformed/results_memory_vs_math.mat']);
            if session == 1
                tmp = zeros(size(results_memory_vs_math.con_array));
                if Scan == 1
                    X = zeros([size(results_memory_vs_math.con_array),nScan]);
                end
            end
            tmp = tmp+results_memory_vs_math.con_array;
            X(:,:,Scan) = tmp./2;
        end
    else
        load([data(Scan).folder '/' data(Scan).name...
            '/memory/z_transformed/results_memory_vs_math.mat']);
        if Scan == 1
            X = zeros([size(results_memory_vs_math.con_array),nScan]);
        end
        X(:,:,Scan) = results_memory_vs_math.con_array;
    end
end

X = X(:,:,[1 2 3 4 5 6 7 8 9]);

%%
p = zeros(30,1);
h = p;
for ii = 10:30
    [p(ii),h(ii)] = signrank(squeeze(X(2,ii,:)),squeeze(X(4,ii,:)),'alpha',0.05);
end

%%
m = squeeze(mean(X([2,4],10:30,:),3)).';
s = squeeze(std(X([2,4],10:30,:),[],3)).';
s = s./sqrt(size(X,3));

figure,
bar(m)
hold on
for ii = 1:2
    if ii == 1
        x = (1:21) -0.15;
    else
        x = (1:21) + 0.15 * (ii- 1);
    end
    er = errorbar(x,squeeze(m(:,ii)),squeeze(s(:,ii)),'LineWidth', 2);
    er(1).Color = [0 0 0];
    er(1).LineStyle = 'none';
end




