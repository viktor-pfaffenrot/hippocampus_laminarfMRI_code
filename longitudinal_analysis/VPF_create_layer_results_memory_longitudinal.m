function out = VPF_create_layer_results_memory_longitudinal(layers,con_idx,idx,SPM,ZTRANS)

if nargin < 6 || isempty(ZTRANS)
    ZTRANS = false;
end


subfields = size(layers,2);
N = size(layers{1,1,1},2);
vols = size(layers{1,1,1},3);
runs = size(layers,3);
Tvols = vols*runs;

W = SPM.xX.W;
GLM_use = SPM.xX.pKX;

con = SPM.xCon(con_idx).c;
con_array = zeros(subfields-1,N);
Y_out = zeros([size(con_array) vols runs]);

for ss = 1:subfields-1
    if ss==subfields-1
        Y = zeros(2,vols,runs);
    else
        Y = zeros(N,vols,runs);
    end

    if ss < subfields-1%all but CA4/DG
        for run = 1:runs
            Y(:,:,run) = squeeze(mean(layers{1,ss,run}(idx{1,ss},:,:),1))...
                + squeeze(mean(layers{2,ss,run}(idx{2,ss},:,:),1));
        end
    else
        % CA4 marks 'outer'. DG marks 'inner'
        % Interpolate final T and beta between inner and outer for visualization purposes
        for run = 1:runs
            %CA4
            tmp = (squeeze(mean(mean(layers{1,end-1,run}(idx{1,end},:,:),2),1))...
                + squeeze(mean(mean(layers{2,end-1,run}(idx{2,end},:,:),2),1))).';
            %DG
            tmp2 = (squeeze(mean(layers{1,end,run},1))...
                + squeeze(mean(layers{2,end,run},1))).';

            Y(:,:,run) = cat(1,tmp2,tmp);
        end
    end

    %the code above adds the signals of both hemispheres. Divide by 2 to average
    Y = Y./2;

    if ZTRANS
        %basline z-transform input. I take the math condition as baseline
        for run = 1:runs
            idx_vol = find(SPM.xX.X(SPM.Sess(run).row,SPM.Sess(run).col(3))>0);
            m = mean(Y(:,idx_vol,run),2);
            s = std(Y(:,idx_vol,run),[],2);
            Y(:,:,run) = (Y(:,:,run)-m)./s;
        end
    end

    if ss==subfields-1
        Y = reshape(Y,[2 Tvols]);
    else
        Y   = reshape(Y,[N Tvols]);
    end
    KWY = spm_filter(SPM.xX.K,W*Y.');
    b   = GLM_use*KWY;

    Y = reshape(Y, [size(Y,1) vols runs]);

    res      = spm_sp('r',SPM.xX.xKXs,KWY);        %-Residuals
    ResSS    = sum(res.^2);                    %-Residual SSQ
    ResMS = ResSS / SPM.xX.trRV;
    if ss == subfields-1
        [~,~,tmpcon,~] = VPF_Tmap_from_SPM(SPM,b,ResMS,con,0.05,'FDR');

        lindex = linspace(1, size(tmpcon,2), N-9);
        con_array(ss,:) = cat(2,nan([1,9]),interp1(1:size(tmpcon,2), tmpcon, lindex, 'linear', 0));
        Y_out(ss,:,:,:) = cat(1,nan([9,size(Y,2),size(Y,3)]),interp1(1:size(tmpcon,2), Y, lindex, 'linear', 0));
    else
        [~,~,con_array(ss,:),~] = VPF_Tmap_from_SPM(SPM,b,ResMS,con,0.05,'FDR');
        Y_out(ss,:,:,:) = Y;
    end
end

con_array([1 4 5],1:9) = nan;
Y = Y_out;
Y([1 4 5],1:9,:,:) = nan;

out = struct('con_array',con_array,'Y',Y);
end