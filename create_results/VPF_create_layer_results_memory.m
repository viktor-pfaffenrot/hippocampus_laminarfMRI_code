function results = VPF_create_layer_results_memory(layers,SPM,con_idx,nonvessel_idx,ZTRANS)

if nargin < 5 || isempty(ZTRANS)
    ZTRANS = false;
end


if nargin < 4 || isempty(nonvessel_idx)
    MASK_VEINS = false;
else
    MASK_VEINS = true;
end



subfields = size(layers,2);
N = size(layers{1,1,1},2);
vols = size(layers{1,1,1},3);
runs = size(layers,3);
Tvols = vols*runs;


W = SPM.xX.W;
GLM_use = SPM.xX.pKX;

con = SPM.xCon(con_idx).c;
T = zeros(subfields-1,N);
Tcrit = zeros(subfields-1,1);
pmax = Tcrit;
con_array = T;
tSNR = T;
for ss = 1:subfields-1
    if ss==subfields-1
        Y = zeros(2,vols,runs);
    else
        Y = zeros(N,vols,runs);
    end
    if ss < subfields-1%all but CA4/DG
        for run = 1:runs
            if MASK_VEINS
                Y(:,:,run) = squeeze(mean(layers{1,ss,run}(nonvessel_idx{1,ss},:,:),1))...
                    + squeeze(mean(layers{2,ss,run}(nonvessel_idx{2,ss},:,:),1));
            else
                Y(:,:,run) = squeeze(mean(layers{1,ss,run},1))...
                    + squeeze(mean(layers{2,ss,run},1));
            end
        end
    else
        % CA4 marks 'outer'. DG marks 'inner'
        % Interpolate final T and beta between inner and outer for visualization purposes
        for run = 1:runs
            %CA4
            if MASK_VEINS
                tmp = (squeeze(mean(mean(layers{1,end-1,run}(nonvessel_idx{1,end-1},:,:),2),1))...
                    + squeeze(mean(mean(layers{2,end-1,run}(nonvessel_idx{2,end-1},:,:),2),1))).';
            else
                tmp = (squeeze(mean(mean(layers{1,end-1,run},2),1))...
                    + squeeze(mean(mean(layers{2,end-1,run},2),1))).';
            end

            %DG
            if MASK_VEINS
                tmp2 = (squeeze(mean(layers{1,end,run}(nonvessel_idx{1,end},:,:),1))...
                    + squeeze(mean(layers{2,end,run}(nonvessel_idx{2,end},:,:),1))).';
            else
                tmp2 = (squeeze(mean(layers{1,end,run},1))...
                    + squeeze(mean(layers{2,end,run},1))).';
            end
            Y(:,:,run) = cat(1,tmp2,tmp);
        end
    end

    %the code above adds the signals of both hemispheres. Divide by 2 to average
    Y = Y./2;

    %calculate tSNR based on first 100 math volumes of the 1st run
    idx = find(SPM.xX.X(SPM.Sess(1).row,SPM.Sess(1).col(3))>0);
    idx = idx(1:100);
    m = mean(Y(:,idx,1),2).';
    s = std(Y(:,idx,1),[],2).';

    if ss == subfields-1
        tmp_tSNR = m./s;
        lindex = linspace(1, size(tmp_tSNR,2), N-9);
        tSNR(ss,:) = cat(2,nan([1,9]),interp1(1:size(tmp_tSNR,2), tmp_tSNR, lindex, 'linear', 0));

    else
        tSNR(ss,:) = m./s;
    end
    tSNR([1 4 5],1:9) = nan;
    
    if ZTRANS
        %basline z-transform input. I take the math condition as baseline
        for run = 1:runs
            idx = find(SPM.xX.X(SPM.Sess(run).row,SPM.Sess(run).col(3))>0);
            m = mean(Y(:,idx,run),2);
            s = std(Y(:,idx,run),[],2);
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

    res      = spm_sp('r',SPM.xX.xKXs,KWY);        %-Residuals
    ResSS    = sum(res.^2);                    %-Residual SSQ
    ResMS = ResSS / SPM.xX.trRV;
    if ss == subfields-1
        [tmpT,Tcrit(ss),tmpcon,pmax(ss)] = VPF_Tmap_from_SPM(SPM,b,ResMS,con,0.05,'FDR');

        lindex = linspace(1, size(tmpT,2), N-9);
        T(ss,:) = cat(2,nan([1,9]),interp1(1:size(tmpT,2), tmpT, lindex, 'linear', 0));
        con_array(ss,:) = cat(2,nan([1,9]),interp1(1:size(tmpcon,2), tmpcon, lindex, 'linear', 0));
    else
        [T(ss,:),Tcrit(ss),con_array(ss,:),pmax(ss)] = VPF_Tmap_from_SPM(SPM,b,ResMS,con,0.05,'FDR');
    end
end

con_array([1 4 5],1:9) = nan;
T([1 4 5],1:9) = nan;

results = struct('con_array',con_array,'T',T,'Tcrit',Tcrit,'pmax',pmax,'tSNR',tSNR);
end