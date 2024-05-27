function VPF_permutation_t_test_memory_layers_SnPM(X,alpha_FWE,bNeg,LAYERS)

if nargin < 4
    LAYERS = true;
end

if nargin < 3
    bNeg = 0;
end

if nargin < 2 
    alpha_FWE = 0.05;
end
%%
if ~LAYERS
    smoothingFactor = 0.75;
    maxiter = 1;
    labels = dir([data(1).folder '/' data(1).name '/hippunfold/surf/sub-*_hemi-L_space-T2w_den-0p5mm_label-hipp_atlas-bigbrain_subfields.label.gii']);
    labels = gifti([labels.folder '/' labels.name]).cdata;
    for iter = 1:maxiter
        for ii = 1:size(X,3)
            for jj = 1:max(labels)
                for SURF = 1:size(X,2)
                    inp = X(labels==jj,SURF,ii);
                    inp(inp==0)= nan;
                    inp = smoothdata(inp,1,'gaussian','smoothingFactor',smoothingFactor);
                    X(labels==jj,SURF,ii) = inp;
                end
            end
        end
    end
end
%%
nScan = size(X,3);
iCond = ones(nScan,1);

if snpm_get_defaults('shuffle_seed')
    % Shuffle seed of random number generator
    try
        rng('shuffle');
    catch
        % Old syntax
        rand('seed',sum(100*clock));
    end
end



PiCond=[];
for i=0:nScan-1
    PiCond=[ones(2^i,1),PiCond;-ones(2^i,1),PiCond];
end

%-Only do half the work, if possible
PiCond=PiCond(PiCond(:,1)==1,:);
nPerm   = size(PiCond,1);		%-# permutations
MaxT  = repmat(-Inf,nPerm,2);	%-Max t

perm = find(all((meshgrid(iCond,1:size(PiCond,1))==PiCond)'));


if (perm<0), PiCond=-PiCond; perm=-perm; end
%-Actual labelling must be at top of PiCond
if (perm~=1)
    PiCond(perm,:)=[];
    PiCond=[iCond;PiCond];
end

%-Randomise order of PiConds, unless already randomized
% Allows interim analysis
PiCond=[PiCond(1,:);PiCond(randperm(size(PiCond,1)-1)+1,:)];


sHCform    = 'spm_DesMtx(PiCond(perm,:),''C'',''Mean'')';
%-Condition partition
[H,Hnames] = spm_DesMtx(iCond,'C','Mean');
%-Contrast of condition effects
% (spm_DesMtx puts condition effects in index order)
CONT       = 1;
%-No block/constant
B=[]; Bnames='';

if size(H,2)>1
    CONT(1:size(H,2)) = CONT(1:size(H,2)) - mean(CONT(1:size(H,2)));
end

[nHCBG,HCBGnames] = spm_DesMtx('Sca',H,Hnames,[],[],[],[],[],[]);


%%

X_in = X;
% for subfield = 1:size(X,1)

X = squeeze(X_in(:,:,:));
%     if subfield ==5
%         X = X([10,30],:);
%     end
sz = size(X);

STalpha = snpm_get_defaults('STalpha');
STprop  = snpm_get_defaults('STprop');

q       = size(H,1);		%-# observations
p       = size(H,2);		%-# predictors
r       = rank(H);		%-Model degrees of freedom
df      = q - r;			%-Residual degrees of freedom

X = reshape(X,[prod(sz(1:end-1)) sz(end)]).';

[~,nonanidx] = find(~isnan(X(2,:)));
X = X(:,nonanidx);

perm = 1;
BETA  = pinv(H)*X;
ResSS = sum((X - H*BETA).^2);

T      = zeros(1,size(BETA,2));
T(1,:) = CONT*BETA./sqrt((ResSS*(CONT*pinv(H'*H)*CONT'))/df);

MaxT(perm,:) = max([ max(T(1,:)), -min(T(1,:));   ...
    MaxT(perm,1),  MaxT(perm,2) ]);

T0 = T;
SnPMt=T;
nP = [];
nPtmp = ones(size(T));
nPtmp = nPtmp + (T0<=0);

for perm = 2:nPerm

    %-Rebuild H C for current permuation
    %-----------------------------------------------------------
    %       HC = eval(sHCform);
    HC = spm_DesMtx(PiCond(perm,:),'C','Mean');
    BETA  = pinv(HC)*X;
    ResSS = sum((X - HC*BETA).^2);
    T      = zeros(1,size(BETA,2));
    % t, as usual
    T(1,:) = CONT*BETA./sqrt((ResSS*(CONT*pinv(HC'*HC)*CONT'))/df);

    MaxT(perm,:) = max([ max(T(1,:)), -min(T(1,:));      ...
        MaxT(perm,1), MaxT(perm,2) ]);

    %-Update nonparametric P-value
    %-----------------------------------------------------------

    nPtmp = nPtmp + (T>=T0) + (-T>=T0);
end
nP = [nP, nPtmp];


nP = nP/(2*nPerm);

tol = 1e-4;	% Tolerance for comparing real numbers
cP_pos=zeros(size(nP));
MaxT_pos=MaxT(:,1);
for t = MaxT_pos'
    %-FEW-corrected p is proportion of randomisation greater or
    % equal to statistic.
    %-Use a > b -tol rather than a >= b to avoid comparing
    % two reals for equality.
    cP_pos = cP_pos + (t > SnPMt -tol);
end
cP_neg=zeros(size(nP));
MaxT_neg=MaxT(:,2);

for t = MaxT_neg'
    cP_neg = cP_neg + (t > -SnPMt -tol);
end

cP_pos = cP_pos / (2* nPerm);
cP_neg = cP_neg / (2* nPerm);


lP_FWE_pos=-log10(cP_pos);
lP_FWE_neg=-log10(cP_neg);


nPermReal = size(MaxT,1); %different with nPerm when bhPerms==1
if bNeg
    SnPMt    = -SnPMt;
    nP = 1+1/nPermReal-nP;
end

[StMaxT, iStMaxT] = sort(MaxT);

iFWE=ceil((1-alpha_FWE)*nPermReal);
Tcrit=StMaxT(iFWE);

if LAYERS==true
    tmp = nan(prod(sz(1,2)),1);
    tmp(nonanidx) = SnPMt;
    T_out = reshape(tmp,sz(1:2));
    T_out(T_out==0) = nan;

    lindex = linspace(1, 2, 30-9);
    T_out(end,:) = cat(2,nan([1,9]),interp1(1:2, T_out(end,[10 30]), lindex, 'linear', 0));
else
    T_out = reshape(SnPMt,sz(1:2));
end

% end
T = T_out;
results = struct('T',T,'Tcrit',Tcrit);
results_json = jsonencode(results,PrettyPrint=true);

results_memory_vs_math_memory_vs_math = results;
outname = 'memory_vs_math';
outpath = '/media/pfaffenrot/Elements/postdoc/projects/data/avg/memory';

if LAYERS==true
    fulloutname = [outpath '/results_' outname '_vessel_masked.mat'];
    fid = fopen([outpath '/results_' outname '_vessel_masked.json'],'w');
else
    fulloutname = [outpath '/results_' outname '_unfolded_vessel_masked.mat'];
    fid = fopen([outpath '/results_' outname '_unfolded_vessel_masked.json'],'w');
end


% save(fulloutname,['results_' outname])
% fprintf(fid,'%s',results_json);
% fclose(fid);


if LAYERS == true
    colorcode = VPF_create_hippocampus_colorcode();
    titles = [cellstr('Subiculum'), cellstr('CA1'), cellstr('CA2'),...
        cellstr('CA3'),cellstr('CA4/DG')];

    plotspecs = struct('FontName','Arial','FontSize',22,'colormap',jet(256),...
        'xtick',[],'xlim',[1 30],'LineWidth',2);
    plotspecs.color = colorcode(:,1);

    fields_to_plot = length(titles);

    plotspecs.ytick = 0:2:14;
    plotspecs.ylim = [0 14];
    %
    %     plotspecs.ytick = -4:2:6;
    %     plotspecs.ylim = [-4 6];

    figure,
    h = VPF_show(@plot,1:30,T(1:fields_to_plot,:).',[],[],[],'T [a.u.]',plotspecs);
    for kk = 1:length(h)
        set(h(kk), 'Color',plotspecs.color{kk});
    end
    lT = line(plotspecs.xlim,[Tcrit(1),Tcrit(1)],'LineWidth',plotspecs.LineWidth,'LineStyle','--','Color','black');
    l = line([10,10],plotspecs.ylim,'LineWidth',plotspecs.LineWidth,'Color','black');
    legend(titles);
    pos = lT(1).Parent.Position;

    %     title('memory > math')
else
    subfs_name = dir([data(1).folder '/' data(1).name '/hippunfold/surf/sub-*_hemi-L_space-T2w_den-0p5mm_label-hipp_atlas-bigbrain_subfields.label.gii']);
    mthick = dir('/media/pfaffenrot/Elements/postdoc/projects/data/avg/hippunfold/canonical_hemi-avg_inner.unfolded.surf.gii');
    SURFS = {'inner','midthickness','outer'};
    for SURF = 1:length(SURFS)

        VPF_plot_hippocampus_unfolded(T_out(:,SURF),[mthick.folder '/' mthick.name],[subfs_name.folder '/' subfs_name.name],'hot',[Tcrit Tcrit+2],0);
        colorbar;
        g = gca;
        VPF_rot_hippocampus_flatmap(g);
        title([SURFS{SURF}])
    end

end
end




