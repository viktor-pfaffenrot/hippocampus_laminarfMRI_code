function [T,Tcrit,con_array,pmax] = VPF_Tmap_from_SPM(SPM,beta,ResMS,con,alpha,flag)

% keyboard
Vc  = con'*SPM.xX.Bcov*con;
SE  = sqrt(ResMS*Vc);    

beta_index = find(abs(con) > 0);
beta_use   = beta(beta_index,:);    

con_array     = zeros(1,size(beta,2));
for j=1:size(beta_use,1)
    con_array = con_array + con(beta_index(j)) * beta_use(j,:);
end


T   = con_array./SE;

switch flag
    case 'FWE'
        Tcrit = spm_uc(alpha,[1 SPM.xX.erdf],'T',SPM.xVol.R,1,numel(beta));
    case 'FDR'
        p = 2 * (1 - spm_Tcdf(T, SPM.xX.erdf));
        p = spm_P_FDR(p);
        Tcrit = min(T(p<alpha));
        if isempty(Tcrit)
            Tcrit = nan;
        end
        pmax = max(p(p<alpha));
        if isempty(pmax)
            pmax = 1;
        end
    case 'none'
        Tcrit = spm_u(alpha,[1 SPM.xX.erdf],'T');
        pmax = alpha;
end


end