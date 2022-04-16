function [Rmean, Rmean_std, Rstd, isodiv, isodiv_std, eapkurt, eapkurt_std] = ...
    extract_DSTATS (DSTATS)
    
% 
% RJ- wrapper to save fields in DSTATS as workspace vars
% 

Rmean = DSTATS.Rmean;
Rmean_std = DSTATS.Rmean_std;
Rstd = DSTATS.Rstd;
isodiv = DSTATS.isodiv;
isodiv_std = DSTATS.isodiv_std;
eapkurt = DSTATS.eapkurt;
eapkurt_std = DSTATS.eapkurt_std;

end
