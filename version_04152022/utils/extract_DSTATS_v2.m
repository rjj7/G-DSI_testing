function [ Rmean, Rmean_std, Rstd, GaussDiv, GaussDiv_std, KurtProfile, ...
    KurtProfile_std, MeanKurt, MeanKurt_std   ] = ...
    extract_DSTATS_v2( DSTATS )
    
% 
% RJ- wrapper to save fields in DSTATS as workspace vars
% 

Rmean = DSTATS.Rmean;
Rmean_std = DSTATS.Rmean_std;
Rstd = DSTATS.Rstd;
GaussDiv = DSTATS.GaussDiv;
GaussDiv_std = DSTATS.GaussDiv_std;
KurtProfile = DSTATS.KurtProfile;
KurtProfile_std = DSTATS.KurtProfile_std;
MeanKurt = DSTATS.MeanKurt;
MeanKurt_std = DSTATS.MeanKurt_std;

end
