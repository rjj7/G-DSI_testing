function [ Rmean, Rmean_std, Rstd, Rstd_std, GaussDiv, GaussDiv_std, KurtProfile, ...
    KurtProfile_std, MeanKurt, MeanKurt_std, MeanEAP, MeanEAP_std, ...
    MeanGA, MeanGA_std, RTOP, RTOP_std, sintheta, sintheta_std, ...
    MSD, MSD_std, MFD, MFD_std, HWHM, HWHM_std, HWHMmin, HWHMmin_std, ...
    HWHMmax, HWHMmax_std, HWHMminmax, HWHMminmax_std, HWHMstd, HWHMstd_std ] = ...
    extract_DSTATS_v3( DSTATS )
    
% 
% RJ- wrapper to save fields in DSTATS as workspace vars
% 

Rmean = DSTATS.Rmean;
Rmean_std = DSTATS.Rmean_std;

Rstd = DSTATS.Rstd;
Rstd_std = DSTATS.Rstd_std;

GaussDiv = DSTATS.GaussDiv;
GaussDiv_std = DSTATS.GaussDiv_std;

KurtProfile = DSTATS.KurtProfile;
KurtProfile_std = DSTATS.KurtProfile_std;

MeanKurt = DSTATS.MeanKurt;
MeanKurt_std = DSTATS.MeanKurt_std;

MeanEAP = DSTATS.MeanEAP;
MeanEAP_std = DSTATS.MeanEAP_std;

MeanGA = DSTATS.MeanGA;
MeanGA_std = DSTATS.MeanGA_std;

RTOP = DSTATS.RTOP;
RTOP_std = DSTATS.RTOP_std;

sintheta = DSTATS.sintheta;
sintheta_std = DSTATS.sintheta_std;

MSD = DSTATS.MSD;
MSD_std = DSTATS.MSD_std;

% MFD = DSTATS.MFD;
% MFD_std = DSTATS.MFD_std;
MFD = [];
MFD_std = [];

HWHM = DSTATS.HWHM;
HWHM_std = DSTATS.HWHM_std;

HWHMmin = DSTATS.HWHMmin;
HWHMmin_std = DSTATS.HWHMmin_std;

HWHMmax = DSTATS.HWHMmax;
HWHMmax_std = DSTATS.HWHMmax_std;

HWHMminmax = DSTATS.HWHMminmax;
HWHMminmax_std = DSTATS.HWHMminmax_std;

HWHMstd = DSTATS.HWHMstd;
HWHMstd_std = DSTATS.HWHMstd_std;

end
