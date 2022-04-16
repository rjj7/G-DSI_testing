function [ DSTATS ] = gdsi_calculate_FS_metrics_Tissue_vHemis( DSTATS, DATAOUT, aparcaseg, fslabels )

% Calculates metrics across different FS tissue labels, save to DSTATS
%  struct
% - mean 1D EAP profile, mean 1D GenAniso profile, etc.
%

%%%% WM
%%% mean 1d eap profile (+ std across labels)
[ temp ] = mask_with_FS_label( DATAOUT.Rvec, aparcaseg, fslabels.wm.all );
DSTATS.Rmean.wm.all = mean(temp,1);
DSTATS.Rmean_std.wm.all = std(temp,0,1);

[ temp ] = mask_with_FS_label( DATAOUT.Rvec, aparcaseg, fslabels.wm.lh );
DSTATS.Rmean.wm.lh = mean(temp,1);
DSTATS.Rmean_std.wm.lh = std(temp,0,1);

[ temp ] = mask_with_FS_label( DATAOUT.Rvec, aparcaseg, fslabels.wm.rh );
DSTATS.Rmean.wm.rh = mean(temp,1);
DSTATS.Rmean_std.wm.rh = std(temp,0,1);


%%% mean 1d GenAniso profile (+ std across labels)
[ temp ] = mask_with_FS_label( DATAOUT.Rstd, aparcaseg, fslabels.wm.all );
DSTATS.Rstd.wm.all = mean(temp,1);
DSTATS.Rstd_std.wm.all = std(temp,0,1);

[ temp ] = mask_with_FS_label( DATAOUT.Rstd, aparcaseg, fslabels.wm.lh );
DSTATS.Rstd.wm.lh = mean(temp,1);
DSTATS.Rstd_std.wm.lh = std(temp,0,1);

[ temp ] = mask_with_FS_label( DATAOUT.Rstd, aparcaseg, fslabels.wm.rh );
DSTATS.Rstd.wm.rh = mean(temp,1);
DSTATS.Rstd_std.wm.rh = std(temp,0,1);


%%% mean GaussDivergence (+ std across labels)
[ temp ] = mask_with_FS_label( (DATAOUT.GaussDiv), aparcaseg, fslabels.wm.all );
DSTATS.GaussDiv.wm.all = mean(temp,1);
DSTATS.GaussDiv_std.wm.all = std(temp,0,1);

[ temp ] = mask_with_FS_label( (DATAOUT.GaussDiv), aparcaseg, fslabels.wm.lh );
DSTATS.GaussDiv.wm.lh = mean(temp,1);
DSTATS.GaussDiv_std.wm.lh = std(temp,0,1);

[ temp ] = mask_with_FS_label( (DATAOUT.GaussDiv), aparcaseg, fslabels.wm.rh );
DSTATS.GaussDiv.wm.rh = mean(temp,1);
DSTATS.GaussDiv_std.wm.rh = std(temp,0,1);


%%% mean 1d EAPKurtosis profile (+ std across labels)
[ temp ] = mask_with_FS_label( (DATAOUT.KurtProfile), aparcaseg, fslabels.wm.all );
DSTATS.KurtProfile.wm.all = mean(temp,1);
DSTATS.KurtProfile_std.wm.all = std(temp,0,1);

[ temp ] = mask_with_FS_label( (DATAOUT.KurtProfile), aparcaseg, fslabels.wm.lh );
DSTATS.KurtProfile.wm.lh = mean(temp,1);
DSTATS.KurtProfile_std.wm.lh = std(temp,0,1);

[ temp ] = mask_with_FS_label( (DATAOUT.KurtProfile), aparcaseg, fslabels.wm.rh );
DSTATS.KurtProfile.wm.rh = mean(temp,1);
DSTATS.KurtProfile_std.wm.rh = std(temp,0,1);


%%% mean EAPKurtosis  (+ std across labels)
[ temp ] = mask_with_FS_label( (DATAOUT.MeanKurt), aparcaseg, fslabels.wm.all );
DSTATS.MeanKurt.wm.all = mean(temp,1);
DSTATS.MeanKurt_std.wm.all = std(temp,0,1);

[ temp ] = mask_with_FS_label( (DATAOUT.MeanKurt), aparcaseg, fslabels.wm.lh );
DSTATS.MeanKurt.wm.lh = mean(temp,1);
DSTATS.MeanKurt_std.wm.lh = std(temp,0,1);

[ temp ] = mask_with_FS_label( (DATAOUT.MeanKurt), aparcaseg, fslabels.wm.rh );
DSTATS.MeanKurt.wm.rh = mean(temp,1);
DSTATS.MeanKurt_std.wm.rh = std(temp,0,1);



%%%% ctx
%%% mean 1d eap profile (+ std across labels)
[ temp ] = mask_with_FS_label( DATAOUT.Rvec, aparcaseg, fslabels.ctx.all );
DSTATS.Rmean.ctx.all = mean(temp,1);
DSTATS.Rmean_std.ctx.all = std(temp,0,1);

[ temp ] = mask_with_FS_label( DATAOUT.Rvec, aparcaseg, fslabels.ctx.lh );
DSTATS.Rmean.ctx.lh = mean(temp,1);
DSTATS.Rmean_std.ctx.lh = std(temp,0,1);

[ temp ] = mask_with_FS_label( DATAOUT.Rvec, aparcaseg, fslabels.ctx.rh );
DSTATS.Rmean.ctx.rh = mean(temp,1);
DSTATS.Rmean_std.ctx.rh = std(temp,0,1);


%%% mean 1d GenAniso profile (+ std across labels)
[ temp ] = mask_with_FS_label( DATAOUT.Rstd, aparcaseg, fslabels.ctx.all );
DSTATS.Rstd.ctx.all = mean(temp,1);
DSTATS.Rstd_std.ctx.all = std(temp,0,1);

[ temp ] = mask_with_FS_label( DATAOUT.Rstd, aparcaseg, fslabels.ctx.lh );
DSTATS.Rstd.ctx.lh = mean(temp,1);
DSTATS.Rstd_std.ctx.lh = std(temp,0,1);

[ temp ] = mask_with_FS_label( DATAOUT.Rstd, aparcaseg, fslabels.ctx.rh );
DSTATS.Rstd.ctx.rh = mean(temp,1);
DSTATS.Rstd_std.ctx.rh = std(temp,0,1);


%%% mean GaussDivergence (+ std across labels)
[ temp ] = mask_with_FS_label( (DATAOUT.GaussDiv), aparcaseg, fslabels.ctx.all );
DSTATS.GaussDiv.ctx.all = mean(temp,1);
DSTATS.GaussDiv_std.ctx.all = std(temp,0,1);

[ temp ] = mask_with_FS_label( (DATAOUT.GaussDiv), aparcaseg, fslabels.ctx.lh );
DSTATS.GaussDiv.ctx.lh = mean(temp,1);
DSTATS.GaussDiv_std.ctx.lh = std(temp,0,1);

[ temp ] = mask_with_FS_label( (DATAOUT.GaussDiv), aparcaseg, fslabels.ctx.rh );
DSTATS.GaussDiv.ctx.rh = mean(temp,1);
DSTATS.GaussDiv_std.ctx.rh = std(temp,0,1);


%%% mean 1d EAPKurtosis profile (+ std across labels)
[ temp ] = mask_with_FS_label( (DATAOUT.KurtProfile), aparcaseg, fslabels.ctx.all );
DSTATS.KurtProfile.ctx.all = mean(temp,1);
DSTATS.KurtProfile_std.ctx.all = std(temp,0,1);

[ temp ] = mask_with_FS_label( (DATAOUT.KurtProfile), aparcaseg, fslabels.ctx.lh );
DSTATS.KurtProfile.ctx.lh = mean(temp,1);
DSTATS.KurtProfile_std.ctx.lh = std(temp,0,1);

[ temp ] = mask_with_FS_label( (DATAOUT.KurtProfile), aparcaseg, fslabels.ctx.rh );
DSTATS.KurtProfile.ctx.rh = mean(temp,1);
DSTATS.KurtProfile_std.ctx.rh = std(temp,0,1);


%%% mean EAPKurtosis  (+ std across labels)
[ temp ] = mask_with_FS_label( (DATAOUT.MeanKurt), aparcaseg, fslabels.ctx.all );
DSTATS.MeanKurt.ctx.all = mean(temp,1);
DSTATS.MeanKurt_std.ctx.all = std(temp,0,1);

[ temp ] = mask_with_FS_label( (DATAOUT.MeanKurt), aparcaseg, fslabels.ctx.lh );
DSTATS.MeanKurt.ctx.lh = mean(temp,1);
DSTATS.MeanKurt_std.ctx.lh = std(temp,0,1);

[ temp ] = mask_with_FS_label( (DATAOUT.MeanKurt), aparcaseg, fslabels.ctx.rh );
DSTATS.MeanKurt.ctx.rh = mean(temp,1);
DSTATS.MeanKurt_std.ctx.rh = std(temp,0,1);



%%%% CSF
%%% mean 1d eap profile (+ std across labels)
[ temp ] = mask_with_FS_label( DATAOUT.Rvec, aparcaseg, fslabels.csf.all );
DSTATS.Rmean.csf.all = mean(temp,1);
DSTATS.Rmean_std.csf.all = std(temp,0,1);

[ temp ] = mask_with_FS_label( DATAOUT.Rvec, aparcaseg, fslabels.csf.orig );
DSTATS.Rmean.csf.orig = mean(temp,1);
DSTATS.Rmean_std.csf.orig = std(temp,0,1);


%%% mean 1d GenAniso profile (+ std across labels)
[ temp ] = mask_with_FS_label( DATAOUT.Rstd, aparcaseg, fslabels.csf.all );
DSTATS.Rstd.csf.all = mean(temp,1);
DSTATS.Rstd_std.csf.all = std(temp,0,1);

[ temp ] = mask_with_FS_label( DATAOUT.Rstd, aparcaseg, fslabels.csf.orig );
DSTATS.Rstd.csf.orig = mean(temp,1);
DSTATS.Rstd_std.csf.orig = std(temp,0,1);


%%% mean GaussDivergence (+ std across labels)
[ temp ] = mask_with_FS_label( (DATAOUT.GaussDiv), aparcaseg, fslabels.csf.all );
DSTATS.GaussDiv.csf.all = mean(temp,1);
DSTATS.GaussDiv_std.csf.all = std(temp,0,1);

[ temp ] = mask_with_FS_label( (DATAOUT.GaussDiv), aparcaseg, fslabels.csf.orig );
DSTATS.GaussDiv.csf.orig = mean(temp,1);
DSTATS.GaussDiv_std.csf.orig = std(temp,0,1);


%%% mean 1d EAPKurtosis profile (+ std across labels)
[ temp ] = mask_with_FS_label( (DATAOUT.KurtProfile), aparcaseg, fslabels.csf.all );
DSTATS.KurtProfile.csf.all = mean(temp,1);
DSTATS.KurtProfile_std.csf.all = std(temp,0,1);

[ temp ] = mask_with_FS_label( (DATAOUT.KurtProfile), aparcaseg, fslabels.csf.orig );
DSTATS.KurtProfile.csf.orig = mean(temp,1);
DSTATS.KurtProfile_std.csf.orig = std(temp,0,1);


%%% mean EAPKurtosis  (+ std across labels)
[ temp ] = mask_with_FS_label( (DATAOUT.MeanKurt), aparcaseg, fslabels.csf.all );
DSTATS.MeanKurt.csf.all = mean(temp,1);
DSTATS.MeanKurt_std.csf.all = std(temp,0,1);

[ temp ] = mask_with_FS_label( (DATAOUT.MeanKurt), aparcaseg, fslabels.csf.orig );
DSTATS.MeanKurt.csf.orig = mean(temp,1);
DSTATS.MeanKurt_std.csf.orig = std(temp,0,1);


end
