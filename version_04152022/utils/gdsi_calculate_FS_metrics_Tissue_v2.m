function [ DSTATS ] = gdsi_calculate_FS_metrics_Tissue_v2( DSTATS, DATAOUT, aparcaseg, fslabels )

% Calculates metrics across different FS tissue labels, save to DSTATS
%  struct
% - mean 1D EAP profile, mean 1D GenAniso profile, etc.
%

%%%% WM
[ mean1d.wm.all ] = mask_with_FS_label( DATAOUT.avg_disp_sph, aparcaseg, fslabels.wm.all );
DSTATS.Rmean.wm = mean(mean1d.wm.all,1);
DSTATS.Rmean_std.wm = std(mean1d.wm.all,0,1);

[ std1d.wm.all ] = mask_with_FS_label( DATAOUT.std_disp_sph, aparcaseg, fslabels.wm.all );
DSTATS.Rstd.wm = mean(std1d.wm.all,1);

[ GaussDiv.wm.all ] = mask_with_FS_label( (DATAOUT.GaussDiv), aparcaseg, fslabels.wm.all );
DSTATS.GaussDiv.wm = mean(GaussDiv.wm.all,1);
DSTATS.GaussDiv_std.wm = std(GaussDiv.wm.all,0,1);

[ KurtProfile.wm.all ] = mask_with_FS_label( (DATAOUT.KurtProfile), aparcaseg, fslabels.wm.all );
DSTATS.KurtProfile.wm = mean(KurtProfile.wm.all,1);
DSTATS.KurtProfile_std.wm = std(KurtProfile.wm.all,0,1);

[ MeanKurt.wm.all ] = mask_with_FS_label( (DATAOUT.MeanKurt), aparcaseg, fslabels.wm.all );
DSTATS.MeanKurt.wm = mean(MeanKurt.wm.all,1);
DSTATS.MeanKurt_std.wm = std(MeanKurt.wm.all,0,1);



%%%% Ctx
[ mean1d.ctx.all ] = mask_with_FS_label( DATAOUT.avg_disp_sph, aparcaseg, fslabels.ctx.all );
DSTATS.Rmean.ctx = mean(mean1d.ctx.all,1);
DSTATS.Rmean_std.ctx = std(mean1d.ctx.all,0,1);

[ std1d.ctx.all ] = mask_with_FS_label( DATAOUT.std_disp_sph, aparcaseg, fslabels.ctx.all );
DSTATS.Rstd.ctx = mean(std1d.ctx.all,1);

[ GaussDiv.ctx.all ] = mask_with_FS_label( (DATAOUT.GaussDiv), aparcaseg, fslabels.ctx.all );
DSTATS.GaussDiv.ctx = mean(GaussDiv.ctx.all,1);
DSTATS.GaussDiv_std.ctx = std(GaussDiv.ctx.all,0,1);

[ KurtProfile.ctx.all ] = mask_with_FS_label( (DATAOUT.KurtProfile), aparcaseg, fslabels.ctx.all );
DSTATS.KurtProfile.ctx = mean(KurtProfile.ctx.all,1);
DSTATS.KurtProfile_std.ctx = std(KurtProfile.ctx.all,0,1);

[ MeanKurt.ctx.all ] = mask_with_FS_label( (DATAOUT.MeanKurt), aparcaseg, fslabels.ctx.all );
DSTATS.MeanKurt.ctx = mean(MeanKurt.ctx.all,1);
DSTATS.MeanKurt_std.ctx = std(MeanKurt.ctx.all,0,1);



%%%% CSF
[ mean1d.csf.all ] = mask_with_FS_label( DATAOUT.avg_disp_sph, aparcaseg, fslabels.csf.all );
DSTATS.Rmean.csf = mean(mean1d.csf.all,1);
DSTATS.Rmean_std.csf = std(mean1d.csf.all,0,1);

[ std1d.csf.all ] = mask_with_FS_label( DATAOUT.std_disp_sph, aparcaseg, fslabels.csf.all );
DSTATS.Rstd.csf = mean(std1d.csf.all,1);

[ GaussDiv.csf.all ] = mask_with_FS_label( real(DATAOUT.GaussDiv), aparcaseg, fslabels.csf.all );
DSTATS.GaussDiv.csf = mean(GaussDiv.csf.all,1);
DSTATS.GaussDiv_std.csf = std(GaussDiv.csf.all,0,1);

[ KurtProfile.csf.all ] = mask_with_FS_label( (DATAOUT.KurtProfile), aparcaseg, fslabels.csf.all );
DSTATS.KurtProfile.csf = mean(KurtProfile.csf.all,1);
DSTATS.KurtProfile_std.csf = std(KurtProfile.csf.all,0,1);

[ MeanKurt.csf.all ] = mask_with_FS_label( (DATAOUT.MeanKurt), aparcaseg, fslabels.csf.all );
DSTATS.MeanKurt.csf = mean(MeanKurt.csf.all,1);
DSTATS.MeanKurt_std.csf = std(MeanKurt.csf.all,0,1);

end
