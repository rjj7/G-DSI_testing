function [ DSTATS ] = gdsi_calculate_TRC_metrics_bundles_v2( DSTATS, DATAOUT, trc )



[ data_masked ] = mask_with_FS_trc( DATAOUT.avg_disp_sph, trc.cst.all );
DSTATS.Rmean.cst = mean(data_masked,1);
DSTATS.Rmean_std.cst = std(data_masked,0,1);

[ data_masked ] = mask_with_FS_trc( DATAOUT.std_disp_sph, trc.cst.all );
DSTATS.Rstd.cst = mean(data_masked,1);

[ data_masked ] = mask_with_FS_trc( (DATAOUT.GaussDiv), trc.cst.all );
DSTATS.GaussDiv.cst = mean(data_masked,1);
DSTATS.GaussDiv_std.cst = std(data_masked,0,1);

[ data_masked ] = mask_with_FS_trc( (DATAOUT.KurtProfile), trc.cst.all );
DSTATS.KurtProfile.cst = mean(data_masked,1);
DSTATS.KurtProfile_std.cst = std(data_masked,0,1);

[ data_masked ] = mask_with_FS_trc( (DATAOUT.MeanKurt), trc.cst.all );
DSTATS.MeanKurt.cst = mean(data_masked,1);
DSTATS.MeanKurt_std.cst = std(data_masked,0,1);




[ data_masked ] = mask_with_FS_trc( DATAOUT.avg_disp_sph, trc.slfp.all );
DSTATS.Rmean.slfp = mean(data_masked,1);
DSTATS.Rmean_std.slfp = std(data_masked,0,1);

[ data_masked ] = mask_with_FS_trc( DATAOUT.std_disp_sph, trc.slfp.all );
DSTATS.Rstd.slfp = mean(data_masked,1);

[ data_masked ] = mask_with_FS_trc( (DATAOUT.GaussDiv), trc.slfp.all );
DSTATS.GaussDiv.slfp = mean(data_masked,1);
DSTATS.GaussDiv_std.slfp = std(data_masked,0,1);

[ data_masked ] = mask_with_FS_trc( (DATAOUT.KurtProfile), trc.slfp.all );
DSTATS.KurtProfile.slfp = mean(data_masked,1);
DSTATS.KurtProfile_std.slfp = std(data_masked,0,1);

[ data_masked ] = mask_with_FS_trc( (DATAOUT.MeanKurt), trc.slfp.all );
DSTATS.MeanKurt.slfp = mean(data_masked,1);
DSTATS.MeanKurt_std.slfp = std(data_masked,0,1);





[ data_masked ] = mask_with_FS_trc( DATAOUT.avg_disp_sph, trc.slft.all );
DSTATS.Rmean.slft = mean(data_masked,1);
DSTATS.Rmean_std.slft = std(data_masked,0,1);

[ data_masked ] = mask_with_FS_trc( DATAOUT.std_disp_sph, trc.slft.all );
DSTATS.Rstd.slft = mean(data_masked,1);

[ data_masked ] = mask_with_FS_trc( (DATAOUT.GaussDiv), trc.slft.all );
DSTATS.GaussDiv.slft = mean(data_masked,1);
DSTATS.GaussDiv_std.slft = std(data_masked,0,1);

[ data_masked ] = mask_with_FS_trc( (DATAOUT.KurtProfile), trc.slft.all );
DSTATS.KurtProfile.slft = mean(data_masked,1);
DSTATS.KurtProfile_std.slft = std(data_masked,0,1);

[ data_masked ] = mask_with_FS_trc( (DATAOUT.MeanKurt), trc.slft.all );
DSTATS.MeanKurt.slft = mean(data_masked,1);
DSTATS.MeanKurt_std.slft = std(data_masked,0,1);






[ data_masked ] = mask_with_FS_trc( DATAOUT.avg_disp_sph, trc.cab.all );
DSTATS.Rmean.cab = mean(data_masked,1);
DSTATS.Rmean_std.cab = std(data_masked,0,1);

[ data_masked ] = mask_with_FS_trc( DATAOUT.std_disp_sph, trc.cab.all );
DSTATS.Rstd.cab = mean(data_masked,1);

[ data_masked ] = mask_with_FS_trc( (DATAOUT.GaussDiv), trc.cab.all );
DSTATS.GaussDiv.cab = mean(data_masked,1);
DSTATS.GaussDiv_std.cab = std(data_masked,0,1);

[ data_masked ] = mask_with_FS_trc( (DATAOUT.KurtProfile), trc.cab.all );
DSTATS.KurtProfile.cab = mean(data_masked,1);
DSTATS.KurtProfile_std.cab = std(data_masked,0,1);

[ data_masked ] = mask_with_FS_trc( (DATAOUT.MeanKurt), trc.cab.all );
DSTATS.MeanKurt.cab = mean(data_masked,1);
DSTATS.MeanKurt_std.cab = std(data_masked,0,1);




[ data_masked ] = mask_with_FS_trc( DATAOUT.avg_disp_sph, trc.ccg.all );
DSTATS.Rmean.ccg = mean(data_masked,1);
DSTATS.Rmean_std.ccg = std(data_masked,0,1);

[ data_masked ] = mask_with_FS_trc( DATAOUT.std_disp_sph, trc.ccg.all );
DSTATS.Rstd.ccg = mean(data_masked,1);

[ data_masked ] = mask_with_FS_trc( (DATAOUT.GaussDiv), trc.ccg.all );
DSTATS.GaussDiv.ccg = mean(data_masked,1);
DSTATS.GaussDiv_std.ccg = std(data_masked,0,1);

[ data_masked ] = mask_with_FS_trc( (DATAOUT.KurtProfile), trc.ccg.all );
DSTATS.KurtProfile.ccg = mean(data_masked,1);
DSTATS.KurtProfile_std.ccg = std(data_masked,0,1);

[ data_masked ] = mask_with_FS_trc( (DATAOUT.MeanKurt), trc.ccg.all );
DSTATS.MeanKurt.ccg = mean(data_masked,1);
DSTATS.MeanKurt_std.ccg = std(data_masked,0,1);





[ data_masked ] = mask_with_FS_trc( DATAOUT.avg_disp_sph, trc.ilf.all );
DSTATS.Rmean.ilf = mean(data_masked,1);
DSTATS.Rmean_std.ilf = std(data_masked,0,1);

[ data_masked ] = mask_with_FS_trc( DATAOUT.std_disp_sph, trc.ilf.all );
DSTATS.Rstd.ilf = mean(data_masked,1);

[ data_masked ] = mask_with_FS_trc( (DATAOUT.GaussDiv), trc.ilf.all );
DSTATS.GaussDiv.ilf = mean(data_masked,1);
DSTATS.GaussDiv_std.ilf = std(data_masked,0,1);

[ data_masked ] = mask_with_FS_trc( (DATAOUT.KurtProfile), trc.ilf.all );
DSTATS.KurtProfile.ilf = mean(data_masked,1);
DSTATS.KurtProfile_std.ilf = std(data_masked,0,1);

[ data_masked ] = mask_with_FS_trc( (DATAOUT.MeanKurt), trc.ilf.all );
DSTATS.MeanKurt.ilf = mean(data_masked,1);
DSTATS.MeanKurt_std.ilf = std(data_masked,0,1);





[ data_masked ] = mask_with_FS_trc( DATAOUT.avg_disp_sph, trc.unc.all );
DSTATS.Rmean.unc = mean(data_masked,1);
DSTATS.Rmean_std.unc = std(data_masked,0,1);

[ data_masked ] = mask_with_FS_trc( DATAOUT.std_disp_sph, trc.unc.all );
DSTATS.Rstd.unc = mean(data_masked,1);

[ data_masked ] = mask_with_FS_trc( (DATAOUT.GaussDiv), trc.unc.all );
DSTATS.GaussDiv.unc = mean(data_masked,1);
DSTATS.GaussDiv_std.unc = std(data_masked,0,1);

[ data_masked ] = mask_with_FS_trc( (DATAOUT.KurtProfile), trc.unc.all );
DSTATS.KurtProfile.unc = mean(data_masked,1);
DSTATS.KurtProfile_std.unc = std(data_masked,0,1);

[ data_masked ] = mask_with_FS_trc( (DATAOUT.MeanKurt), trc.unc.all );
DSTATS.MeanKurt.unc = mean(data_masked,1);
DSTATS.MeanKurt_std.unc = std(data_masked,0,1);





[ data_masked ] = mask_with_FS_trc( DATAOUT.avg_disp_sph, trc.fminor );
DSTATS.Rmean.fminor = mean(data_masked,1);
DSTATS.Rmean_std.fminor = std(data_masked,0,1);

[ data_masked ] = mask_with_FS_trc( DATAOUT.std_disp_sph, trc.fminor );
DSTATS.Rstd.fminor = mean(data_masked,1);

[ data_masked ] = mask_with_FS_trc( (DATAOUT.GaussDiv), trc.fminor );
DSTATS.GaussDiv.fminor = mean(data_masked,1);
DSTATS.GaussDiv_std.fminor = std(data_masked,0,1);

[ data_masked ] = mask_with_FS_trc( (DATAOUT.KurtProfile), trc.fminor );
DSTATS.KurtProfile.fminor = mean(data_masked,1);
DSTATS.KurtProfile_std.fminor = std(data_masked,0,1);

[ data_masked ] = mask_with_FS_trc( (DATAOUT.MeanKurt), trc.fminor );
DSTATS.MeanKurt.fminor = mean(data_masked,1);
DSTATS.MeanKurt_std.fminor = std(data_masked,0,1);






[ data_masked ] = mask_with_FS_trc( DATAOUT.avg_disp_sph, trc.fmajor );
DSTATS.Rmean.fmajor = mean(data_masked,1);
DSTATS.Rmean_std.fmajor = std(data_masked,0,1);

[ data_masked ] = mask_with_FS_trc( DATAOUT.std_disp_sph, trc.fmajor );
DSTATS.Rstd.fmajor = mean(data_masked,1);

[ data_masked ] = mask_with_FS_trc( (DATAOUT.GaussDiv), trc.fmajor );
DSTATS.GaussDiv.fmajor = mean(data_masked,1);
DSTATS.GaussDiv_std.fmajor = std(data_masked,0,1);

[ data_masked ] = mask_with_FS_trc( (DATAOUT.KurtProfile), trc.fmajor );
DSTATS.KurtProfile.fmajor = mean(data_masked,1);
DSTATS.KurtProfile_std.fmajor = std(data_masked,0,1);

[ data_masked ] = mask_with_FS_trc( (DATAOUT.MeanKurt), trc.fmajor );
DSTATS.MeanKurt.fmajor = mean(data_masked,1);
DSTATS.MeanKurt_std.fmajor = std(data_masked,0,1);


end


