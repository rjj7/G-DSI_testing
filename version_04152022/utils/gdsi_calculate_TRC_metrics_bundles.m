function [ DSTATS ] = gdsi_calculate_TRC_metrics_bundles( DSTATS, DATAOUT, trc )



[ data_masked ] = mask_with_FS_trc( DATAOUT.avg_disp_sph, trc.cst.all );
DSTATS.Rmean.cst = mean(data_masked,1);
DSTATS.Rmean_std.cst = std(data_masked,0,1);

[ data_masked ] = mask_with_FS_trc( DATAOUT.std_disp_sph, trc.cst.all );
DSTATS.Rstd.cst = mean(data_masked,1);

[ data_masked ] = mask_with_FS_trc( DATAOUT.isokldiv, trc.cst.all );
DSTATS.isodiv.cst = mean(data_masked,1);
DSTATS.isodiv_std.cst = std(data_masked,0,1);

[ data_masked ] = mask_with_FS_trc( DATAOUT.eap_kurtosis, trc.cst.all );
DSTATS.eapkurt.cst = mean(data_masked,1);
DSTATS.eapkurt_std.cst = std(data_masked,0,1);




[ data_masked ] = mask_with_FS_trc( DATAOUT.avg_disp_sph, trc.slfp.all );
DSTATS.Rmean.slfp = mean(data_masked,1);
DSTATS.Rmean_std.slfp = std(data_masked,0,1);

[ data_masked ] = mask_with_FS_trc( DATAOUT.std_disp_sph, trc.slfp.all );
DSTATS.Rstd.slfp = mean(data_masked,1);

[ data_masked ] = mask_with_FS_trc( DATAOUT.isokldiv, trc.slfp.all );
DSTATS.isodiv.slfp = mean(data_masked,1);
DSTATS.isodiv_std.slfp = std(data_masked,0,1);

[ data_masked ] = mask_with_FS_trc( DATAOUT.eap_kurtosis, trc.slfp.all );
DSTATS.eapkurt.slfp = mean(data_masked,1);
DSTATS.eapkurt_std.slfp = std(data_masked,0,1);





[ data_masked ] = mask_with_FS_trc( DATAOUT.avg_disp_sph, trc.slft.all );
DSTATS.Rmean.slft = mean(data_masked,1);
DSTATS.Rmean_std.slft = std(data_masked,0,1);

[ data_masked ] = mask_with_FS_trc( DATAOUT.std_disp_sph, trc.slft.all );
DSTATS.Rstd.slft = mean(data_masked,1);

[ data_masked ] = mask_with_FS_trc( DATAOUT.isokldiv, trc.slft.all );
DSTATS.isodiv.slft = mean(data_masked,1);
DSTATS.isodiv_std.slft = std(data_masked,0,1);

[ data_masked ] = mask_with_FS_trc( DATAOUT.eap_kurtosis, trc.slft.all );
DSTATS.eapkurt.slft = mean(data_masked,1);
DSTATS.eapkurt_std.slft = std(data_masked,0,1);





[ data_masked ] = mask_with_FS_trc( DATAOUT.avg_disp_sph, trc.cab.all );
DSTATS.Rmean.cab = mean(data_masked,1);
DSTATS.Rmean_std.cab = std(data_masked,0,1);

[ data_masked ] = mask_with_FS_trc( DATAOUT.std_disp_sph, trc.cab.all );
DSTATS.Rstd.cab = mean(data_masked,1);

[ data_masked ] = mask_with_FS_trc( DATAOUT.isokldiv, trc.cab.all );
DSTATS.isodiv.cab = mean(data_masked,1);
DSTATS.isodiv_std.cab = std(data_masked,0,1);

[ data_masked ] = mask_with_FS_trc( DATAOUT.eap_kurtosis, trc.cab.all );
DSTATS.eapkurt.cab = mean(data_masked,1);
DSTATS.eapkurt_std.cab = std(data_masked,0,1);




[ data_masked ] = mask_with_FS_trc( DATAOUT.avg_disp_sph, trc.ccg.all );
DSTATS.Rmean.ccg = mean(data_masked,1);
DSTATS.Rmean_std.ccg = std(data_masked,0,1);

[ data_masked ] = mask_with_FS_trc( DATAOUT.std_disp_sph, trc.ccg.all );
DSTATS.Rstd.ccg = mean(data_masked,1);

[ data_masked ] = mask_with_FS_trc( DATAOUT.isokldiv, trc.ccg.all );
DSTATS.isodiv.ccg = mean(data_masked,1);
DSTATS.isodiv_std.ccg = std(data_masked,0,1);

[ data_masked ] = mask_with_FS_trc( DATAOUT.eap_kurtosis, trc.ccg.all );
DSTATS.eapkurt.ccg = mean(data_masked,1);
DSTATS.eapkurt_std.ccg = std(data_masked,0,1);





[ data_masked ] = mask_with_FS_trc( DATAOUT.avg_disp_sph, trc.ilf.all );
DSTATS.Rmean.ilf = mean(data_masked,1);
DSTATS.Rmean_std.ilf = std(data_masked,0,1);

[ data_masked ] = mask_with_FS_trc( DATAOUT.std_disp_sph, trc.ilf.all );
DSTATS.Rstd.ilf = mean(data_masked,1);

[ data_masked ] = mask_with_FS_trc( DATAOUT.isokldiv, trc.ilf.all );
DSTATS.isodiv.ilf = mean(data_masked,1);
DSTATS.isodiv_std.ilf = std(data_masked,0,1);

[ data_masked ] = mask_with_FS_trc( DATAOUT.eap_kurtosis, trc.ilf.all );
DSTATS.eapkurt.ilf = mean(data_masked,1);
DSTATS.eapkurt_std.ilf = std(data_masked,0,1);





[ data_masked ] = mask_with_FS_trc( DATAOUT.avg_disp_sph, trc.unc.all );
DSTATS.Rmean.unc = mean(data_masked,1);
DSTATS.Rmean_std.unc = std(data_masked,0,1);

[ data_masked ] = mask_with_FS_trc( DATAOUT.std_disp_sph, trc.unc.all );
DSTATS.Rstd.unc = mean(data_masked,1);

[ data_masked ] = mask_with_FS_trc( DATAOUT.isokldiv, trc.unc.all );
DSTATS.isodiv.unc = mean(data_masked,1);
DSTATS.isodiv_std.unc = std(data_masked,0,1);

[ data_masked ] = mask_with_FS_trc( DATAOUT.eap_kurtosis, trc.unc.all );
DSTATS.eapkurt.unc = mean(data_masked,1);
DSTATS.eapkurt_std.unc = std(data_masked,0,1);





[ data_masked ] = mask_with_FS_trc( DATAOUT.avg_disp_sph, trc.fminor );
DSTATS.Rmean.fminor = mean(data_masked,1);
DSTATS.Rmean_std.fminor = std(data_masked,0,1);

[ data_masked ] = mask_with_FS_trc( DATAOUT.std_disp_sph, trc.fminor );
DSTATS.Rstd.fminor = mean(data_masked,1);

[ data_masked ] = mask_with_FS_trc( DATAOUT.isokldiv, trc.fminor );
DSTATS.isodiv.fminor = mean(data_masked,1);
DSTATS.isodiv_std.fminor = std(data_masked,0,1);

[ data_masked ] = mask_with_FS_trc( DATAOUT.eap_kurtosis, trc.fminor );
DSTATS.eapkurt.fminor = mean(data_masked,1);
DSTATS.eapkurt_std.fminor = std(data_masked,0,1);






[ data_masked ] = mask_with_FS_trc( DATAOUT.avg_disp_sph, trc.fmajor );
DSTATS.Rmean.fmajor = mean(data_masked,1);
DSTATS.Rmean_std.fmajor = std(data_masked,0,1);

[ data_masked ] = mask_with_FS_trc( DATAOUT.std_disp_sph, trc.fmajor );
DSTATS.Rstd.fmajor = mean(data_masked,1);

[ data_masked ] = mask_with_FS_trc( DATAOUT.isokldiv, trc.fmajor );
DSTATS.isodiv.fmajor = mean(data_masked,1);
DSTATS.isodiv_std.fmajor = std(data_masked,0,1);

[ data_masked ] = mask_with_FS_trc( DATAOUT.eap_kurtosis, trc.fminor );
DSTATS.eapkurt.fmajor = mean(data_masked,1);
DSTATS.eapkurt_std.fmajor = std(data_masked,0,1);


end


