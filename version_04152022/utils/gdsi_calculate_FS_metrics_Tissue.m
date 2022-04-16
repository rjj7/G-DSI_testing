function [ DSTATS ] = gdsi_calculate_FS_metrics_Tissue( DSTATS, DATAOUT, aparcaseg, fslabels )


%%%% WM
[ mean1d.wm.all ] = mask_with_FS_label( DATAOUT.avg_disp_sph, aparcaseg, fslabels.wm.all );
DSTATS.Rmean.wm = mean(mean1d.wm.all,1);
DSTATS.Rmean_std.wm = std(mean1d.wm.all,0,1);

[ std1d.wm.all ] = mask_with_FS_label( DATAOUT.std_disp_sph, aparcaseg, fslabels.wm.all );
DSTATS.Rstd.wm = mean(std1d.wm.all,1);

[ kld1d.wm.all ] = mask_with_FS_label( real(DATAOUT.isokldiv), aparcaseg, fslabels.wm.all );
tmpmin = min(kld1d.wm.all,[],2);
negind = tmpmin<0;
kld1d.wm.all(negind,:)=[];
DSTATS.isodiv.wm = mean(kld1d.wm.all,1);
DSTATS.isodiv_std.wm = std(kld1d.wm.all,0,1);

[ kurt1d.wm.all ] = mask_with_FS_label( (DATAOUT.eap_kurtosis), aparcaseg, fslabels.wm.all );
DSTATS.eapkurt.wm = mean(kurt1d.wm.all,1);
DSTATS.eapkurt_std.wm = std(kurt1d.wm.all,0,1);

%%%% Ctx
[ mean1d.ctx.all ] = mask_with_FS_label( DATAOUT.avg_disp_sph, aparcaseg, fslabels.ctx.all );
DSTATS.Rmean.ctx = mean(mean1d.ctx.all,1);
DSTATS.Rmean_std.ctx = std(mean1d.ctx.all,0,1);

[ std1d.ctx.all ] = mask_with_FS_label( DATAOUT.std_disp_sph, aparcaseg, fslabels.ctx.all );
DSTATS.Rstd.ctx = mean(std1d.ctx.all,1);

[ kld1d.ctx.all ] = mask_with_FS_label( real(DATAOUT.isokldiv), aparcaseg, fslabels.ctx.all );
tmpmin = min(kld1d.ctx.all,[],2);
negind = tmpmin<0;
kld1d.ctx.all(negind,:)=[];
DSTATS.isodiv.ctx = mean(kld1d.ctx.all,1);
DSTATS.isodiv_std.ctx = std(kld1d.ctx.all,0,1);

[ kurt1d.ctx.all ] = mask_with_FS_label( (DATAOUT.eap_kurtosis), aparcaseg, fslabels.ctx.all );
DSTATS.eapkurt.ctx = mean(kurt1d.ctx.all,1);
DSTATS.eapkurt_std.ctx = std(kurt1d.ctx.all,0,1);

%%%% CSF
[ mean1d.csf.all ] = mask_with_FS_label( DATAOUT.avg_disp_sph, aparcaseg, fslabels.csf.all );
DSTATS.Rmean.csf = mean(mean1d.csf.all,1);
DSTATS.Rmean_std.csf = std(mean1d.csf.all,0,1);

[ std1d.csf.all ] = mask_with_FS_label( DATAOUT.std_disp_sph, aparcaseg, fslabels.csf.all );
DSTATS.Rstd.csf = mean(std1d.csf.all,1);

[ kld1d.csf.all ] = mask_with_FS_label( real(DATAOUT.isokldiv), aparcaseg, fslabels.csf.all );
tmpmin = min(kld1d.csf.all,[],2);
negind = tmpmin<0;
kld1d.csf.all(negind,:)=[];
DSTATS.isodiv.csf = mean(kld1d.csf.all,1);
DSTATS.isodiv_std.csf = std(kld1d.csf.all,0,1);

[ kurt1d.csf.all ] = mask_with_FS_label( (DATAOUT.eap_kurtosis), aparcaseg, fslabels.csf.all );
DSTATS.eapkurt.csf = mean(kurt1d.csf.all,1);
DSTATS.eapkurt_std.csf = std(kurt1d.csf.all,0,1);


end
