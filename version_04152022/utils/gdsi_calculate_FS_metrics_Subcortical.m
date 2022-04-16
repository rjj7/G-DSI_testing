function [ DSTATS ] = gdsi_calculate_FS_metrics_Subcortical( DSTATS, DATAOUT, aparcaseg, fslabels )


%%%% brainstem
[ mean1d.brainstem.all ] = mask_with_FS_label( DATAOUT.avg_disp_sph, aparcaseg, fslabels.brainstem );
DSTATS.Rmean.brainstem = mean(mean1d.brainstem.all,1);
DSTATS.Rmean_std.brainstem = std(mean1d.brainstem.all,0,1);

[ std1d.brainstem.all ] = mask_with_FS_label( DATAOUT.std_disp_sph, aparcaseg, fslabels.brainstem );
DSTATS.Rstd.brainstem = mean(std1d.brainstem.all,1);

[ kld1d.brainstem.all ] = mask_with_FS_label( real(DATAOUT.isokldiv), aparcaseg, fslabels.brainstem );
DSTATS.isodiv.brainstem = mean(kld1d.brainstem.all,1);
DSTATS.isodiv_std.brainstem = std(kld1d.brainstem.all,0,1);

[ kurt1d.brainstem.all ] = mask_with_FS_label( (DATAOUT.eap_kurtosis), aparcaseg, fslabels.brainstem );
DSTATS.eapkurt.brainstem = mean(kurt1d.brainstem.all,1);
DSTATS.eapkurt_std.brainstem = std(kurt1d.brainstem.all,0,1);


%%%% cc
[ mean1d.cc.all ] = mask_with_FS_label( DATAOUT.avg_disp_sph, aparcaseg, fslabels.cc );
DSTATS.Rmean.cc = mean(mean1d.cc.all,1);
DSTATS.Rmean_std.cc = std(mean1d.cc.all,0,1);

[ std1d.cc.all ] = mask_with_FS_label( DATAOUT.std_disp_sph, aparcaseg, fslabels.cc );
DSTATS.Rstd.cc = mean(std1d.cc.all,1);

[ kld1d.cc.all ] = mask_with_FS_label( real(DATAOUT.isokldiv), aparcaseg, fslabels.cc );
DSTATS.isodiv.cc = mean(kld1d.cc.all,1);
DSTATS.isodiv_std.cc = std(kld1d.cc.all,0,1);

[ kurt1d.cc.all ] = mask_with_FS_label( (DATAOUT.eap_kurtosis), aparcaseg, fslabels.cc );
DSTATS.eapkurt.cc = mean(kurt1d.cc.all,1);
DSTATS.eapkurt_std.cc = std(kurt1d.cc.all,0,1);


%%%% cerebellum
[ mean1d.cerebellum.all ] = mask_with_FS_label( DATAOUT.avg_disp_sph, aparcaseg, fslabels.cerebellum.all );
DSTATS.Rmean.cerebellum = mean(mean1d.cerebellum.all,1);
DSTATS.Rmean_std.cerebellum = std(mean1d.cerebellum.all,0,1);

[ std1d.cerebellum.all ] = mask_with_FS_label( DATAOUT.std_disp_sph, aparcaseg, fslabels.cerebellum.all );
DSTATS.Rstd.cerebellum = mean(std1d.cerebellum.all,1);

[ kld1d.cerebellum.all ] = mask_with_FS_label( real(DATAOUT.isokldiv), aparcaseg, fslabels.cerebellum.all );
DSTATS.isodiv.cerebellum = mean(kld1d.cerebellum.all,1);
DSTATS.isodiv_std.cerebellum = std(kld1d.cerebellum.all,0,1);

[ kurt1d.cerebellum.all ] = mask_with_FS_label( (DATAOUT.eap_kurtosis), aparcaseg, fslabels.cerebellum.all );
DSTATS.eapkurt.cerebellum = mean(kurt1d.cerebellum.all,1);
DSTATS.eapkurt_std.cerebellum = std(kurt1d.cerebellum.all,0,1);


%%%% thalamus
[ mean1d.thalamus.all ] = mask_with_FS_label( DATAOUT.avg_disp_sph, aparcaseg, fslabels.subcort.thalamus.all );
DSTATS.Rmean.thalamus = mean(mean1d.thalamus.all,1);
DSTATS.Rmean_std.thalamus = std(mean1d.thalamus.all,0,1);

[ std1d.thalamus.all ] = mask_with_FS_label( DATAOUT.std_disp_sph, aparcaseg, fslabels.subcort.thalamus.all );
DSTATS.Rstd.thalamus = mean(std1d.thalamus.all,1);

[ kld1d.thalamus.all ] = mask_with_FS_label( real(DATAOUT.isokldiv), aparcaseg, fslabels.subcort.thalamus.all );
DSTATS.isodiv.thalamus = mean(kld1d.thalamus.all,1);
DSTATS.isodiv_std.thalamus = std(kld1d.thalamus.all,0,1);

[ kurt1d.thalamus.all ] = mask_with_FS_label( (DATAOUT.eap_kurtosis), aparcaseg, fslabels.subcort.thalamus.all );
DSTATS.eapkurt.thalamus = mean(kurt1d.thalamus.all,1);
DSTATS.eapkurt_std.thalamus = std(kurt1d.thalamus.all,0,1);


%%%% putamen
[ mean1d.putamen.all ] = mask_with_FS_label( DATAOUT.avg_disp_sph, aparcaseg, fslabels.subcort.putamen.all );
DSTATS.Rmean.putamen = mean(mean1d.putamen.all,1);
DSTATS.Rmean_std.putamen = std(mean1d.putamen.all,0,1);

[ std1d.putamen.all ] = mask_with_FS_label( DATAOUT.std_disp_sph, aparcaseg, fslabels.subcort.putamen.all );
DSTATS.Rstd.putamen = mean(std1d.putamen.all,1);

[ kld1d.putamen.all ] = mask_with_FS_label( real(DATAOUT.isokldiv), aparcaseg, fslabels.subcort.putamen.all );
DSTATS.isodiv.putamen = mean(kld1d.putamen.all,1);
DSTATS.isodiv_std.putamen = std(kld1d.putamen.all,0,1);

[ kurt1d.putamen.all ] = mask_with_FS_label( (DATAOUT.eap_kurtosis), aparcaseg, fslabels.subcort.putamen.all );
DSTATS.eapkurt.putamen = mean(kurt1d.putamen.all,1);
DSTATS.eapkurt_std.putamen = std(kurt1d.putamen.all,0,1);


%%%% caudate
[ mean1d.caudate.all ] = mask_with_FS_label( DATAOUT.avg_disp_sph, aparcaseg, fslabels.subcort.caudate.all );
DSTATS.Rmean.caudate = mean(mean1d.caudate.all,1);
DSTATS.Rmean_std.caudate = std(mean1d.caudate.all,0,1);

[ std1d.caudate.all ] = mask_with_FS_label( DATAOUT.std_disp_sph, aparcaseg, fslabels.subcort.caudate.all );
DSTATS.Rstd.caudate = mean(std1d.caudate.all,1);

[ kld1d.caudate.all ] = mask_with_FS_label( real(DATAOUT.isokldiv), aparcaseg, fslabels.subcort.caudate.all );
DSTATS.isodiv.caudate = mean(kld1d.caudate.all,1);
DSTATS.isodiv_std.caudate = std(kld1d.caudate.all,0,1);

[ kurt1d.caudate.all ] = mask_with_FS_label( (DATAOUT.eap_kurtosis), aparcaseg, fslabels.subcort.caudate.all );
DSTATS.eapkurt.caudate = mean(kurt1d.caudate.all,1);
DSTATS.eapkurt_std.caudate = std(kurt1d.caudate.all,0,1);


%%%% pallidum
[ mean1d.pallidum.all ] = mask_with_FS_label( DATAOUT.avg_disp_sph, aparcaseg, fslabels.subcort.pallidum.all );
DSTATS.Rmean.pallidum = mean(mean1d.pallidum.all,1);
DSTATS.Rmean_std.pallidum = std(mean1d.pallidum.all,0,1);

[ std1d.pallidum.all ] = mask_with_FS_label( DATAOUT.std_disp_sph, aparcaseg, fslabels.subcort.pallidum.all );
DSTATS.Rstd.pallidum = mean(std1d.pallidum.all,1);

[ kld1d.pallidum.all ] = mask_with_FS_label( real(DATAOUT.isokldiv), aparcaseg, fslabels.subcort.pallidum.all );
DSTATS.isodiv.pallidum = mean(kld1d.pallidum.all,1);
DSTATS.isodiv_std.pallidum = std(kld1d.pallidum.all,0,1);

[ kurt1d.pallidum.all ] = mask_with_FS_label( (DATAOUT.eap_kurtosis), aparcaseg, fslabels.subcort.pallidum.all );
DSTATS.eapkurt.pallidum = mean(kurt1d.pallidum.all,1);
DSTATS.eapkurt_std.pallidum = std(kurt1d.pallidum.all,0,1);


%%%% amygdala
[ mean1d.amygdala.all ] = mask_with_FS_label( DATAOUT.avg_disp_sph, aparcaseg, fslabels.subcort.amygdala.all );
DSTATS.Rmean.amygdala = mean(mean1d.amygdala.all,1);
DSTATS.Rmean_std.amygdala = std(mean1d.amygdala.all,0,1);

[ std1d.amygdala.all ] = mask_with_FS_label( DATAOUT.std_disp_sph, aparcaseg, fslabels.subcort.amygdala.all );
DSTATS.Rstd.amygdala = mean(std1d.amygdala.all,1);

[ kld1d.amygdala.all ] = mask_with_FS_label( real(DATAOUT.isokldiv), aparcaseg, fslabels.subcort.amygdala.all );
DSTATS.isodiv.amygdala = mean(kld1d.amygdala.all,1);
DSTATS.isodiv_std.amygdala = std(kld1d.amygdala.all,0,1);

[ kurt1d.amygdala.all ] = mask_with_FS_label( (DATAOUT.eap_kurtosis), aparcaseg, fslabels.subcort.amygdala.all );
DSTATS.eapkurt.amygdala = mean(kurt1d.amygdala.all,1);
DSTATS.eapkurt_std.amygdala = std(kurt1d.amygdala.all,0,1);


%%%% hippocampus
[ mean1d.hippocampus.all ] = mask_with_FS_label( DATAOUT.avg_disp_sph, aparcaseg, fslabels.subcort.hippocampus.all );
DSTATS.Rmean.hippocampus = mean(mean1d.hippocampus.all,1);
DSTATS.Rmean_std.hippocampus = std(mean1d.hippocampus.all,0,1);

[ std1d.hippocampus.all ] = mask_with_FS_label( DATAOUT.std_disp_sph, aparcaseg, fslabels.subcort.hippocampus.all );
DSTATS.Rstd.hippocampus = mean(std1d.hippocampus.all,1);

[ kld1d.hippocampus.all ] = mask_with_FS_label( real(DATAOUT.isokldiv), aparcaseg, fslabels.subcort.hippocampus.all );
DSTATS.isodiv.hippocampus = mean(kld1d.hippocampus.all,1);
DSTATS.isodiv_std.hippocampus = std(kld1d.hippocampus.all,0,1);

[ kurt1d.hippocampus.all ] = mask_with_FS_label( (DATAOUT.eap_kurtosis), aparcaseg, fslabels.subcort.hippocampus.all );
DSTATS.eapkurt.hippocampus = mean(kurt1d.hippocampus.all,1);
DSTATS.eapkurt_std.hippocampus = std(kurt1d.hippocampus.all,0,1);


%%%% accumbens
[ mean1d.accumbens.all ] = mask_with_FS_label( DATAOUT.avg_disp_sph, aparcaseg, fslabels.subcort.accumbens.all );
DSTATS.Rmean.accumbens = mean(mean1d.accumbens.all,1);
DSTATS.Rmean_std.accumbens = std(mean1d.accumbens.all,0,1);

[ std1d.accumbens.all ] = mask_with_FS_label( DATAOUT.std_disp_sph, aparcaseg, fslabels.subcort.accumbens.all );
DSTATS.Rstd.accumbens = mean(std1d.accumbens.all,1);

[ kld1d.accumbens.all ] = mask_with_FS_label( real(DATAOUT.isokldiv), aparcaseg, fslabels.subcort.accumbens.all );
DSTATS.isodiv.accumbens = mean(kld1d.accumbens.all,1);
DSTATS.isodiv_std.accumbens = std(kld1d.accumbens.all,0,1);

[ kurt1d.accumbens.all ] = mask_with_FS_label( (DATAOUT.eap_kurtosis), aparcaseg, fslabels.subcort.accumbens.all );
DSTATS.eapkurt.accumbens = mean(kurt1d.accumbens.all,1);
DSTATS.eapkurt_std.accumbens = std(kurt1d.accumbens.all,0,1);


%%%% ventral dc
[ mean1d.ventraldc.all ] = mask_with_FS_label( DATAOUT.avg_disp_sph, aparcaseg, fslabels.subcort.ventraldc.all );
DSTATS.Rmean.ventraldc = mean(mean1d.ventraldc.all,1);
DSTATS.Rmean_std.ventraldc = std(mean1d.ventraldc.all,0,1);

[ std1d.ventraldc.all ] = mask_with_FS_label( DATAOUT.std_disp_sph, aparcaseg, fslabels.subcort.ventraldc.all );
DSTATS.Rstd.ventraldc = mean(std1d.ventraldc.all,1);

[ kld1d.ventraldc.all ] = mask_with_FS_label( real(DATAOUT.isokldiv), aparcaseg, fslabels.subcort.ventraldc.all );
DSTATS.isodiv.ventraldc = mean(kld1d.ventraldc.all,1);
DSTATS.isodiv_std.ventraldc = std(kld1d.ventraldc.all,0,1);

[ kurt1d.ventraldc.all ] = mask_with_FS_label( (DATAOUT.eap_kurtosis), aparcaseg, fslabels.subcort.ventraldc.all );
DSTATS.eapkurt.ventraldc = mean(kurt1d.ventraldc.all,1);
DSTATS.eapkurt_std.ventraldc = std(kurt1d.ventraldc.all,0,1);



end