function [ DSTATS ] = gdsi_calculate_FS_metrics_Subcortical_vHemis( DSTATS, DATAOUT, aparcaseg, fslabels )

% Calculates metrics across different FS subcortical labels, save to DSTATS
%  struct
% - mean 1D EAP profile, mean 1D GenAniso profile, etc.
%

%%%% brainstem
[ curr ] = mask_with_FS_label( DATAOUT.avg_disp_sph, aparcaseg, fslabels.brainstem );
DSTATS.Rmean.brainstem = mean(curr,1);
DSTATS.Rmean_std.brainstem = std(curr,0,1);

[ curr ] = mask_with_FS_label( DATAOUT.std_disp_sph, aparcaseg, fslabels.brainstem );
DSTATS.Rstd.brainstem = mean(curr,1);

[ curr ] = mask_with_FS_label( (DATAOUT.GaussDiv), aparcaseg, fslabels.brainstem );
DSTATS.GaussDiv.brainstem = mean(curr,1);
DSTATS.GaussDiv_std.brainstem = std(curr,0,1);

[ curr ] = mask_with_FS_label( (DATAOUT.KurtProfile), aparcaseg, fslabels.brainstem );
DSTATS.KurtProfile.brainstem = mean(curr,1);
DSTATS.KurtProfile_std.brainstem = std(curr,0,1);

[ curr ] = mask_with_FS_label( (DATAOUT.MeanKurt), aparcaseg, fslabels.brainstem );
DSTATS.MeanKurt.brainstem = mean(curr,1);
DSTATS.MeanKurt_std.brainstem = std(curr,0,1);



%%%% cc
[ curr ] = mask_with_FS_label( DATAOUT.avg_disp_sph, aparcaseg, fslabels.cc );
DSTATS.Rmean.cc = mean(curr,1);
DSTATS.Rmean_std.cc = std(curr,0,1);

[ curr ] = mask_with_FS_label( DATAOUT.std_disp_sph, aparcaseg, fslabels.cc );
DSTATS.Rstd.cc = mean(curr,1);

[ curr ] = mask_with_FS_label( (DATAOUT.GaussDiv), aparcaseg, fslabels.cc );
DSTATS.GaussDiv.cc = mean(curr,1);
DSTATS.GaussDiv_std.cc = std(curr,0,1);

[ curr ] = mask_with_FS_label( (DATAOUT.KurtProfile), aparcaseg, fslabels.cc );
DSTATS.KurtProfile.cc = mean(curr,1);
DSTATS.KurtProfile_std.cc = std(curr,0,1);

[ curr ] = mask_with_FS_label( (DATAOUT.MeanKurt), aparcaseg, fslabels.cc );
DSTATS.MeanKurt.cc = mean(curr,1);
DSTATS.MeanKurt_std.cc = std(curr,0,1);



%%%% cerebellum
currlabel = fslabels.cerebellum.all;
[ curr ] = mask_with_FS_label( DATAOUT.avg_disp_sph, aparcaseg, currlabel );
DSTATS.Rmean.cerebellum = mean(curr,1);
DSTATS.Rmean_std.cerebellum = std(curr,0,1);

[ curr ] = mask_with_FS_label( DATAOUT.std_disp_sph, aparcaseg, currlabel );
DSTATS.Rstd.cerebellum = mean(curr,1);

[ curr ] = mask_with_FS_label( (DATAOUT.GaussDiv), aparcaseg, currlabel );
DSTATS.GaussDiv.cerebellum = mean(curr,1);
DSTATS.GaussDiv_std.cerebellum = std(curr,0,1);

[ curr ] = mask_with_FS_label( (DATAOUT.KurtProfile), aparcaseg, currlabel );
DSTATS.KurtProfile.cerebellum = mean(curr,1);
DSTATS.KurtProfile_std.cerebellum = std(curr,0,1);

[ curr ] = mask_with_FS_label( (DATAOUT.MeanKurt), aparcaseg, currlabel );
DSTATS.MeanKurt.cerebellum = mean(curr,1);
DSTATS.MeanKurt_std.cerebellum = std(curr,0,1);



%%%% thalamus
currlabel = fslabels.subcort.thalamus.all;
[ curr ] = mask_with_FS_label( DATAOUT.avg_disp_sph, aparcaseg, currlabel );
DSTATS.Rmean.thalamus = mean(curr,1);
DSTATS.Rmean_std.thalamus = std(curr,0,1);

[ curr ] = mask_with_FS_label( DATAOUT.std_disp_sph, aparcaseg, currlabel );
DSTATS.Rstd.thalamus = mean(curr,1);

[ curr ] = mask_with_FS_label( (DATAOUT.GaussDiv), aparcaseg, currlabel );
DSTATS.GaussDiv.thalamus = mean(curr,1);
DSTATS.GaussDiv_std.thalamus = std(curr,0,1);

[ curr ] = mask_with_FS_label( (DATAOUT.KurtProfile), aparcaseg, currlabel );
DSTATS.KurtProfile.thalamus = mean(curr,1);
DSTATS.KurtProfile_std.thalamus = std(curr,0,1);

[ curr ] = mask_with_FS_label( (DATAOUT.MeanKurt), aparcaseg, currlabel );
DSTATS.MeanKurt.thalamus = mean(curr,1);
DSTATS.MeanKurt_std.thalamus = std(curr,0,1);



%%%% putamen
currlabel = fslabels.subcort.putamen.all;
[ curr ] = mask_with_FS_label( DATAOUT.avg_disp_sph, aparcaseg, currlabel );
DSTATS.Rmean.putamen = mean(curr,1);
DSTATS.Rmean_std.putamen = std(curr,0,1);

[ curr ] = mask_with_FS_label( DATAOUT.std_disp_sph, aparcaseg, currlabel );
DSTATS.Rstd.putamen = mean(curr,1);

[ curr ] = mask_with_FS_label( (DATAOUT.GaussDiv), aparcaseg, currlabel );
DSTATS.GaussDiv.putamen = mean(curr,1);
DSTATS.GaussDiv_std.putamen = std(curr,0,1);

[ curr ] = mask_with_FS_label( (DATAOUT.KurtProfile), aparcaseg, currlabel );
DSTATS.KurtProfile.putamen = mean(curr,1);
DSTATS.KurtProfile_std.putamen = std(curr,0,1);

[ curr ] = mask_with_FS_label( (DATAOUT.MeanKurt), aparcaseg, currlabel );
DSTATS.MeanKurt.putamen = mean(curr,1);
DSTATS.MeanKurt_std.putamen = std(curr,0,1);



%%%% caudate
currlabel = fslabels.subcort.caudate.all;
[ curr ] = mask_with_FS_label( DATAOUT.avg_disp_sph, aparcaseg, currlabel );
DSTATS.Rmean.caudate = mean(curr,1);
DSTATS.Rmean_std.caudate = std(curr,0,1);

[ curr ] = mask_with_FS_label( DATAOUT.std_disp_sph, aparcaseg, currlabel );
DSTATS.Rstd.caudate = mean(curr,1);

[ curr ] = mask_with_FS_label( (DATAOUT.GaussDiv), aparcaseg, currlabel );
DSTATS.GaussDiv.caudate = mean(curr,1);
DSTATS.GaussDiv_std.caudate = std(curr,0,1);

[ curr ] = mask_with_FS_label( (DATAOUT.KurtProfile), aparcaseg, currlabel );
DSTATS.KurtProfile.caudate = mean(curr,1);
DSTATS.KurtProfile_std.caudate = std(curr,0,1);

[ curr ] = mask_with_FS_label( (DATAOUT.MeanKurt), aparcaseg, currlabel );
DSTATS.MeanKurt.caudate = mean(curr,1);
DSTATS.MeanKurt_std.caudate = std(curr,0,1);




%%%% pallidum
currlabel = fslabels.subcort.pallidum.all;
[ curr ] = mask_with_FS_label( DATAOUT.avg_disp_sph, aparcaseg, currlabel );
DSTATS.Rmean.pallidum = mean(curr,1);
DSTATS.Rmean_std.pallidum = std(curr,0,1);

[ curr ] = mask_with_FS_label( DATAOUT.std_disp_sph, aparcaseg, currlabel );
DSTATS.Rstd.pallidum = mean(curr,1);

[ curr ] = mask_with_FS_label( (DATAOUT.GaussDiv), aparcaseg, currlabel );
DSTATS.GaussDiv.pallidum = mean(curr,1);
DSTATS.GaussDiv_std.pallidum = std(curr,0,1);

[ curr ] = mask_with_FS_label( (DATAOUT.KurtProfile), aparcaseg, currlabel );
DSTATS.KurtProfile.pallidum = mean(curr,1);
DSTATS.KurtProfile_std.pallidum = std(curr,0,1);

[ curr ] = mask_with_FS_label( (DATAOUT.MeanKurt), aparcaseg, currlabel );
DSTATS.MeanKurt.pallidum = mean(curr,1);
DSTATS.MeanKurt_std.pallidum = std(curr,0,1);



%%%% amygdala
currlabel = fslabels.subcort.amygdala.all;
[ curr ] = mask_with_FS_label( DATAOUT.avg_disp_sph, aparcaseg, currlabel );
DSTATS.Rmean.amygdala = mean(curr,1);
DSTATS.Rmean_std.amygdala = std(curr,0,1);

[ curr ] = mask_with_FS_label( DATAOUT.std_disp_sph, aparcaseg, currlabel );
DSTATS.Rstd.amygdala = mean(curr,1);

[ curr ] = mask_with_FS_label( (DATAOUT.GaussDiv), aparcaseg, currlabel );
DSTATS.GaussDiv.amygdala = mean(curr,1);
DSTATS.GaussDiv_std.amygdala = std(curr,0,1);

[ curr ] = mask_with_FS_label( (DATAOUT.KurtProfile), aparcaseg, currlabel );
DSTATS.KurtProfile.amygdala = mean(curr,1);
DSTATS.KurtProfile_std.amygdala = std(curr,0,1);

[ curr ] = mask_with_FS_label( (DATAOUT.MeanKurt), aparcaseg, currlabel );
DSTATS.MeanKurt.amygdala = mean(curr,1);
DSTATS.MeanKurt_std.amygdala = std(curr,0,1);



%%%% hippocampus
currlabel = fslabels.subcort.hippocampus.all;
[ curr ] = mask_with_FS_label( DATAOUT.avg_disp_sph, aparcaseg, currlabel );
DSTATS.Rmean.hippocampus = mean(curr,1);
DSTATS.Rmean_std.hippocampus = std(curr,0,1);

[ curr ] = mask_with_FS_label( DATAOUT.std_disp_sph, aparcaseg, currlabel );
DSTATS.Rstd.hippocampus = mean(curr,1);

[ curr ] = mask_with_FS_label( (DATAOUT.GaussDiv), aparcaseg, currlabel );
DSTATS.GaussDiv.hippocampus = mean(curr,1);
DSTATS.GaussDiv_std.hippocampus = std(curr,0,1);

[ curr ] = mask_with_FS_label( (DATAOUT.KurtProfile), aparcaseg, currlabel );
DSTATS.KurtProfile.hippocampus = mean(curr,1);
DSTATS.KurtProfile_std.hippocampus = std(curr,0,1);

[ curr ] = mask_with_FS_label( (DATAOUT.MeanKurt), aparcaseg, currlabel );
DSTATS.MeanKurt.hippocampus = mean(curr,1);
DSTATS.MeanKurt_std.hippocampus = std(curr,0,1);



%%%% hippocampus
currlabel = fslabels.subcort.accumbens.all;
[ curr ] = mask_with_FS_label( DATAOUT.avg_disp_sph, aparcaseg, currlabel );
DSTATS.Rmean.accumbens = mean(curr,1);
DSTATS.Rmean_std.accumbens = std(curr,0,1);

[ curr ] = mask_with_FS_label( DATAOUT.std_disp_sph, aparcaseg, currlabel );
DSTATS.Rstd.accumbens = mean(curr,1);

[ curr ] = mask_with_FS_label( (DATAOUT.GaussDiv), aparcaseg, currlabel );
DSTATS.GaussDiv.accumbens = mean(curr,1);
DSTATS.GaussDiv_std.accumbens = std(curr,0,1);

[ curr ] = mask_with_FS_label( (DATAOUT.KurtProfile), aparcaseg, currlabel );
DSTATS.KurtProfile.accumbens = mean(curr,1);
DSTATS.KurtProfile_std.accumbens = std(curr,0,1);

[ curr ] = mask_with_FS_label( (DATAOUT.MeanKurt), aparcaseg, currlabel );
DSTATS.MeanKurt.accumbens = mean(curr,1);
DSTATS.MeanKurt_std.accumbens = std(curr,0,1);




%%%% ventraldc
currlabel = fslabels.subcort.ventraldc.all;
[ curr ] = mask_with_FS_label( DATAOUT.avg_disp_sph, aparcaseg, currlabel );
DSTATS.Rmean.ventraldc = mean(curr,1);
DSTATS.Rmean_std.ventraldc = std(curr,0,1);

[ curr ] = mask_with_FS_label( DATAOUT.std_disp_sph, aparcaseg, currlabel );
DSTATS.Rstd.ventraldc = mean(curr,1);

[ curr ] = mask_with_FS_label( (DATAOUT.GaussDiv), aparcaseg, currlabel );
DSTATS.GaussDiv.ventraldc = mean(curr,1);
DSTATS.GaussDiv_std.ventraldc = std(curr,0,1);

[ curr ] = mask_with_FS_label( (DATAOUT.KurtProfile), aparcaseg, currlabel );
DSTATS.KurtProfile.ventraldc = mean(curr,1);
DSTATS.KurtProfile_std.ventraldc = std(curr,0,1);

[ curr ] = mask_with_FS_label( (DATAOUT.MeanKurt), aparcaseg, currlabel );
DSTATS.MeanKurt.ventraldc = mean(curr,1);
DSTATS.MeanKurt_std.ventraldc = std(curr,0,1);



end
