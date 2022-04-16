function [  ] = gdsi_save_nii_vols( DATAOUT, outniidir, vol, verbose )

if nargin<3
    disp('USAGE: ');
    disp(' [  ] = gdsi_save_nii_vols( DATAOUT, outniidir, vol[, verbose] )');
    return
end
if nagrin<4
    verbose = 0;
end

if verbose, disp(' --- Saving nifti volumes --- '); end

if verbose, disp(' -saving rtop'); end
vol.img = DATAOUT.scalarmaps.rtop1;
foutnii = [outniidir filesep 'rtop.nii.gz'];
save_untouch_nii(vol,foutnii);

if verbose, disp(' -saving mean kurtosis'); end
vol.img = DATAOUT.scalarmaps.kurtosis;
foutnii = [outniidir filesep 'kurtosis.nii.gz'];
save_untouch_nii(vol,foutnii);

if verbose, disp(' -saving mean KL div'); end
vol.img = DATAOUT.scalarmaps.isokldiv;
foutnii = [outniidir filesep 'isokldiv.nii.gz'];
save_untouch_nii(vol,foutnii);

if verbose, disp(' -saving mean GFA'); end
vol.img = DATAOUT.scalarmaps.std_disp_sph_mean;
foutnii = [outniidir filesep 'mean_gfa.nii.gz'];
save_untouch_nii(vol,foutnii);

if verbose, disp(' -saving gra profile'); end
bigvol = vol;
bigvol.hdr.dime.dim([1 5]) = [4 100];
bigvol.img = DATAOUT.std_disp_sph;
foutnii = [outniidir filesep 'profile_gfa.nii.gz'];
save_untouch_nii(bigvol,foutnii);

if verbose, disp(' -saving mean eap profile'); end
bigvol.img = DATAOUT.avg_disp_sph;
foutnii = [outniidir filesep 'profile_eap.nii.gz'];
save_untouch_nii(bigvol,foutnii);

if verbose, disp(' -saving KL div profile'); end
bigvol.img = DATAOUT.isokldiv;
foutnii = [outniidir filesep 'profile_isokldiv.nii.gz'];
save_untouch_nii(bigvol,foutnii);

end