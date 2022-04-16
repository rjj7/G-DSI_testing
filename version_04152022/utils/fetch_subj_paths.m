function [info] = fetch_subj_paths( subjdir )

%dwi
info.fdwi = [subjdir filesep 'dwi.nii.gz'];
info.flowb = [subjdir filesep 'lowb.nii.gz'];

%brain mask
info.fmask = [subjdir filesep 'nodif_brain_mask.nii.gz'];

%bvec bval
info.fbvec = [subjdir filesep 'bvecs_col_corr'];
info.fbval = [subjdir filesep 'bvals_col_corr'];

%aparc+aseg vol
info.aparcaseg = [subjdir filesep 'aparc+aseg.nii.gz'];

%output dir
info.outdir = [subjdir filesep 'gdsi_recon'];
if ~exist(info.outdir,'dir')
    fprintf('\n --- Making output dir:\n%s\n',info.outdir);
    mkdir(info.outdir);
end

%acq specs file
info.fspecs = [subjdir filesep 'acqspecs.txt'];

%file/proc dir (subj dir)
info.procdir = subjdir;

%dtifit dir
info.dtifitdir = [subjdir filesep 'dtifit'];

%segm mask dir
info.segmdir = [subjdir filesep 'segm'];

%tracula dir
info.trcdir = [subjdir filesep 'trc'];
if ~exist(info.trcdir,'dir')
    info.trcdir = [];
end


end