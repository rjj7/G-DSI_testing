function [ SEGM ] = load_subject_segmasks( INFO )
% [ SEGM ] = load_subject_segmasks( INFO )

if ~isdeployed
    addpath('/autofs/space/hemera_001/users/rjjones/code/toolboxes/NIfTI_20140122');
end


%%%%
SEGM = [];

%%% Load these masks: (made with: mri_binarize --i aparc_aseg.nii.gz ... )
% mask.all-wm.nii.gz
% mask.ctx-wm.nii.gz  
% mask.gm.nii.gz  
% mask.subcort-gm.nii.gz  
% mask.ventricles.nii.gz  
% mask.wm+vcsf.nii.gz

%%%% all-wm
fin = [INFO.segmdir filesep 'mask.all-wm.nii.gz'];
vol = load_untouch_nii(fin);
SEGM.allwm = double(vol.img);

%%%% ctx-wm
fin = [INFO.segmdir filesep 'mask.ctx-wm.nii.gz'];
vol = load_untouch_nii(fin);
SEGM.ctxwm = double(vol.img);

%%%% gm
fin = [INFO.segmdir filesep 'mask.gm.nii.gz'];
vol = load_untouch_nii(fin);
SEGM.gm = double(vol.img);

%%%% subcort-gm
fin = [INFO.segmdir filesep 'mask.subcort-gm.nii.gz'];
vol = load_untouch_nii(fin);
SEGM.subcortgm = double(vol.img);

%%%% ventricles
fin = [INFO.segmdir filesep 'mask.ventricles.nii.gz'];
vol = load_untouch_nii(fin);
SEGM.ventricles = double(vol.img);

%%%% wm+vcsf
fin = [INFO.segmdir filesep 'mask.wm+vcsf.nii.gz'];
vol = load_untouch_nii(fin);
SEGM.wm_vcsf = double(vol.img);


end



