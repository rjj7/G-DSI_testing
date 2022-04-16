function [ fslabels ] = fetch_fs_lut_labels
%
% [ fslabels ] = fetch_fs_lut_labels
%
% Create cell array with FS label numbers for various tissues/structures.
%
% NOTE: 
%   - .all is .lh+.rh, i couldn't use + in field name
%   - .bh (ventricles) is "both hemis"; two vent.s arent assigned lh or rh
% 
% Labels in fslabels:
% 
% .wm (.lh, .rh, .all)
% .ctx (.lh, .rh, .all)
% .ventricles (.lh, .rh, .bh, .all)
% .csf (.orig = single LUT CSF label; .all = .csf.orig + .ventricles.all)
% .cerebellum (.wm.lh, .wm.rh, .wm.all, .ctx.lh, .ctx.rh, .ctx.all, .all)
% .subcort (.thalamus, .putamen, .caudate, .pallidum, .hippocampus,
%           .amygdala, .accumbens, .ventraldc; all have .lh, .rh, .all)
% .brainstem
% .cc
% 
% 


%%% A little info
fslabels = {};
fslabels.info.name = 'FS LUT aparc+aseg';
fslabels.info.timestamp = datetime;
[a,b]=system('echo $FREESURFER');
if a==0, fslabels.info.version = b(1:end-1); end
%%% Some notes on notation
fslabels.info.readme.lhrh = '.lh+.rh are left/right hemis';
fslabels.info.readme.ventricles = '.bh = ventricles without lh/rh label';
fslabels.info.readme.csf = '.orig = FS LUT CSF label; .all = .ventricles.all + .csf.orig';
fslabels.info.readme.subcort = '.lh, .rh, .all for all structures';
fslabels.info.readme.cerebellum = '.lh, .rh, .all; also all 3 for .wm and .ctx';


%%% Cerebelar white matter:
fslabels.wm.lh = 2;
fslabels.wm.rh = 41;
fslabels.wm.all = cat(2,fslabels.wm.lh,fslabels.wm.rh);


%%% Cerebelar cortex:
fslabels.ctx.lh = 1000:1035;
fslabels.ctx.rh = 2000:2035;
fslabels.ctx.all = cat(2,fslabels.ctx.lh,fslabels.ctx.rh);


%%% Ventricles + CSF
fslabels.ventricles.lh = [4,5];
fslabels.ventricles.rh = [43,44];
fslabels.ventricles.bh = [14,15];
fslabels.ventricles.all = cat(2,fslabels.ventricles.lh,...
    fslabels.ventricles.rh,fslabels.ventricles.bh);

fslabels.csf.orig = 24;
fslabels.csf.all = cat(2,fslabels.ventricles.all,fslabels.csf.orig);


%%% Cerebellum:
fslabels.cerebellum.wm.lh = 7;
fslabels.cerebellum.wm.rh = 46;
fslabels.cerebellum.wm.all = cat(2,fslabels.cerebellum.wm.lh,...
    fslabels.cerebellum.wm.rh);
fslabels.cerebellum.ctx.lh = 8;
fslabels.cerebellum.ctx.rh = 47;
fslabels.cerebellum.ctx.all = cat(2,fslabels.cerebellum.ctx.lh,...
    fslabels.cerebellum.ctx.rh);
fslabels.cerebellum.all = cat(2,fslabels.cerebellum.wm.all,...
    fslabels.cerebellum.ctx.all);
fslabels.cerebellum.lh.all = cat(2,fslabels.cerebellum.wm.lh,...
    fslabels.cerebellum.ctx.lh);
fslabels.cerebellum.rh.all = cat(2,fslabels.cerebellum.wm.rh,...
    fslabels.cerebellum.ctx.rh);


%%% Subcortical structures/nuclei:
fslabels.subcort.thalamus.lh = 10;
fslabels.subcort.thalamus.rh = 49;
fslabels.subcort.thalamus.all = cat(2,fslabels.subcort.thalamus.lh,...
    fslabels.subcort.thalamus.rh);
fslabels.subcort.putamen.lh = 12;
fslabels.subcort.putamen.rh = 51;
fslabels.subcort.putamen.all = cat(2,fslabels.subcort.putamen.lh,...
    fslabels.subcort.putamen.rh);
fslabels.subcort.caudate.lh = 11;
fslabels.subcort.caudate.rh = 50;
fslabels.subcort.caudate.all = cat(2,fslabels.subcort.caudate.lh,...
    fslabels.subcort.caudate.rh);
fslabels.subcort.pallidum.lh = 13;
fslabels.subcort.pallidum.rh = 52;
fslabels.subcort.pallidum.all = cat(2,fslabels.subcort.pallidum.lh,...
    fslabels.subcort.pallidum.rh);
fslabels.subcort.hippocampus.lh = 17;
fslabels.subcort.hippocampus.rh = 53;
fslabels.subcort.hippocampus.all = cat(2,fslabels.subcort.hippocampus.lh,...
    fslabels.subcort.hippocampus.rh);
fslabels.subcort.amygdala.lh = 18;
fslabels.subcort.amygdala.rh = 54;
fslabels.subcort.amygdala.all = cat(2,fslabels.subcort.amygdala.lh,...
    fslabels.subcort.amygdala.rh);
fslabels.subcort.accumbens.lh = 26;
fslabels.subcort.accumbens.rh = 58;
fslabels.subcort.accumbens.all = cat(2,fslabels.subcort.accumbens.lh,...
    fslabels.subcort.accumbens.rh);
fslabels.subcort.ventraldc.lh = 28;
fslabels.subcort.ventraldc.rh = 60;
fslabels.subcort.ventraldc.all = cat(2,fslabels.subcort.ventraldc.lh,...
    fslabels.subcort.ventraldc.rh);

fslabels.subcort.lh = cat(2,fslabels.subcort.thalamus.lh,...
    fslabels.subcort.putamen.lh,fslabels.subcort.caudate.lh,...
    fslabels.subcort.pallidum.lh,fslabels.subcort.hippocampus.lh,...
    fslabels.subcort.amygdala.lh,fslabels.subcort.accumbens.lh,...
    fslabels.subcort.ventraldc.lh);
fslabels.subcort.rh = cat(2,fslabels.subcort.thalamus.rh,...
    fslabels.subcort.putamen.rh,fslabels.subcort.caudate.rh,...
    fslabels.subcort.pallidum.rh,fslabels.subcort.hippocampus.rh,...
    fslabels.subcort.amygdala.rh,fslabels.subcort.accumbens.rh,...
    fslabels.subcort.ventraldc.rh);
fslabels.subcort.all = cat(2,fslabels.subcort.lh,fslabels.subcort.rh);


%%% Other labels that might be of interest:
fslabels.fornix = 250;
fslabels.brainstem = 15;
fslabels.cc = 251:255;
fslabels.cc_post = 251;
fslabels.cc_midpost = 252;
fslabels.cc_central = 253;
fslabels.cc_midant = 254;
fslabels.cc_ant = 255;


end
