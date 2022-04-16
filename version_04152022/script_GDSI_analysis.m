% function [ DATAOUT ] = script_GDSI_recon %( verbose ) %( fconfig, verbose )
% 
% GDSI reconstruction of RSI data
%
% original script by (c) Qiyuan Tian, Stanford RSL
% Modified by RJ, Feb 2021 - Feb 2022
%

% % --- EXAMPLE ---
% % fconfig = '/autofs/space/nyx_002/users/rjones/Bay8/code/gdsi/version_020122/hcp_mgh.config';
% fconfig = '/autofs/space/nyx_002/users/rjones/Bay8/CS_DSI_10_01_2021/analysis/gdsi/dsi.config';
% verbose = true;


%% Setup

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isdeployed
    disp('--RUNNING DEPLOYED APP--');
else
    addpath('/autofs/space/hemera_001/users/rjjones/code/toolboxes/NIfTI_20140122');
    addpath(genpath('/autofs/space/nyx_002/users/rjones/Bay8/code/gdsi/version_021822'));
end


% clear; close all; clc;

% if nargin<1, verbose = true; end

verbose = true;

isshelled = true;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- File definition ---
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Hard-coded variables ---
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%diffusivity of free water in tissue at 37C
%   (this is to set MDD_water, for F_mtx)
D_water = 2.5e-3;

%use spherical vertices from DSI studio (8-fold tesselation)
%   (also used to create F_mtx)
Use_DSI_Studio_dirs     = true;

%if 2+ lowb vols, concat and place in front of 4D dwi stack.
%   (and do the same for bvals/bvecs - 
%    this is needed to create qvecs, then precompute F_mtx)
concatlowb              = true;     %concat + avg lowb volumes


%other parameters
Remove_HighBValue       = false;    %discard high-bvals
ClipRmatrix             = false;    %clip negative values of propagator
Norm2PDF                = false;    %normalize EAP so sum(EAP(:))=1
ReconODFs               = false;    %recon dODFs
ReconZDPs               = false;    %recon zero-displacement maps (RTOP, RTAP, RTPP)
ReconMDs                = false;    %recon mean displacements (1st, 2nd, 4th moments)
MakePlots.odf           = false;    %make plot of ODF
MakePlots.rgbprofiles   = false;    %make direction-encoded plot of all 1D EAP profiles
StoreRmatrix            = false;     %store full EAP in each voxel (ndirs-by-nr mtx)


%   inputs/paths
dtroot      = '/autofs/space/nyx_002/users/rjones/Bay8/CS_DSI_10_01_2021';
filedir     = [dtroot filesep 'diff/preproc/mri/eddy/results_def'];
% filedir     = [dtroot filesep 'diff/preproc/mri/eddy/results_def/nufft/hcp_3sh_bmax5k'];

% fdwi        = [filedir filesep 'data.nii.gz'];
% fmask       = [filedir filesep 'nodif_brain_f0.25_mask.nii.gz'];
% fbvec       = [filedir filesep 'bvecs'];
% fbval       = [filedir filesep 'bvals'];
% flowb       = [filedir filesep 'lowb.nii.gz']; %dtroot filesep 'preproc' filesep 'lowb_eddy_mean.nii.gz'];

outdir      = [filedir filesep 'gdsi']; %[dtroot filesep 'preproc' filesep 'analysis' filesep 'gdsi' filesep 'results-testing'];
if ~exist(outdir,'dir')
    mkdir(outdir); fprintf(' making outdir: %s\n',outdir); end



procdir     = [];
acqspecs    = [dtroot filesep 'analysis' filesep 'gdsi' filesep 'acqspecs.txt'];




big_delta   = 25.5 * 1e-3;    %ms
small_delta = 10.8 * 1e-3;    %ms

%%% Calculate mean diffusion distance of water (micron)
mdd_water = sqrt(6 * D_water * (big_delta-(small_delta/3)) ) * 1000;
% mdd_water = 10;
fprintf('   Diffusivity \t = %g mm^2/s\n', D_water);
fprintf('   MDD water\t = %.2f micron\n', mdd_water);


DATAIN = [];
DATAIN.outdir = outdir;
INFO.outdir = outdir;


fprintf('\n-------------\n - Loading GDSI recon - \n-------------\n');

foutmat = [DATAIN.outdir filesep 'gdsi_vol_recon_prof+scal.mat']
load(foutmat);

foutmat = [DATAIN.outdir filesep 'gdsi_vol_recon_prof+scal_recpar.mat']
load(foutmat);




DATAIN.mask = mask;





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% FREESURFER aparc+aseg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Load FS aparc_aseg labels mapped to diffusion space
anat2diffdir = '/autofs/space/nyx_002/users/rjones/Bay8/CS_DSI_10_01_2021/diff/preproc/mri/eddy/results_def/anat';

INFO.faparcaseg = [anat2diffdir filesep 'aparc+aseg2diff_las.nii.gz'];

% Load aparc+aseg mapped to diffusion space
if ~exist(INFO.faparcaseg,'file')
    disp(' --- COULD NOT FIND APARC+ASEG FILE! ---');
    disp(INFO.faparcaseg);
    return
else
    disp(' - Loading FS aparc+aseg file...');
    temp = load_untouch_nii(INFO.faparcaseg);
    aparcaseg = temp.img;
end




%Zero voxels outside diff data (diffusion data had inferior cut off)
aparcaseg_orig = aparcaseg;
aparcaseg(DATAIN.mask==0)=0;

uthr = prctile(DATAOUT.RTOP(DATAIN.mask==1),97.5);
HyperIntMask = DATAOUT.RTOP>uthr;
aparcaseg(HyperIntMask==1) = 0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Load FS LUT as struct, with fields for nb, name, rgb color
disp(' -Fetching FS aparc+aseg LUT...');
[ fslut ] = load_fs_lut;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Load mapping from labels to tissue/structure type
disp(' -Fetching FS aparc+aseg label numbers...');
[ fslabels ] = fetch_fs_lut_labels;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% 4 plots

[linespecs] = fetch_line_specs;
x_disp_vec = dispvec;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% LOAD ANAT AND DTIFIT VOLS
ft1 = [anat2diffdir filesep 'orig2diff_las.nii.gz'];
ft2 = [anat2diffdir filesep 'T2.prenorm2diff_las.nii.gz'];
ft1n = [anat2diffdir filesep 'norm2diff_las.nii.gz'];
ft2n = [anat2diffdir filesep 'T22diff_las.nii.gz'];

vol=load_untouch_nii(ft1);
t1w = double(vol.img);
vol = load_untouch_nii(ft2);
t2w = double(vol.img);
vol=load_untouch_nii(ft1n);
t1wn = double(vol.img);
vol = load_untouch_nii(ft2n);
t2wn = double(vol.img);


flowb = '/autofs/space/nyx_002/users/rjones/Bay8/CS_DSI_10_01_2021/diff/preproc/mri/eddy/results_def/dmri/lowb.mean.nii.gz';
fdtidir = '/autofs/space/nyx_002/users/rjones/Bay8/CS_DSI_10_01_2021/diff/preproc/mri/eddy/results_def/dtifit';

vol=load_untouch_nii(flowb);
lowb = double(vol.img);

%%%%% ADD function for loading dtifit dir files
vol=load_untouch_nii([fdtidir filesep 'dtifit_FA.nii.gz']);
dti.fa = double(vol.img);
vol=load_untouch_nii([fdtidir filesep 'dtifit_MD.nii.gz']);
dti.md = double(vol.img);
vol=load_untouch_nii([fdtidir filesep 'dtifit_L1.nii.gz']);
dti.l1 = double(vol.img);
vol=load_untouch_nii([fdtidir filesep 'dtifit_L2.nii.gz']);
dti.l2 = double(vol.img);
vol=load_untouch_nii([fdtidir filesep 'dtifit_L3.nii.gz']);
dti.l3 = double(vol.img);
vol=load_untouch_nii([fdtidir filesep 'dtifit_S0.nii.gz']);
dti.s0 = double(vol.img);


DATAOUT.t1w = t1w;
DATAOUT.t2w = t2w;
DATAOUT.t1wn = t1wn;
DATAOUT.t2wn = t2wn;

DATAOUT.lowb = lowb;

DATAOUT.dtifa = dti.fa;
DATAOUT.dtimd = dti.md;
DATAOUT.dtil1 = dti.l1;
DATAOUT.dtil2 = dti.l2;
DATAOUT.dtil3 = dti.l3;
DATAOUT.dtis0 = dti.s0;

    
% VOLS.t1 = t1w;
% VOLS.t2 = t2w;
% VOLS.fa = dtifa;
% VOLS.md = dtimd;
% VOLS.rtop = DATAOUT.RTOP;
% VOLS.ga = DATAOUT.MeanGA;
% VOLS.kurt = DATAOUT.MeanKurt;
% VOLS.gdiv = DATAOUT.GaussDiv;
% VOLS.ng = DATAOUT.sintheta;
% VOLS.hwhm_mean = DATAOUT.HWHM;
% VOLS.hwhm_min = DATAOUT.HWHMmin;
% VOLS.hwhm_max = DATAOUT.HWHMmax;
% VOLS.hwhm_minmax = DATAOUT.HWHMminmax;
% VOLS.hwhm_std = DATAOUT.HWHMstd;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Mask DATAOUT with FS labels

GetDstats = true;
LoadDstats = xor(GetDstats,1);

fdstats = [DATAIN.outdir filesep 'DSTATS.mat'];

if GetDstats
    disp(' - Masking output results with aparc+aseg --SUBCORT-- ...');

    DSTATS = [];

    MetricsList = {'Rvec','Rstd', ...
        'EntropyEAP','MeanGA','RTOP', 'NG','KLDi', 'KLDg',...
        'MeanKurt', 'NLrms', ...
        'HWHM','HWHMmin','HWHMmax','HWHMstd',...
        'qiv','qvar','qentropy',...
        };


    LabelsList = {
        'subcort.thalamus.lh','subcort.thalamus.rh','subcort.thalamus.all', ...
        'subcort.putamen.lh','subcort.putamen.rh','subcort.putamen.all', ...
        'subcort.caudate.lh','subcort.caudate.rh','subcort.caudate.all', ...
        'subcort.pallidum.lh','subcort.pallidum.rh','subcort.pallidum.all', ...
        'subcort.hippocampus.lh','subcort.hippocampus.rh','subcort.hippocampus.all', ...
        'subcort.amygdala.lh','subcort.amygdala.rh','subcort.amygdala.all', ...
        'subcort.accumbens.lh','subcort.accumbens.rh','subcort.accumbens.all', ...
        'subcort.ventraldc.lh','subcort.ventraldc.rh','subcort.ventraldc.all', ...
        'subcort.lh','subcort.rh','subcort.all'};
    % no fornix labels found!!!

    [ DSTATS ] = gdsi_calculate_FS_metrics_Subcortical_v3( DSTATS, DATAOUT, aparcaseg, fslabels, MetricsList, LabelsList );


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Mask DATAOUT with FS labels - for anat + dtifit volumes
    MetricsList = {'t1w','t2w','t1wn','t2wn',...
        'dtifa','dtimd','dtil1','dtil2','dtil3','dtis0'};

    LabelsList = {
        'subcort.thalamus.lh','subcort.thalamus.rh','subcort.thalamus.all', ...
        'subcort.putamen.lh','subcort.putamen.rh','subcort.putamen.all', ...
        'subcort.caudate.lh','subcort.caudate.rh','subcort.caudate.all', ...
        'subcort.pallidum.lh','subcort.pallidum.rh','subcort.pallidum.all', ...
        'subcort.hippocampus.lh','subcort.hippocampus.rh','subcort.hippocampus.all', ...
        'subcort.amygdala.lh','subcort.amygdala.rh','subcort.amygdala.all', ...
        'subcort.accumbens.lh','subcort.accumbens.rh','subcort.accumbens.all', ...
        'subcort.ventraldc.lh','subcort.ventraldc.rh','subcort.ventraldc.all', ...
        'subcort.lh','subcort.rh','subcort.all'};

    [ DSTATS ] = gdsi_calculate_FS_metrics_Subcortical_v3( DSTATS, DATAOUT, aparcaseg, fslabels, MetricsList, LabelsList );

%     save(fdstats,'DSTATS','linespecs','x_disp_vec','INFO','DATAIN');
    


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% FOR TISSUES + OTHER FS-LABELS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp(' - Masking output results with aparc+aseg --TISSUE-- ...');

    MetricsList = {'Rvec','Rstd', ...
        'EntropyEAP','MeanGA','RTOP', 'NG','KLDi', 'KLDg',...
        'MeanKurt', 'NLrms', ...
        'HWHM','HWHMmin','HWHMmax','HWHMstd',...
        'qiv','qvar','qentropy',...
        };

    LabelsList = {'wm.all','wm.lh','wm.rh',...
        'ctx.all','ctx.rh','ctx.lh', ...
        'csf.all','csf.orig',...
        'ventricles.all','ventricles.lh','ventricles.rh', ...
        'cerebellum.wm.lh','cerebellum.wm.rh','cerebellum.wm.all', ...
        'cerebellum.ctx.lh','cerebellum.ctx.rh','cerebellum.ctx.all', ...
        'brainstem','cc', 'cc_post', 'cc_midpost', 'cc_central', 'cc_midant', 'cc_ant', ...
        }; 


    [ DSTATS ] = gdsi_calculate_FS_metrics_Subcortical_v3( DSTATS, DATAOUT, aparcaseg, fslabels, MetricsList, LabelsList );




    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Mask DATAOUT with FS labels - for anat + dtifit volumes

    MetricsList = {'t1w','t2w','t1wn','t2wn',...
        'dtifa','dtimd','dtil1','dtil2','dtil3','dtis0'};

    LabelsList = {'wm.all','wm.lh','wm.rh',...
        'ctx.all','ctx.rh','ctx.lh', ...
        'csf.all','csf.orig',...
        'ventricles.all','ventricles.lh','ventricles.rh', ...
        'cerebellum.wm.lh','cerebellum.wm.rh','cerebellum.wm.all', ...
        'cerebellum.ctx.lh','cerebellum.ctx.rh','cerebellum.ctx.all', ...
        'brainstem','cc', 'cc_post', 'cc_midpost', 'cc_central', 'cc_midant', 'cc_ant', ...
        }; 

    [ DSTATS ] = gdsi_calculate_FS_metrics_Subcortical_v3( DSTATS, DATAOUT, aparcaseg, fslabels, MetricsList, LabelsList );

    
    save(fdstats,'DSTATS','linespecs','x_disp_vec','INFO','DATAIN','aparcaseg');       
    
end



if LoadDstats
    disp(' - Loading DSTATS struct...');
    load(fdstats,'DSTATS','linespecs','x_disp_vec','INFO');
end





% % % % WAITING ON TRACULA RESULTS.........
if 0 == 1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% TRACULA wm tracts
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Load the subjects tracts as logical (binary) arrays
    disp('Loading tracula masks...');
    [ trc ] = load_subject_trc_masks( INFO.FSdir);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Mask GDSI results with trc tract masks
    disp('Masking with tracula bundles masks...');

    [ DSTATS ] = gdsi_calculate_TRC_metrics_bundles_v2( DSTATS, DATAOUT, trc );
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% extract fields from DSTATS struct
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% [ Rmean, Rmean_std, Rstd, Rstd_std, GaussDiv, GaussDiv_std, KurtProfile, ...
%     KurtProfile_std, MeanKurt, MeanKurt_std, MeanEAP, MeanEAP_std, ...
%     MeanGA, MeanGA_std, RTOP, RTOP_std, sintheta, sintheta_std, ...
%     MSD, MSD_std, HWHM, HWHM_std, HWHMmin, HWHMmin_std, ...
%     HWHMmax, HWHMmax_std, HWHMminmax, HWHMminmax_std, HWHMstd, HWHMstd_std ] = ...
%     extract_DSTATS_v3( DSTATS );

fns=fieldnames(DSTATS);
for fni=1:length(fns)
    fnn = fns{fni};
    eval([fnn '=DSTATS.' fnn ';']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Now, do some analysis and make some plots. Fun!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% [linespecs] = fetch_line_specs;
% x_disp_vec = dispvec;

INFO.plotdir = [INFO.outdir filesep 'plots'];
if ~exist(INFO.plotdir,'dir')
    fprintf(' - Making plotdir:\n %s\n',INFO.plotdir);
    mkdir(INFO.plotdir);
end

plotdir = INFO.plotdir;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compare Rmean of WM GM CSF and bundles









%% TAKEN FROM MakePlots4ISMRM.m script



f = figure('position',[1780 320 642 554],'color','w','InvertHardcopy','off');

hold on;
plot(x_disp_vec,Rmean.wm.lh,'Color',linespecs.wm.color,'LineWidth',3,'LineStyle',linespecs.lh.style);
plot(x_disp_vec,Rmean.wm.rh,'Color',linespecs.wm.color,'LineWidth',3,'LineStyle',linespecs.rh.style);

plot(x_disp_vec,Rmean.ctx.lh,'Color',linespecs.ctx.color,'LineWidth',3,'LineStyle',linespecs.lh.style);
plot(x_disp_vec,Rmean.ctx.rh,'Color',linespecs.ctx.color,'LineWidth',3,'LineStyle',linespecs.rh.style);

plot(x_disp_vec,Rmean.ventricles.lh,'Color',linespecs.csf.color,'LineWidth',3,'LineStyle',linespecs.lh.style);
plot(x_disp_vec,Rmean.ventricles.rh,'Color',linespecs.csf.color,'LineWidth',3,'LineStyle',linespecs.rh.style);

legend('WM lh','WM rh','GM-ctx lh','GM-ctx rh','Vent. lh','Vent. rh');
title('Mean 1D EAP profile');
ylabel('Mean probability density');
xlabel('Displacement (micron)');
set(gca,'FontSize',15);

xlim([0 max(x_disp_vec)]);

fout = [plotdir filesep 'Rvec_3tissue_lineplot.png'];
print(f,fout,'-dpng');
% 
% 
% 
% 
f = figure('position',[1780 320 642 554],'color','w','InvertHardcopy','off');

hold on;
plot(x_disp_vec,Rmean.wm.lh,'Color',linespecs.wm.color,'LineWidth',3,'LineStyle',linespecs.lh.style);
plot(x_disp_vec,Rmean.wm.rh,'Color',linespecs.wm.color,'LineWidth',3,'LineStyle',linespecs.rh.style);

plot(x_disp_vec,Rmean.subcort.lh,'Color','m','LineWidth',3,'LineStyle',linespecs.lh.style);
plot(x_disp_vec,Rmean.subcort.rh,'Color','m','LineWidth',3,'LineStyle',linespecs.rh.style);

plot(x_disp_vec,Rmean.ctx.lh,'Color',linespecs.ctx.color,'LineWidth',3,'LineStyle',linespecs.lh.style);
plot(x_disp_vec,Rmean.ctx.rh,'Color',linespecs.ctx.color,'LineWidth',3,'LineStyle',linespecs.rh.style);

plot(x_disp_vec,Rmean.ventricles.lh,'Color',linespecs.csf.color,'LineWidth',3,'LineStyle',linespecs.lh.style);
plot(x_disp_vec,Rmean.ventricles.rh,'Color',linespecs.csf.color,'LineWidth',3,'LineStyle',linespecs.rh.style);

legend('WM lh','WM rh','GM-subctx lh','GM-subctx rh','GM-ctx lh','GM-ctx rh','Vent. lh','Vent. rh');
title('Mean 1D EAP profile');
ylabel('Mean probability density');
xlabel('Displacement (micron)');
set(gca,'FontSize',15);

xlim([0 max(x_disp_vec)]);

fout = [plotdir filesep 'Rvec_4tissue_lineplot.png'];
print(f,fout,'-dpng');

fout = [plotdir filesep 'Rvec_4tissue_lineplot_hi-res.png'];
print(f,fout,'-dpng','-r300');





f = figure('position',[1780 320 642 554],'color','w','InvertHardcopy','off');

hold on;
errorbar(x_disp_vec,Rmean.wm.all,Rmean_std.wm.all,'Color',linespecs.wm.color,'LineWidth',linespecs.all.width,'LineStyle',linespecs.all.style);

errorbar(x_disp_vec,Rmean.ctx.all,Rmean_std.ctx.all,'Color',linespecs.ctx.color,'LineWidth',linespecs.all.width,'LineStyle',linespecs.all.style);

errorbar(x_disp_vec,Rmean.ventricles.all,Rmean_std.ventricles.all,'Color',linespecs.csf.color,'LineWidth',linespecs.all.width,'LineStyle',linespecs.all.style);

plot(x_disp_vec,zeros(size(x_disp_vec)),'-k');

legend('WM','GM-ctx','CSF');
title('Mean 1D EAP profile');
ylabel('Mean probability density');
xlabel('Displacement (micron)');
set(gca,'FontSize',15);

xlim([0 max(x_disp_vec)]);
ylyl=ylim;
ylim([-0.005 ylyl(2)]);

fout = [plotdir filesep 'Rvec_3tissue_errorbars.png'];
print(f,fout,'-dpng');



f = figure('position',[1780 320 642 554],'color','w','InvertHardcopy','off');

hold on;
errorbar(x_disp_vec,Rmean.wm.all,Rmean_std.wm.all,'Color',linespecs.wm.color,'LineWidth',linespecs.all.width,'LineStyle',linespecs.all.style);

errorbar(x_disp_vec,Rmean.subcort.all,Rmean_std.subcort.all,'Color','m','LineWidth',linespecs.all.width,'LineStyle',linespecs.all.style);

errorbar(x_disp_vec,Rmean.ctx.all,Rmean_std.ctx.all,'Color',linespecs.ctx.color,'LineWidth',linespecs.all.width,'LineStyle',linespecs.all.style);

errorbar(x_disp_vec,Rmean.ventricles.all,Rmean_std.ventricles.all,'Color',linespecs.csf.color,'LineWidth',linespecs.all.width,'LineStyle',linespecs.all.style);

plot(x_disp_vec,zeros(size(x_disp_vec)),'-k');

legend('WM','GM-subctx','GM-ctx','CSF');
title('Mean 1D EAP profile');
ylabel('Mean probability density');
xlabel('Displacement (micron)');
set(gca,'FontSize',15);

xlim([0 max(x_disp_vec)]);
ylyl=ylim;
ylim([-0.005 ylyl(2)]);

% 
fout = [plotdir filesep 'Rvec_4tissue_errorbars.png'];
print(f,fout,'-dpng');






linespecs.hippocampus.color = [0.1 0.8 0.1];
linespecs.putamen.color = 0.6*linespecs.putamen.color;

% f = figure('position',[1780 320 642 554],'color','w','InvertHardcopy','off');
% f = figure('position',[322 47 1092 903],'color','w','InvertHardcopy','off');
f = figure('position',[1728 360 743 675],'color','w','InvertHardcopy','off');

hold on;
plot(x_disp_vec,Rmean.subcort.putamen.lh,'Color',linespecs.putamen.color,'LineWidth',linespecs.lh.width,'LineStyle',linespecs.lh.style);
plot(x_disp_vec,Rmean.subcort.putamen.rh,'Color',linespecs.putamen.color,'LineWidth',linespecs.rh.width,'LineStyle',linespecs.rh.style);

plot(x_disp_vec,Rmean.subcort.thalamus.lh,'Color',linespecs.thalamus.color,'LineWidth',linespecs.lh.width,'LineStyle',linespecs.lh.style);
plot(x_disp_vec,Rmean.subcort.thalamus.rh,'Color',linespecs.thalamus.color,'LineWidth',linespecs.rh.width,'LineStyle',linespecs.rh.style);

plot(x_disp_vec,Rmean.subcort.caudate.lh,'Color',linespecs.caudate.color,'LineWidth',linespecs.lh.width,'LineStyle',linespecs.lh.style);
plot(x_disp_vec,Rmean.subcort.caudate.rh,'Color',linespecs.caudate.color,'LineWidth',linespecs.rh.width,'LineStyle',linespecs.rh.style);

plot(x_disp_vec,Rmean.subcort.pallidum.lh,'Color',linespecs.pallidum.color,'LineWidth',linespecs.lh.width,'LineStyle',linespecs.lh.style);
plot(x_disp_vec,Rmean.subcort.pallidum.rh,'Color',linespecs.pallidum.color,'LineWidth',linespecs.rh.width,'LineStyle',linespecs.rh.style);

plot(x_disp_vec,Rmean.subcort.amygdala.lh,'Color',linespecs.amygdala.color,'LineWidth',linespecs.lh.width,'LineStyle',linespecs.lh.style);
plot(x_disp_vec,Rmean.subcort.amygdala.rh,'Color',linespecs.amygdala.color,'LineWidth',linespecs.rh.width,'LineStyle',linespecs.rh.style);

plot(x_disp_vec,Rmean.subcort.hippocampus.lh,'Color',linespecs.hippocampus.color,'LineWidth',linespecs.lh.width,'LineStyle',linespecs.lh.style);
plot(x_disp_vec,Rmean.subcort.hippocampus.rh,'Color',linespecs.hippocampus.color,'LineWidth',linespecs.rh.width,'LineStyle',linespecs.rh.style);

plot(x_disp_vec,Rmean.subcort.accumbens.lh,'Color',linespecs.accumbens.color,'LineWidth',linespecs.lh.width,'LineStyle',linespecs.lh.style);
plot(x_disp_vec,Rmean.subcort.accumbens.rh,'Color',linespecs.accumbens.color,'LineWidth',linespecs.rh.width,'LineStyle',linespecs.rh.style);

plot(x_disp_vec,Rmean.subcort.ventraldc.lh,'Color',linespecs.ventraldc.color,'LineWidth',linespecs.lh.width,'LineStyle',linespecs.lh.style);
plot(x_disp_vec,Rmean.subcort.ventraldc.rh,'Color',linespecs.ventraldc.color,'LineWidth',linespecs.rh.width,'LineStyle',linespecs.rh.style);

plot([min(x_disp_vec) max(x_disp_vec)],[0 0],'-k');

legend('put lh','put rh','thal lh','thal rh','caud lh','caud rh','pall lh','pall rh',...
    'amyg lh','amyg rh','hipp lh','hipp rh','accmb lh','accmb rh','ventdc lh','ventdc rh');

title('Subcortical - Mean 1D EAP profiles');
ylabel('Mean probability density');
xlabel('Displacement (micron)');
set(gca,'FontSize',15);

xlim([min(x_disp_vec) max(x_disp_vec)]);
ylyl=ylim;
ylim([-0.005 ylyl(2)]);
whyes = ylim;

fplot = [plotdir filesep 'Rmean_subcort_lh+rh.png'];
print(f,fplot,'-dpng');

fplot = [plotdir filesep 'Rmean_subcort_lh+rh_hi-res.png'];
print(f,fplot,'-dpng','-r300');


plot([x_disp_vec(1) x_disp_vec(1)],whyes,':r','LineWidth',2);
plot([x_disp_vec(31) x_disp_vec(31)],whyes,':r','LineWidth',2);
plot([x_disp_vec(46) x_disp_vec(46)],whyes,':r','LineWidth',2);

legend('put lh','put rh','thal lh','thal rh','caud lh','caud rh','pall lh','pall rh',...
    'amyg lh','amyg rh','hipp lh','hipp rh','accmb lh','accmb rh','ventdc lh','ventdc rh');


fplot = [plotdir filesep 'Rmean_subcort_lh+rh_withVertLines.png'];
print(f,fplot,'-dpng');

fplot = [plotdir filesep 'Rmean_subcort_lh+rh_withVertLines_hi-res.png'];
print(f,fplot,'-dpng','-r300');


clf;











f = figure('position',[1780 320 642 554],'color','w','InvertHardcopy','off');

hold on;
plot(x_disp_vec,Rstd.wm.lh,'Color',linespecs.wm.color,'LineWidth',3,'LineStyle',linespecs.lh.style);
plot(x_disp_vec,Rstd.wm.rh,'Color',linespecs.wm.color,'LineWidth',3,'LineStyle',linespecs.rh.style);

plot(x_disp_vec,Rstd.ctx.lh,'Color',linespecs.ctx.color,'LineWidth',3,'LineStyle',linespecs.lh.style);
plot(x_disp_vec,Rstd.ctx.rh,'Color',linespecs.ctx.color,'LineWidth',3,'LineStyle',linespecs.rh.style);

plot(x_disp_vec,Rstd.ventricles.lh,'Color',linespecs.csf.color,'LineWidth',3,'LineStyle',linespecs.lh.style);
plot(x_disp_vec,Rstd.ventricles.rh,'Color',linespecs.csf.color,'LineWidth',3,'LineStyle',linespecs.rh.style);

legend('WM lh','WM rh','GM lh','GM rh','Vent. lh','Vent. rh');
title('Mean 1D Gen. Aniso. profile');
ylabel('Generalized Anisotropy ');
xlabel('Displacement (micron)');
set(gca,'FontSize',15);

xlim([0 max(x_disp_vec)]);

fout = [plotdir filesep 'Rstd_3tissue_lineplot.png'];
print(f,fout,'-dpng');
% 
% 
% 
% 
f = figure('position',[1780 320 642 554],'color','w','InvertHardcopy','off');

hold on;
plot(x_disp_vec,Rstd.wm.lh,'Color',linespecs.wm.color,'LineWidth',3,'LineStyle',linespecs.lh.style);
plot(x_disp_vec,Rstd.wm.rh,'Color',linespecs.wm.color,'LineWidth',3,'LineStyle',linespecs.rh.style);

plot(x_disp_vec,Rstd.subcort.lh,'Color','m','LineWidth',3,'LineStyle',linespecs.lh.style);
plot(x_disp_vec,Rstd.subcort.rh,'Color','m','LineWidth',3,'LineStyle',linespecs.rh.style);

plot(x_disp_vec,Rstd.ctx.lh,'Color',linespecs.ctx.color,'LineWidth',3,'LineStyle',linespecs.lh.style);
plot(x_disp_vec,Rstd.ctx.rh,'Color',linespecs.ctx.color,'LineWidth',3,'LineStyle',linespecs.rh.style);

plot(x_disp_vec,Rstd.ventricles.lh,'Color',linespecs.csf.color,'LineWidth',3,'LineStyle',linespecs.lh.style);
plot(x_disp_vec,Rstd.ventricles.rh,'Color',linespecs.csf.color,'LineWidth',3,'LineStyle',linespecs.rh.style);

legend('WM lh','WM rh','GM-subctx lh','GM-subctx rh','GM-ctx lh','GM-ctx rh','Vent. lh','Vent. rh');
title('Mean 1D Gen. Aniso. profile');
ylabel('Generalized Anisotropy ');
xlabel('Displacement (micron)');
set(gca,'FontSize',15);

xlim([0 max(x_disp_vec)]);

fout = [plotdir filesep 'Rstd_4tissue_lineplot.png'];
print(f,fout,'-dpng');

fout = [plotdir filesep 'Rstd_4tissue_lineplot_hi-res.png'];
print(f,fout,'-dpng','-r300');
% 
% 
% 
% 
f = figure('position',[1780 320 642 554],'color','w','InvertHardcopy','off');

hold on;
errorbar(x_disp_vec,Rstd.wm.all,Rstd_std.wm.all,'Color',linespecs.wm.color,'LineWidth',linespecs.all.width,'LineStyle',linespecs.all.style);

errorbar(x_disp_vec,Rstd.ctx.all,Rstd_std.ctx.all,'Color',linespecs.ctx.color,'LineWidth',linespecs.all.width,'LineStyle',linespecs.all.style);

errorbar(x_disp_vec,Rstd.ventricles.all,Rstd_std.ventricles.all,'Color',linespecs.csf.color,'LineWidth',linespecs.all.width,'LineStyle',linespecs.all.style);

plot(x_disp_vec,zeros(size(x_disp_vec)),'-k');

legend('WM','GM','CSF');
title('Mean 1D Gen. Aniso. profile');
ylabel('Generalized Anisotropy ');
xlabel('Displacement (micron)');
set(gca,'FontSize',15);

xlim([0 max(x_disp_vec)]);
ylyl=ylim;
ylim([-0.001 ylyl(2)]);
% 
fout = [plotdir filesep 'Rstd_3tissue_errorbars.png'];
print(f,fout,'-dpng');





% f = figure('position',[1780 320 642 554],'color','w','InvertHardcopy','off');
f = figure('position',[1728 360 743 675],'color','w','InvertHardcopy','off');

hold on;
plot(x_disp_vec,Rstd.subcort.putamen.lh,'Color',linespecs.putamen.color,'LineWidth',linespecs.lh.width,'LineStyle',linespecs.lh.style);
plot(x_disp_vec,Rstd.subcort.putamen.rh,'Color',linespecs.putamen.color,'LineWidth',linespecs.rh.width,'LineStyle',linespecs.rh.style);

plot(x_disp_vec,Rstd.subcort.thalamus.lh,'Color',linespecs.thalamus.color,'LineWidth',linespecs.lh.width,'LineStyle',linespecs.lh.style);
plot(x_disp_vec,Rstd.subcort.thalamus.rh,'Color',linespecs.thalamus.color,'LineWidth',linespecs.rh.width,'LineStyle',linespecs.rh.style);

plot(x_disp_vec,Rstd.subcort.caudate.lh,'Color',linespecs.caudate.color,'LineWidth',linespecs.lh.width,'LineStyle',linespecs.lh.style);
plot(x_disp_vec,Rstd.subcort.caudate.rh,'Color',linespecs.caudate.color,'LineWidth',linespecs.rh.width,'LineStyle',linespecs.rh.style);

plot(x_disp_vec,Rstd.subcort.pallidum.lh,'Color',linespecs.pallidum.color,'LineWidth',linespecs.lh.width,'LineStyle',linespecs.lh.style);
plot(x_disp_vec,Rstd.subcort.pallidum.rh,'Color',linespecs.pallidum.color,'LineWidth',linespecs.rh.width,'LineStyle',linespecs.rh.style);

plot(x_disp_vec,Rstd.subcort.amygdala.lh,'Color',linespecs.amygdala.color,'LineWidth',linespecs.lh.width,'LineStyle',linespecs.lh.style);
plot(x_disp_vec,Rstd.subcort.amygdala.rh,'Color',linespecs.amygdala.color,'LineWidth',linespecs.rh.width,'LineStyle',linespecs.rh.style);

plot(x_disp_vec,Rstd.subcort.hippocampus.lh,'Color',linespecs.hippocampus.color,'LineWidth',linespecs.lh.width,'LineStyle',linespecs.lh.style);
plot(x_disp_vec,Rstd.subcort.hippocampus.rh,'Color',linespecs.hippocampus.color,'LineWidth',linespecs.rh.width,'LineStyle',linespecs.rh.style);

plot(x_disp_vec,Rstd.subcort.accumbens.lh,'Color',linespecs.accumbens.color,'LineWidth',linespecs.lh.width,'LineStyle',linespecs.lh.style);
plot(x_disp_vec,Rstd.subcort.accumbens.rh,'Color',linespecs.accumbens.color,'LineWidth',linespecs.rh.width,'LineStyle',linespecs.rh.style);

plot(x_disp_vec,Rstd.subcort.ventraldc.lh,'Color',linespecs.ventraldc.color,'LineWidth',linespecs.lh.width,'LineStyle',linespecs.lh.style);
plot(x_disp_vec,Rstd.subcort.ventraldc.rh,'Color',linespecs.ventraldc.color,'LineWidth',linespecs.rh.width,'LineStyle',linespecs.rh.style);

% plot([min(x_disp_vec) max(x_disp_vec)],[0 0],'-k');

legend('put lh','put rh','thal lh','thal rh','caud lh','caud rh','pall lh','pall rh',...
    'amyg lh','amyg rh','hipp lh','hipp rh','accmb lh','accmb rh','ventdc lh','ventdc rh');

title('Subcortical - Mean 1D Gen. Aniso. profile');
ylabel('Generalized Anisotropy');
xlabel('Displacement (micron)');
set(gca,'FontSize',15);

xlim([min(x_disp_vec) max(x_disp_vec)]);
% ylim([-0.01 0.25]);
whyes = ylim;

fplot = [plotdir filesep 'Rstd_subcort_lh+rh.png'];
print(f,fplot,'-dpng');

fplot = [plotdir filesep 'Rstd_subcort_lh+rh_hi-res.png'];
print(f,fplot,'-dpng','-r300');

plot([x_disp_vec(2) x_disp_vec(2)],whyes,':r','LineWidth',2);
plot([x_disp_vec(31) x_disp_vec(31)],whyes,':r','LineWidth',2);
plot([x_disp_vec(46) x_disp_vec(46)],whyes,':r','LineWidth',2);

legend('put lh','put rh','thal lh','thal rh','caud lh','caud rh','pall lh','pall rh',...
    'amyg lh','amyg rh','hipp lh','hipp rh','accmb lh','accmb rh','ventdc lh','ventdc rh');


fplot = [plotdir filesep 'Rstd_subcort_lh+rh_withVertLines.png'];
print(f,fplot,'-dpng');
fplot = [plotdir filesep 'Rstd_subcort_lh+rh_withVertLines_hi-res.png'];
print(f,fplot,'-dpng','-r300');


clf;










