% function [  ] = script_GDSI_analysis_nufft( dtroot, filedir ) %( fconfig, verbose )
% dtroot      = '/autofs/space/nyx_002/users/rjones/Bay8/CS_DSI_10_01_2021';
% % filedir     = [dtroot filesep 'diff/preproc/mri/eddy/results_def/nufft/hcp_bmax7k'];
% filedir     = [dtroot filesep 'diff/preproc/mri/eddy/results_def/nufft/hcp_3sh_bmax_5k'];
%
% 
% GDSI analysis of nufft-shelled data
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



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Mask DATAOUT with FS labels

dtroot      = '/autofs/space/nyx_002/users/rjones/Bay8/CS_DSI_10_01_2021';

% filedir     = [dtroot filesep 'diff/preproc/mri/eddy/results_def/nufft/hcp_bmax7k'];
filedir     = [dtroot filesep 'diff/preproc/mri/eddy/results_def/nufft/hcp_3sh_bmax_5k'];
outdir      = [filedir filesep 'gdsi']; 

fdstats = [outdir filesep 'DSTATS.mat'];
disp(' - Loading DSTATS struct...');
shell3 = load(fdstats,'DSTATS','linespecs','x_disp_vec','INFO','DATAIN','aparcaseg');


filedir     = [dtroot filesep 'diff/preproc/mri/eddy/results_def/nufft/hcp_bmax7k'];
% filedir     = [dtroot filesep 'diff/preproc/mri/eddy/results_def/nufft/hcp_3sh_bmax_5k'];
outdir      = [filedir filesep 'gdsi']; 

fdstats = [outdir filesep 'DSTATS.mat'];
disp(' - Loading DSTATS struct...');
shell4 = load(fdstats,'DSTATS','linespecs','x_disp_vec','INFO','DATAIN','aparcaseg');



filedir     = [dtroot filesep 'diff/preproc/mri/eddy/results_def'];
% filedir     = [dtroot filesep 'diff/preproc/mri/eddy/results_def/nufft/hcp_3sh_bmax_5k'];
outdir      = [filedir filesep 'gdsi'];
fdstats = [outdir filesep 'DSTATS.mat'];
disp(' - Loading DSTATS struct...');
grid = load(fdstats,'DSTATS','linespecs','x_disp_vec','INFO','DATAIN','aparcaseg');


% filedir     = [dtroot filesep 'diff/preproc/mri/eddy/results_def'];
filedir     = [dtroot filesep 'diff/preproc/mri/eddy/results_def/nufft/hcp_2sh_1k,5k'];
outdir      = [filedir filesep 'gdsi'];
fdstats = [outdir filesep 'DSTATS.mat'];
disp(' - Loading DSTATS struct...');
shell1b = load(fdstats,'DSTATS','linespecs','x_disp_vec','INFO','DATAIN','aparcaseg');


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
% ylim([0 0.22]);

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
% ylim([0 0.22]);
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
ylim([-0.0025 ylyl(2)]);
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
% ylim([-0.001 ylyl(2)]);
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






close all;



