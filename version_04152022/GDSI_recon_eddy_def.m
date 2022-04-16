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

isshelled = false;


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
filedir     = [dtroot filesep 'diff/preproc/mri/eddy/results_def/dmri'];

fdwi        = [filedir filesep 'data.nii'];
fmask       = [filedir filesep 'nodif_brain_f0.25_mask.nii.gz'];
fbvec       = [filedir filesep 'bvecs'];
fbval       = [filedir filesep 'bvals'];
flowb       = [filedir filesep 'lowb.mean.nii.gz']; %dtroot filesep 'preproc' filesep 'lowb_eddy_mean.nii.gz'];

outdir      = [dtroot filesep 'diff/preproc/mri/eddy/results_def/gdsi']; %[dtroot filesep 'preproc' filesep 'analysis' filesep 'gdsi' filesep 'results-testing'];
if ~exist(outdir,'dir')
    mkdir(outdir); fprintf(' making outdir: %s\n',outdir); end


procdir     = [];
acqspecs    = [dtroot filesep 'analysis' filesep 'gdsi' filesep 'acqspecs.txt'];




big_delta   = 25.5 * 1e-3;    %ms
small_delta = 10.8 * 1e-3;    %ms




%% Load data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Load subject information
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if verbose
%     disp('Loading config file...');
% end
% [ INFO ] = read_cfg_file( fconfig );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Load acquisition information
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if verbose
%     disp('Loading acq specs file...');
% end
% [ ACQ ] = read_acqspecs_file( INFO.fspecs );


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Load diffusion data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if verbose
    disp('Loading subject data...');
end
% [ DATAIN ] = load_subject_data( INFO, verbose ); %opts.loaddwi=true
DATAIN = [];

%%%% DWI (4D VOL)
if verbose
    disp('Loading dwi...');
end
vol = load_untouch_nii(fdwi);
DATAIN.dwi = double(vol.img);
% DATAIN.lowb0 = DATAIN.dwi(:,:,:,1);

%%%% BRAIN MASK (3D VOL)
if verbose
    disp('Loading mask...');
end
vol = load_untouch_nii(fmask);
DATAIN.mask = double(vol.img);

%%%% BVEC/BVAL
if verbose
    disp('Loading bvecs/bvals...');
end
DATAIN.bvecs = dlmread(fbvec);
DATAIN.bvals = dlmread(fbval);

% Check that dims of bvecs/bvals/dwi match up
validate_diff_data( DATAIN );
if verbose
    disp(' -- Diffusion data validated --');
end

%%%% Find b=0 indices
index_b0 = DATAIN.bvals==0;

%%%% MEAN LOWB VOLUME
DATAIN.lowb = mean(DATAIN.dwi(:,:,:,index_b0),4);

%%% Remove duplicate lowb volumes
if verbose
    disp('Creating qvec...');
end
if size(DATAIN.bvals,2)>size(DATAIN.bvals,1)
    DATAIN.bvals = DATAIN.bvals';
end
if size(DATAIN.bvecs,2)>size(DATAIN.bvecs,1)
    DATAIN.bvecs = DATAIN.bvecs';
end
if concatlowb && nnz(index_b0)>1
    disp(' Concatenating b=0s in bval, bvec...');
    DATAIN.bvals_orig = DATAIN.bvals;
    DATAIN.bvecs_orig = DATAIN.bvecs;
    DATAIN.bvals = DATAIN.bvals(~index_b0);
    DATAIN.bvals = cat(1,0,DATAIN.bvals);
    DATAIN.bvecs = DATAIN.bvecs(~index_b0,:);
    DATAIN.bvecs = cat(1,[0 0 0],DATAIN.bvecs);
end

%%% Create qvec
qvec = repmat(sqrt( 6 * D_water * DATAIN.bvals ), [1, 3]) .* DATAIN.bvecs;
nqvec = round(qvec * (5/round(max(qvec(:)))));
mag_qvec = vecnorm(nqvec,2,2);

DATAIN.qvec = qvec;

DATAIN.outdir = outdir;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Prepare for GDSI recon
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%
%%%%% NOTE: 
%%%%%   Orginally not sure what Delta/delta were, so used values from Fan
%%%%%   et al. (i think, or related bay 8 paper). 
%%%%%    - was using Delta/delta = 21.8/12.9 ms --> MDDwater = 16.20 um
%%%%%    - actually  Delta/delta = 25.5/10.8 ms --> MDDwater = 18.12 um    
%%%%%

%%% Calculate mean diffusion distance of water (micron)
mdd_water = sqrt(6 * D_water * (big_delta-(small_delta/3)) ) * 1000;
% mdd_water = 10;
fprintf('   Diffusivity \t = %g mm^2/s\n', D_water);
fprintf('   MDD water\t = %.2f micron\n', mdd_water);


%% Generate vertices on sphere
%%% Define EAP directions on sphere
if verbose
    disp('Creating pdf_dirs...');
end
if Use_DSI_Studio_dirs
    % Use DSI Studio dirs
    load('odf8.mat');
    diff_disp_sphere.vertices = odf_vertices';
    diff_disp_sphere.faces = odf_faces';
    nb_pdf_dirs_full           = length(diff_disp_sphere.vertices);
    nb_pdf_dirs                = nb_pdf_dirs_full/2;
    pdf_dirs              = diff_disp_sphere.vertices(1:nb_pdf_dirs,:);
    
    odf_verticesI = odf_vertices';
    odf_facesI = odf_faces';
else
    % Use matlab sphere function 
%     [ diff_disp_sphere ]  = create_sphere( 25 );
    [ diff_disp_sphere ]  = create_sphere_nneg( 25 );
    nb_pdf_dirs_full      = length(diff_disp_sphere.vertices);
    nb_pdf_dirs           = nb_pdf_dirs_full;
    pdf_dirs              = diff_disp_sphere.vertices(1:nb_pdf_dirs,:);

end



% f=figure('color','w','position',[900 400 900 500],'InvertHardcopy','off'); 
% ax = axes('position',[0.05 0.05 0.4 0.9]);
% scatter3(nqvec(:,1),nqvec(:,2),nqvec(:,3), ...
%     30,vecnorm(nqvec,2,2),'filled'); 
% set(ax,'Colormap',jet);
% axis equal;
% ax = axes('position',[0.5 0.05 0.4 0.9]);
% scatter3(odf_verticesI(:,1),odf_verticesI(:,2),odf_verticesI(:,3), ...
%     30,vecnorm(odf_verticesI,2,2),'filled'); 
% set(ax,'Colormap',jet); caxis([0 1])
% axis equal;



%% Generate radial displacement vector

%%% Define displacement vector
if verbose
    disp('Creating displacement vec...');
end
% Define radial points to evaluate EAP at
nr = 100; %number of points
rs = linspace(0, 1, nr)';  
dispvec  = linspace(0, mdd_water, length(rs));
r_step_size           = dispvec(2)-dispvec(1);
dispvec2 = dispvec.^2;
dispvec4 = dispvec.^4;


%% For ODFs
if ReconODFs
    if verbose
        disp('Setting odf params...');
    end
    odfpar.idx_start = 1;
    odfpar.idx_end   = 60;
    odfpar.r_weight  = 2;
    odfpar.Wmtx      = repmat(dispvec.^odfpar.r_weight, ...
        nb_pdf_dirs, 1);
    odfpar.sep       = 25; %degrees, peak separation
    odfpar.clip      = true; %clip negative eap
    
    TR2 = triangulation((odf_faces+1)',odf_vertices');
end


%% For ZDPs
if ReconZDPs
    if verbose
        disp('Setting RTPP params...');
    end
    v = [1;0;0];
    a_circle = (0:pi/50:2*pi);
    a_circle = a_circle(2:end);
    a_plane = cat(1,zeros(size(a_circle)),cos(a_circle),sin(a_circle));
    % rtap_matrix = permute(repmat(rs(2:end),1,3,100),[3 2 1]);
    % [ rtap_Fmtx ] = precompute_Fmtx( qvec, a_plane, rs );
    zpds.v = v;
    zpds.a_plane = a_plane;
    zpds.qvec = qvec;
    zpds.rs = rs;
end
        

%% Precompute Fmtx
if verbose
    disp('Creating Fmtx...');
end
[ Fmtx ] = precompute_Fmtx( qvec, pdf_dirs, rs );          

DATAIN.Fmtx = Fmtx;




%% Create output structure
if verbose
    disp('Creating output matrices...');
end
voldims = size(DATAIN.mask);

% [ DATAOUT, ~ ] = allocate_gdsi_v4( voldims, nr ); 
DATAOUT = [];

if StoreRmatrix
    DATAOUT.Rmatrix = zeros([voldims nr nb_pdf_dirs]);
end

DATAOUT.Rvec = zeros([voldims nr]);
DATAOUT.Rstd = zeros([voldims nr]);
% DATAOUT.Rrms = zeros([voldims nr]);
DATAOUT.MeanGA = zeros(voldims);
% DATAOUT.MeanEAP = zeros(voldims);
% DATAOUT.StdEAP = zeros(voldims);
% DATAOUT.RmsEAP = zeros(voldims);
DATAOUT.EntropyEAP = zeros(voldims);

DATAOUT.MeanKurt = zeros(voldims);
DATAOUT.NG = zeros(voldims);
DATAOUT.KLDg = zeros(voldims);
DATAOUT.KLDi = zeros(voldims);
DATAOUT.NLrms = zeros(voldims);

DATAOUT.RTOP = zeros(voldims);
% DATAOUT.RTAP = zeros(voldims);
% DATAOUT.RTPP = zeros(voldims);
% DATAOUT.V1 = zeros([voldims 3]);

DATAOUT.HWHM = zeros(voldims);
DATAOUT.HWHMstd = zeros(voldims);
DATAOUT.HWHMmin = zeros(voldims);
DATAOUT.HWHMmax = zeros(voldims);
DATAOUT.HWHMmin_vec = zeros([voldims 3]);
DATAOUT.HWHMmax_vec = zeros([voldims 3]);

DATAOUT.qiv = zeros(voldims);
DATAOUT.qvar = zeros(voldims);
% DATAOUT.qmsd = zeros(voldims);
% DATAOUT.qrtop = zeros(voldims);
DATAOUT.qentropy = zeros(voldims);


%% BEGIN RECON!
if verbose
    disp(' --- Starting GDSI recon --- ');
end


%%%% % % % CHANGE MASK HERE!!!! % % % %%%%
% DATAIN.mask(:,:,[1:34,36:end]) = 0;
% totalvoxels = nnz(DATAIN.mask);


%%% Voxels of interest:
% 
% CS_DSI_10_01_2021: cc --- ix=65, iy=49, iz=35
% CS_DSI_10_01_2021: cross --- ix=56, iy=45, iz=35


%%% Run the reconstruction
% 
% [ DATAOUT ] = EvalPropagator_vol( DATAIN, params, opts ) 
%
% DATAIN has:   dwi, mask, index_b0, qvol, pdf_dirs, isshelled, nr, ndirs
% DATAOUT has:  empty matrices to store measures/stats from recon
% opts has:     StoreRmatrix, SaveDATAOUT

DATAIN.index_b0 = index_b0;
DATAIN.isshelled = isshelled;
if isshelled
    DATAIN.qvol = qvol;
else
    DATAIN.qvol = [];
end
DATAIN.nr = nr;
DATAIN.ndirs = nb_pdf_dirs;
DATAIN.pdf_dirs = pdf_dirs;
DATAIN.Fmtx = Fmtx;
DATAIN.dispvec2 = dispvec2;
DATAIN.dispvec = dispvec;

opts.StoreRmatrix = false;
opts.SaveDATAOUT = false;

DATAIN.mask(DATAIN.lowb<=0)=0;

[ DATAOUT ] = EvalPropagator_vol( DATAIN, DATAOUT, opts ) ;



%%% Voxels of interest:

% foutdir = '/autofs/space/hemera_002/users/rjjones/gdap/subjects/qt';
% foutmat = [foutdir filesep 'z0037_Rmatrix.mat'];
% save(foutmat,'DATAOUT');


disp(' --- Finished reconstructing voxels --- ');
ftoc=toc(ftic);
fprintf('%.1f seconds | %d voxels\n',ftoc,totalvoxels);


disp(' - Saving output .mat file...');
% outdir = [procdir filesep 'gdsi_recon_vol'];

if ~exist(DATAIN.outdir,'dir')
    fprintf('Making outdir:\n %s\n',DATAIN.outdir);
    mkdir(DATAIN.outdir);
end

bvals               = DATAIN.bvals;
bvecs               = DATAIN.bvecs;
mask                = DATAIN.mask;
lowb                = DATAIN.lowb;

% foutmat = [outdir filesep 'gdsi_recon_csdsi_subcortical-only.mat'];
foutmat = [DATAIN.outdir filesep 'gdsi_vol_recon_prof+scal_recpar.mat']

fprintf('Saving recon to file:\n %s\n',foutmat);
save(foutmat, ... %'DATAOUT', ... 
   'pdf_dirs', 'dispvec', 'dispvec2', ...
    'diff_disp_sphere', 'qvec', 'D_water', 'mdd_water', ...
    'bvals', 'bvecs', ...
    'mask', '-v7.3');


fprintf('\n-------------\n - Finished GDSI recon - \n-------------\n');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%  HWHM min, max DEC maps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% b=squeeze(DATAOUT.HWHMmin(:,:,35));
% a=squeeze(DATAOUT.HWHMmin_vec(:,:,35,:));
% dec_hwhm_min = abs(a).*repmat((b-3),1,1,3);
% figure; imshow(dec_hwhm_min); title('DEC HWHM_[min}');
% 
% d=squeeze(DATAOUT.HWHMmax(:,:,35));
% c=squeeze(DATAOUT.HWHMmax_vec(:,:,35,:));
% dec_hwhm_max = abs(c).*repmat((d-4)/1.25,1,1,3);
% % dec_hwhm_max = abs(DATAOUT.HWHMmax_vec).*repmat((DATAOUT.HWHMmax-4)/1.25,1,1,3);
% figure; imshow(dec_hwhm_max); title('DEC HWHM_{max}','Interpreter','tex');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%  WM, GM basic stats
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% vol=load_untouch_nii('/autofs/space/nyx_002/users/rjones/Bay8/CS_DSI_10_01_2021/preproc/anat2diff/fs-tissue-masks/all-wmmask.nii.gz');
% wmmask=double(vol.img);
% wmslice = (wmmask(:,:,35));
% 
% vol=load_untouch_nii('/autofs/space/nyx_002/users/rjones/Bay8/CS_DSI_10_01_2021/preproc/anat2diff/fs-tissue-masks/ctx-gmmask.nii.gz');
% gmmask=double(vol.img);
% gmslice = (gmmask(:,:,35));
% 
% 
% curr = DATAOUT.Rvec;
% for rr=1:nr
%     tmp=curr(:,:,rr);
%     tmp(wmslice==0)=0;
%     a=mean(tmp(tmp>0));
%     b=std(tmp(tmp>0));
%     dat.rvec.wm.mean(rr)=a;
%     dat.rvec.wm.std(rr)=b;
% end
% 
% curr = DATAOUT.Rvec;
% for rr=1:nr
%     tmp=curr(:,:,rr);
%     tmp(gmslice==0)=0;
%     a=mean(tmp(tmp>0));
%     b=std(tmp(tmp>0));
%     dat.rvec.gm.mean(rr)=a;
%     dat.rvec.gm.std(rr)=b;
% end
% 
% curr = DATAOUT.Rstd;
% for rr=1:nr
%     tmp=curr(:,:,rr);
%     tmp(wmslice==0)=0;
%     a=mean(tmp(tmp>0));
%     b=std(tmp(tmp>0));
%     dat.rstd.wm.mean(rr)=a;
%     dat.rstd.wm.std(rr)=b;
% end
% 
% curr = DATAOUT.Rstd;
% for rr=1:nr
%     tmp=curr(:,:,rr);
%     tmp(gmslice==0)=0;
%     a=mean(tmp(tmp>0));
%     b=std(tmp(tmp>0));
%     dat.rstd.gm.mean(rr)=a;
%     dat.rstd.gm.std(rr)=b;
% end
% 
% figure; hold on; plot(dat.rvec.wm.mean); plot(dat.rvec.gm.mean);title('Mean EAP - dsi 10k');
% figure; hold on; plot(dat.rstd.wm.mean); plot(dat.rstd.gm.mean);title('GA profile - dsi 10k');





% % Save volumes in .nii format
fafile='/autofs/space/nyx_002/users/rjones/Bay8/CS_DSI_10_01_2021/preproc/dtifit/dtifit_FA.nii.gz';
vol = load_untouch_nii(fafile);
% vol.hdr.dime.datatype = 16;
% vol.hdr.dime.bitpix = 32;

niioutdir = [outdir filesep 'gdsi_recon_csdsi_021022_slice_z35'];
if ~exist(niioutdir,'dir')
    mkdir(niioutdir);
end




filelist = {'RTOP','StdEAP','MeanEAP','MeanGA','MeanKurt','KLDg',...
            'KLDi','NLrms','EntropyEAP','RmsEAP','NG',...
            'HWHM','HWHMstd','HWHMmin','HWHMmax',...
            'qiv','qvar','qmsd','qentropy'};
for f=1:length(filelist)
    F = filelist{f};
    eval(['curr=DATAOUT.' F ';']);
    tmp = zeros(size(vol.img));
    tmp(:,:,35)=curr;
    outfile = [niioutdir filesep  F '.nii.gz'];
    write2nifti( vol, tmp, outfile );
end

filelist = {'HWHMmin_vec','HWHMmax_vec'  }; 
for f=1:length(filelist)
    F = filelist{f};
    eval(['curr=DATAOUT.' F ';']);
    tmp = zeros([size(vol.img) 3]);
    tmp(:,:,35,:)=curr;
    outfile = [niioutdir filesep  F '.nii.gz'];
    write2nifti( vol, tmp, outfile );
end


filelist = {'Rvec','Rstd'};
for f=1:length(filelist)
    F = filelist{f};
    eval(['curr=DATAOUT.' F ';']);
    tmp = zeros([size(vol.img) 100]);
    tmp(:,:,35,:) = curr;
    outfile = [niioutdir filesep  F '.nii.gz'];
    write2nifti( vol, tmp, outfile );
end



% % Save volumes in .nii format
niidir = [outdir filesep 'nii'];
if ~exist(niidir,'dir'), mkdir(niidir); end

% fafile='/autofs/space/nyx_002/users/rjones/Bay8/CS_DSI_10_01_2021/preproc/dtifit/dtifit_FA.nii.gz';
% fafile = [filedir filesep 'dtifit' filesep 'dtifit_FA.nii.gz'];
fafile = '/autofs/space/nyx_002/users/rjones/Bay8/CS_DSI_10_01_2021/diff/preproc/mri/eddy/results_def/dtifit/dtifit_FA.nii.gz';

vol = load_untouch_nii(fafile);
% vol.hdr.dime.datatype = 16;
% vol.hdr.dime.bitpix = 32;


filelist = {'RTOP', 'MeanGA', 'MeanKurt', ...
            'KLDg', 'KLDi', 'NG', 'NLrms', ...
            'HWHM', 'HWHMstd', 'HWHMmin', 'HWHMmax', ...
            'qiv', 'qentropy', 'qvar'};
for f=1:length(filelist)
    F = filelist{f};
    eval(['curr=DATAOUT.' F ';']);
    outfile = [niidir filesep F '.nii.gz'];
    write2nifti( vol, curr, outfile );
end

filelist = {'Rvec','Rstd','HWHMmin_vec','HWHMmax_vec'};
for f=1:length(filelist)
    F = filelist{f};
    eval(['curr=DATAOUT.' F ';']);
    outfile = [niidir filesep F '.nii.gz'];
    write2nifti( vol, curr, outfile );
end

