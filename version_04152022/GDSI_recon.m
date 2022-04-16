% function [ DATAOUT ] = GDSI_recon( fconfig, verbose )
% 
% GDSI reconstruction of RSI data
%
% original script by (c) Qiyuan Tian, Stanford RSL
% Modified by RJ, Feb 2021 - Feb 2022
%

%%%%%%%%%%%%%%
% Things to try:

% --- EXAMPLE ---
fconfig = '/autofs/space/nyx_002/users/rjones/Bay8/code/gdsi/version_020122/hcp_mgh.config';
verbose = true;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Hard-coded variables ---
Dwater = 2.5e-3;
Use_DSI_Studio_dirs     = true;

Remove_HighBValue       = false;
ClipRmatrix             = false;
Norm2PDF                = false;
ReconODFs               = false;
ReconZDPs               = false;
ReconMDs                = false;
MakePlots.odf           = false;
MakePlots.rgbprofiles   = false;
StoreRmatrix            = false;
concatlowb              = false;


% %%% For later, to clean argins
% Use_DSI_Studio_dirs = opts.Use_DSI_Studio_dirs;
% Remove_HighBValue = opts.Remove_HighBValue;
% ClipRmatrix = opts.ClipRmatrix;
% Norm2PDF = opts.Norm2PDF;
% ReconODFs = opts.ReconODFs;
% ReconZDPs = opts.ReconZDPs;
% ReconMDs  = opts.ReconMDs;
% MakePlots.odf = opts.MakePlots.odf;
% MakePlots.rgbprofiles = opts.MakePlots.rgbprofiles;
% StoreRmatrix            = opts.StoreRmatrix;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isdeployed
    disp('--RUNNING DEPLOYED APP--');
else
    addpath('/autofs/space/hemera_001/users/rjjones/code/toolboxes/NIfTI_20140122');
    addpath(genpath('/autofs/space/nyx_002/users/rjones/Bay8/code/gdsi/version_020122'));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if nargin<2, verbose=0; end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Load subject information
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if verbose
    disp('Loading config file...');
end
[ INFO ] = read_cfg_file( fconfig );

INFO.flowb = [INFO.procdir filesep 'lowb.nii.gz'];
if ~exist(INFO.flowb,'file'), INFO.flowb=[]; end

INFO.ODFoutdir = [INFO.outdir filesep 'gdsi_odf'];
if ~exist(INFO.ODFoutdir,'dir'), mkdir(INFO.ODFoutdir); end

INFO.fvecs = [INFO.procdir filesep 'dtifit/dtifit_V1.nii.gz'];
if ~exist(INFO.fvecs,'file'), INFO.fvecs=[]; end



a=load_untouch_nii('/autofs/space/nyx_002/users/rjones/Bay8/CS_DSI_10_01_2021/preproc/brain_mask.nii.gz');
fullmask=double(a.img);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% TO ONLY DO SUBCORTICAL!!!!!!

% INFO.fmask = '/autofs/space/nyx_002/users/rjones/Bay8/CS_DSI_10_01_2021/preproc/subcortical_brain_mask.nii.gz';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% TO ONLY DO SUBCORTICAL!!!!!!

INFO.fmask = '/autofs/space/nyx_002/users/rjones/Bay8/CS_DSI_10_01_2021/preproc/non_subcortical_brain_mask.nii.gz';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%







if verbose
    disp('Loading acq specs file...');
end
[ ACQ ] = read_acqspecs_file( INFO.fspecs );


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Load diffusion data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if verbose
    disp('Loading subject data...');
end
[ DATAIN ] = load_subject_data( INFO, verbose ); %opts.loaddwi=true
%%%%

% DATAIN = [];
% 
% %%%% DWI
% vol = load_untouch_nii(INFO.fdwi);
% DATAIN.dwi = double(vol.img);
% DATAIN.lowb = DATAIN.dwi(:,:,:,1);
% 
% %%%% BRAIN MASK
% if verbose
%     disp('Loading mask...');
% end
% vol = load_untouch_nii(INFO.fmask);
% DATAIN.mask = double(vol.img);
% 
% 
% %%%% BVEC/BVAL
% if verbose
%     disp('Loading bvecs/bvals...');
% end
% DATAIN.bvecs = dlmread(INFO.fbvec);
% DATAIN.bvals = dlmread(INFO.fbval);



DATAIN.mask(fullmask==0)=0;



% Check that dims of bvecs/bvals/dwi match up
validate_diff_data( DATAIN );
if verbose
    disp(' -- Diffusion data validated --');
end

% Find b=0 indices
index_b0 = DATAIN.bvals==0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Prepare for GDSI recon
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define diffusivity of free water (in vivo)
% Dwater = 2.5e-3;    % mm^2*s^-1


%%%%%%%%%%%%%
%%%%% NOTE: 
%%%%%   Orginally not sure what Delta/delta were, so used values from Fan
%%%%%   et al. (i think, or related bay 8 paper). 
%%%%%    - was using Delta/delta = 21.8/12.9 ms --> MDDwater = 16.20 um
%%%%%    - actually  Delta/delta = 25.5/10.8 ms --> MDDwater = 18.12 um    
%%%%%

%%% Calculate mean diffusion distance of water (micron)
mdd_water = sqrt(6 * Dwater * (ACQ.big_delta-(ACQ.small_delta/3)) ) * 1000;
% mdd_water = 10;
fprintf(' - Diffusivity water = %g mm^2/s\n', Dwater);
fprintf(' - MDD water = %.2f micron\n', mdd_water);

%%% Create qvec
if verbose
    disp('Creating qvec...');
end
if size(DATAIN.bvals,2)>size(DATAIN.bvals,1)
    DATAIN.bvals = DATAIN.bvals';
end
if size(DATAIN.bvecs,2)>size(DATAIN.bvecs,1)
    DATAIN.bvecs = DATAIN.bvecs';
end
if concatlowb
    DATAIN.bvals_orig = DATAIN.bvals;
    DATAIN.bvecs_orig = DATAIN.bvecs;
    DATAIN.bvals = DATAIN.bvals(~index_b0);
    DATAIN.bvals = cat(1,0,DATAIN.bvals);
    DATAIN.bvecs = DATAIN.bvecs(~index_b0,:);
    DATAIN.bvecs = cat(1,[0 0 0],DATAIN.bvecs);
end
qvec = repmat(sqrt( 6 * Dwater * DATAIN.bvals ), [1, 3]) .* DATAIN.bvecs;


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
    
else
    % Use matlab sphere function 
%     [ diff_disp_sphere ]  = create_sphere( 25 );
    [ diff_disp_sphere ]  = create_sphere_nneg( 25 );
    nb_pdf_dirs_full      = length(diff_disp_sphere.vertices);
    nb_pdf_dirs           = nb_pdf_dirs_full;
    pdf_dirs              = diff_disp_sphere.vertices(1:nb_pdf_dirs,:);

end

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

%%% For ODFs
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

%%% For RTPP
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
        

%%% Precompute Fmtx
if verbose
    disp('Creating Fmtx...');
end
[ Fmtx ] = precompute_Fmtx( qvec, pdf_dirs, rs );          


%%% Create output structure
if verbose
    disp('Creating output matrices...');
end
voldims = size(DATAIN.mask);
[ DATAOUT, ~ ] = allocate_gdsi_v4( voldims, nr ); 
if StoreRmatrix
    DATAOUT.Rmatrix = zeros([voldims nr nb_pdf_dirs]);
end
DATAOUT.HWHMmin = zeros(voldims);
DATAOUT.HWHMmax = zeros(voldims);
DATAOUT.HWHMmin_vec = zeros([voldims 3]);
DATAOUT.HWHMmax_vec = zeros([voldims 3]);

%%% BEGIN RECON!
if verbose
    disp(' --- Starting GDSI recon --- ');
end

%%%% % % % CHANGE MASK HERE!!!! % % % %%%%
% DATAIN.mask(:,:,[1:24,26:end]) = 0;

totalvoxels = nnz(DATAIN.mask);

%%% Voxels of interest:
% 
% CS_DSI_10_01_2021: cc --- ix=65, iy=49, iz=35
% CS_DSI_10_01_2021: cross --- ix=56, iy=45, iz=35

cnt = 0;
ftic = tic;
rectimer = tic; timerstepsize = 100;
for ix=1:size(DATAIN.mask,1)
    for iy=1:size(DATAIN.mask,2)
        for iz=1:size(DATAIN.mask,3)
            if DATAIN.mask(ix,iy,iz)>0
                
                if verbose && cnt>0 && mod(cnt,timerstepsize)==0
                    currtimer = toc(rectimer);
                    fprintf(' - %d/%d vox = %.02f pct\n',cnt,totalvoxels,(cnt/totalvoxels)*100);
                    fprintf('  - Elapsed time:             %.01f sec \n',currtimer);
                    secremaining = round((currtimer/cnt)*((totalvoxels-cnt)));
                    fprintf('  - Estimated time remaining: %d h %d m %d s \n',...
                        floor(secremaining/3600),...
                        floor(mod(secremaining,3600)/60),...
                        mod(secremaining,60));
                end

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % RECONSTRUCT EAP
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                % Get voxel qspace data as col vector
                data = reshape(DATAIN.dwi(ix,iy,iz,:),[],1);
                b0s = data(index_b0);       % b0 volumes
                dwis = data(~index_b0);     % dwi volumes
                data = [mean(b0s); dwis];   % average b0 volumes
                data = data / mean(b0s);    % normalize diff signal
                
                %%% Compute 1D EAP along each pdf_dir
                [ Rmatrix ] = compute_1d_eap_mtx( data, Fmtx );
                Rvec = mean(Rmatrix,2);
                Rstd = std(Rmatrix,0,2);
                
                DATAOUT.MeanEAP(ix,iy,iz)        = mean(Rmatrix(:));  %mean of eap
                Rsigma     = std(Rmatrix(:));   %std of eap
                DATAOUT.StdEAP(ix,iy,iz)        = Rsigma;
                DATAOUT.RTOP(ix,iy,iz)          = Rmatrix(1,1);      %rtop of eap
                DATAOUT.Rvec(ix,iy,iz,:)          = Rvec;      
                DATAOUT.Rstd(ix,iy,iz,:)          = Rstd;      
                DATAOUT.MeanGA(ix,iy,iz)          = mean(Rstd);    
                
                %%%%%%%% (NON)GAUSSIAN METRICS  %%%%%%%%
                % Divergence of Mean 1D EAP from Gaussian
                gaussmtx = exp(-dispvec2'*Rsigma);        %make gaussian
                DATAOUT.GaussDiv(ix,iy,iz)   = abs(KLDiv(Rvec',gaussmtx'));   %compute KL divergence
                
                % Calculate non-Gaussianity (NG) [Ning 2015]
                costheta = dot(Rvec,gaussmtx)/sqrt(dot(Rvec,Rvec)*dot(gaussmtx,gaussmtx));
                DATAOUT.sintheta(ix,iy,iz) = sqrt(1-costheta^2);
                
                %%%%%%%% KURTOSIS %%%%%%%%
                % Calculate kurtosis at each distance r
                DATAOUT.KurtProfile(ix,iy,iz,:) = kurtosis(Rmatrix');
                
                % Calculate average eap kurtosis
                DATAOUT.MeanKurt(ix,iy,iz) = kurtosis(Rmatrix(:));
                
                % Mean squared displacement
                DATAOUT.MSD(ix,iy,iz) = sum(dispvec2' .* Rvec);
                DATAOUT.MFD(ix,iy,iz) = sum(dispvec4' .* Rvec);
                
                % HWHM radius (only half width, not full!)
                [ hwhm ] = compute_ralpha( Rmatrix, 0.5, dispvec );
                DATAOUT.HWHM(ix,iy,iz) = hwhm.mean;
                DATAOUT.HWHMstd(ix,iy,iz) = hwhm.std;
                DATAOUT.HWHMmin(ix,iy,iz) = hwhm.min;
                DATAOUT.HWHMmax(ix,iy,iz) = hwhm.max;
                
                tmp_mvec = pdf_dirs((hwhm.vec==hwhm.min),:);
                if size(tmp_mvec,1)>1
                    tmp_mvec = tmp_mvec(1,:);
                end
                DATAOUT.HWHMmin_vec(ix,iy,iz,:) = tmp_mvec;
                
                tmp_mvec = pdf_dirs((hwhm.vec==hwhm.max),:);
                if size(tmp_mvec,1)>1
                    tmp_mvec = tmp_mvec(1,:);
                end
                DATAOUT.HWHMmax_vec(ix,iy,iz,:) = tmp_mvec;
                
                if StoreRmatrix
                    DATAOUT.Rmatrix(ix,iy,iz,:,:) = Rmatrix;
                end
                
                % Normalize EAP to a PDF [i.e., sum(eap(:))=1]
                if Norm2PDF
                    normRmatrix = Rmatrix/sum(Rmatrix(:));
                    DATAOUT.GAProfile_norm(ix,iy,iz,:) = std(normRmatrix,0,2);
                    DATAOUT.RTOP_norm(ix,iy,iz) = normRmatrix(1,1);
                end
                
                
                if ClipRmatrix
                    % Clip after first zero crossing
                    [ nRmatrix ] = clip_negative_Rmatrix( Rmatrix );
%                     [ rRmatrix ] = clip_ringing_Rmatrix( Rmatrix );
%                     [ nrRmatrix ] = clip_negative_Rmatrix( rRmatrix );
                    nRvec = mean(nRmatrix, 2);
                    nRstd = std( nRmatrix, 0, 2);

                    if MakePlots.rgbprofiles
                        [ f1 ] = plot_all_profiles_voxel( dispvec, Rmatrix, ...
                            abs(diff_disp_sphere.vertices)' );
                         hold on; plot(dispvec,zeros(size(dispvec)),'-k');
                         
                        pause(1);
                        
                        [ f1 ] = plot_all_profiles_voxel( dispvec, nRmatrix, ...
                            abs(diff_disp_sphere.vertices)' ); 
                         hold on; plot(dispvec,zeros(size(dispvec)),'-k');
                    end
                end
                

                
                if ReconODFs
                    %%%%%%%% ODF %%%%%%%%
                    % Recon dODF
                    Omatrix = Rmatrix.*odfpar.Wmtx';
                    if Use_DSI_Studio_dirs
                        Omatrix = repmat(Omatrix,1,2);
                    end                
                    ODF = sum( Omatrix( odfpar.idx_start:odfpar.idx_end, :) );
                    ODF_points = diff_disp_sphere.vertices .* repmat(ODF', [1, 3]);
                    
                    [ f1 ] = plot_all_profiles_voxel( dispvec, Omatrix, ...
                            abs(diff_disp_sphere.vertices)' ); 
                         hold on; plot(dispvec,zeros(size(dispvec)),'-k');
                         
                         
                         
                    nOmatrix = nRmatrix.*odfpar.Wmtx';
                    if Use_DSI_Studio_dirs
                        nOmatrix = repmat(nOmatrix,1,2);
                    end                
                    nODF = sum( nOmatrix ); %( odfpar.idx_start:odfpar.idx_end, :) );
                    nODF_points = diff_disp_sphere.vertices .* repmat(nODF', [1, 3]);
                    
                    [ f1 ] = plot_all_profiles_voxel( dispvec, nOmatrix, ...
                            abs(diff_disp_sphere.vertices)' ); 
                         hold on; plot(dispvec,zeros(size(dispvec)),'-k');
                         
                    [pks2, allpks2] = gdsi_find_peak(nODF', diff_disp_sphere, odfpar.sep);     
                         

                    % find ODF peaks
                    [pks, allpks] = gdsi_find_peak(ODF', diff_disp_sphere, odfpar.sep);
                    npeaks.D19 = length(pks.D19.pkvals);
                    totalpeaks.D19 = npeaks.D19;
                    if npeaks.D19>3                    
                        npeaks.D19=3;
                    end
                    
                    %%% Compute RTAP, RTPP on major fiber axis
                    pdd = pks.pkvecs(1,:);
                    [ rtap ] = compute_RTAP( data, v, pdd, qvec, rs, dispvec2 );
                    [ rtpp ] = compute_RTPP( data, v, a_plane, pdd, qvec, rs, dispvec2 );

                    %%% Compute fixel metrics
                    [fixels] = compute_fixel_metrics_v2( pks, diff_disp_sphere.vertices, ...
                        Rmatrix, data, dispvec, zpds); 


                    if MakePlots.odf
                        %%%%% plot ODF
    %                     plotODF = ODF;
                        plotODF = ODF/max(ODF);
        % %                 plotODF = (ODF-min(ODF))/(max(ODF)-min(ODF));
                        plotPeaks = pks.pkvecs;

                        [ f ] = gdsi_plot_odf( diff_disp_sphere, plotODF, plotPeaks );
                        
                        
                        
                        
                        plotODF = nODF/max(nODF);
        % %                 plotODF = (ODF-min(ODF))/(max(ODF)-min(ODF));
                        plotPeaks = pks2.pkvecs;

                        [ f ] = gdsi_plot_odf( diff_disp_sphere, plotODF, plotPeaks );
                        
                        
                    end
                
                    %%% STORE RESULTS
                    % ODF GFA, RMS, STD, npeaks, iso, mass, max
                    DATAODF.GFA(ix,iy,iz) = std(nODFhalf)/rms(nODFhalf);
                    DATAODF.RMS(ix,iy,iz) = rms(nODFhalf);
                    DATAODF.STD(ix,iy,iz) = std(nODFhalf);
                    DATAODF.npeaks(ix,iy,iz) = n_totalpeaks;
                    % ODF iso (min value)
                    DATAODF.iso(ix,iy,iz) = min(nODF);
                    % ODF mass (sum of all odf values)
                    DATAODF.mass(ix,iy,iz) = sum(nODF);

                    %%% Fixels
                    for nn=1:n_npeaks
                        eval(['curr = n_fixels.v' num2str(nn) ';']);
                        currind = curr.diridx;
                        currO = nOmatrix(odfpar.idx_start:odfpar.idx_end, currind);
                        currfa = sum(currO);
    %                     currnfa = currfa/sum(currO);

                        eval(['DATAODF.ind' num2str(nn) '(ix,iy,iz) = curr.diridx;']);
                        eval(['DATAODF.vec' num2str(nn) '(ix,iy,iz,:) = curr.dir;']);

                        eval(['DATAODF.f' num2str(nn) '(ix,iy,iz) = curr.f;']);
                        eval(['DATAODF.fa' num2str(nn) '(ix,iy,iz) = curr.fa;']);

                        eval(['DATAODF.Rvec' num2str(nn) '(ix,iy,iz,:) = curr.Rvec;']);
                        eval(['DATAODF.rtap' num2str(nn) '(ix,iy,iz) = curr.Rmean;']);
                        eval(['DATAODF.Rstd' num2str(nn) '(ix,iy,iz) = curr.Rstd;']);
                        eval(['DATAODF.Rmsd' num2str(nn) '(ix,iy,iz) = curr.Rmsd;']);            
                        eval(['DATAODF.kurt' num2str(nn) '(ix,iy,iz) = curr.Kurtosis;']);                    
                        eval(['DATAODF.ip' num2str(nn) '(ix,iy,iz) = curr.inflpnt;']);
                        eval(['DATAODF.r_50_' num2str(nn) '(ix,iy,iz) = curr.r_50;']);
                        eval(['DATAODF.r_25_' num2str(nn) '(ix,iy,iz) = curr.r_25;']);

                        eval(['DATAODF.rtpp.Rvec' num2str(nn) '(ix,iy,iz,:) = curr.rtpp.Rvec;']);                    
                        eval(['DATAODF.rtpp.Rstd' num2str(nn) '(ix,iy,iz,:) = curr.rtpp.Rstd;']);                    
                        eval(['DATAODF.rtpp.rtpp' num2str(nn) '(ix,iy,iz,:) = curr.rtpp.RTPP;']);                    
                        eval(['DATAODF.rtpp.rtpv' num2str(nn) '(ix,iy,iz,:) = curr.rtpp.RTPV;']);                    
                        eval(['DATAODF.rtpp.msd' num2str(nn) '(ix,iy,iz,:) = curr.rtpp.msd;']);                    
                        eval(['DATAODF.rtpp.kurt' num2str(nn) '(ix,iy,iz,:) = curr.rtpp.kurtosis;']);                    

                    end
                end
                
                if ReconZDPs
                    % ZDPs
                    DATAOUT.RTOP(ix,iy,iz) = Rvec(1);
                    DATAOUT.RTAP(ix,iy,iz) = n_fixels.v1.Rmean;
                    DATAOUT.RTPP(ix,iy,iz) = n_fixels.v1.rtpp.RTPP;
                    DATAOUT.RTAV(ix,iy,iz) = rtap.RTAV;
                    DATAOUT.RTAP_vec(ix,iy,iz,:) = rtap.Rvec;
                    DATAOUT.RTAV_vec(ix,iy,iz,:) = rtap.Rstd;
                end
                if ReconMDs
                    % MD, MSD, MFD - raw EAP 
                    DATAOUT.MD1(ix,iy,iz) = Rmd1;
                    DATAOUT.MD2(ix,iy,iz) = Rmd2;
                    DATAOUT.MD4(ix,iy,iz) = Rmd4;
                end
                
                cnt=cnt+1;
            end
        end
    end
end

% foutdir = '/autofs/space/hemera_002/users/rjjones/gdap/subjects/qt';
% foutmat = [foutdir filesep 'z0037_Rmatrix.mat'];
% save(foutmat,'DATAOUT');


disp(' --- Finished reconstructing voxels --- ');
ftoc=toc(ftic);
fprintf('%.1f seconds | %d voxels\n',ftoc,totalvoxels);


disp(' - Saving output .mat file...');
% INFO.outdir = [INFO.procdir filesep 'gdsi_recon_vol'];

if ~exist(INFO.outdir,'dir')
    fprintf('Making outdir:\n %s\n',INFO.outdir);
    mkdir(INFO.outdir);
end

bvals               = DATAIN.bvals;
bvecs               = DATAIN.bvecs;
nodif_brain_mask    = DATAIN.mask;
lowb                = DATAIN.lowb;

% foutmat = [INFO.outdir filesep 'gdsi_recon_csdsi_subcortical-only.mat'];
foutmat = [INFO.outdir filesep 'gdsi_recon_csdsi_not-subcortical-only.mat'];

fprintf('Saving recon to file:\n %s\n',foutmat);
save(foutmat, 'DATAOUT', 'INFO', 'ACQ',...
   'pdf_dirs', 'dispvec', ...
    'diff_disp_sphere', 'qvec', 'Dwater', 'mdd_water', ...
    'bvals', 'bvecs', 'totalvoxels', 'ftoc', ...
    'nodif_brain_mask', '-v7.3');


fprintf('\n-------------\n - Finished GDSI recon - \n-------------\n');


if 0 == 1
    
    segmasks = [];
    
    maskdir = [INFO.procdir filesep 'segm'];
    
    fwmmask = [maskdir filesep 'mask.all-wm.nii.gz'];
    vol = load_untouch_nii(fwmmask);
    segmasks.wm = double(vol.img);
    
    fgmmask = [maskdir filesep 'mask.gm.nii.gz'];
    vol = load_untouch_nii(fgmmask);
    segmasks.gm = double(vol.img);
    
    fsubcortmask = [maskdir filesep 'mask.subcort-gm.nii.gz'];
    vol = load_untouch_nii(fsubcortmask);
    segmasks.subcort = double(vol.img);
    
    fcsfmask = [maskdir filesep 'mask.ventricles.nii.gz'];
    vol = load_untouch_nii(fcsfmask);
    segmasks.csf = double(vol.img);
    
    
    
    slices = [];

    filelist = {'RTOP','MeanEAP','MeanGA','MSD',...
                'HWHM','HWHMstd','HWHMminmax', ...
                'RmatrixSum0','RTOP0','MeanEAP0','StdEAP','StdEAP0',...
                'GaussDiv','NG','sintheta','MD','MeanKurt'};
    for f=1:length(filelist)
        F = filelist{f};
        eval(['curr=DATAOUT.' F ';']);
        currslice = squeeze(curr(:,:,92));
        eval(['slices.' F ' = currslice;']);
        
    end
    
    
    mask = DATAIN.mask;
    lowb = DATAIN.lowb;
    wmmask = segmasks.wm;
    gmmask = segmasks.gm;
    csfmask = segmasks.csf;
    subcortmask = segmasks.subcort;
    
    filelist = {'mask','lowb','wmmask','gmmask','csfmask','subcortmask'};
    for f=1:length(filelist)
        F = filelist{f};
        eval(['curr=' F ';']);
        currslice = squeeze(curr(:,:,92));
        eval(['slices.' F ' = currslice;']);
        
    end
    
    
    save(foutmat,'slices','segmasks','-append');
    
    
    
    
    wmRvec = zeros(1,100);
    gmRvec = zeros(1,100);
    csfRvec = zeros(1,100);
    scRvec = zeros(1,100);
    
    wmRvec0 = zeros(1,100);
    gmRvec0 = zeros(1,100);
    csfRvec0 = zeros(1,100);
    scRvec0 = zeros(1,100);
    
    wmRstd = zeros(1,100);
    gmRstd = zeros(1,100);
    csfRstd = zeros(1,100);
    scRstd = zeros(1,100);
    wmKurt = zeros(1,100);
    gmKurt = zeros(1,100);
    csfKurt = zeros(1,100);
    scKurt = zeros(1,100);
    
    for rr=1:100
        curr = DATAOUT.Rvec(:,:,92,rr);
        wmRvec(rr) = mean(curr(slices.wmmask==1));
        gmRvec(rr) = mean(curr(slices.gmmask==1));
        csfRvec(rr) = mean(curr(slices.csfmask==1));
        scRvec(rr) = mean(curr(slices.subcortmask==1));
        curr = DATAOUT.Rstd(:,:,92,rr);
        wmRstd(rr) = mean(curr(slices.wmmask==1));
        gmRstd(rr) = mean(curr(slices.gmmask==1));
        csfRstd(rr) = mean(curr(slices.csfmask==1));
        scRstd(rr) = mean(curr(slices.subcortmask==1));
        curr = DATAOUT.KurtProfile(:,:,92,rr);
        wmKurt(rr) = mean(curr(slices.wmmask==1));
        gmKurt(rr) = mean(curr(slices.gmmask==1));
        csfKurt(rr) = mean(curr(slices.csfmask==1));
        scKurt(rr) = mean(curr(slices.subcortmask==1));
    end
    
    
    
    
    
%     temp1 = DATAOUT.RmatrixSum0;
%     temp2 = DATAOUT.Rvec;
%     temp3 = temp2 .* repmat(temp1,1,1,1,100);
%     
%     for rr=1:100
%         curr = temp3(:,:,92,rr);
%         wmRvec0(rr) = mean(curr(slices.wmmask==1));
%         gmRvec0(rr) = mean(curr(slices.gmmask==1));
%         csfRvec0(rr) = mean(curr(slices.csfmask==1));
%         scRvec0(rr) = mean(curr(slices.subcortmask==1));
%     end
%         
%         
%     Rvec0 = temp3;
%     save(foutmat,'slices','segmasks','Rvec0','-append');
    
    
    INFO.plotdir = [INFO.outdir filesep 'plots'];
    if ~exist(INFO.plotdir,'dir'), mkdir(INFO.plotdir); end
    
    if 1 == 1
        f=figure; f.Color='w'; f.InvertHardcopy='off'; f.Position=[670 659 492 441];
        hold on;   
        plot(dispvec,csfRvec,'-g'); plot(dispvec,scRvec,'-c');
        plot(dispvec,wmRvec,'-r'); plot(dispvec,gmRvec,'-b');
        aaa=findall(gcf,'type','line'); set(aaa,'LineWidth',2);
        plot(dispvec,zeros(size(dispvec)),'-k');
        legend('Ventr','Subcort','WM','GM','location','best');
        title('Rvec'); xlabel('Disp.'); ylabel('Prob. dens.');
        fout = [INFO.plotdir filesep 'no-norm_Rvec_tissues.jpg'];
        print(f,fout,'-djpeg');

    %     f=figure; f.Color='w'; f.InvertHardcopy='off'; f.Position=[670 659 492 441];
    %     hold on;  
    %     plot(dispvec,csfRvec0,'-g'); plot(dispvec,scRvec0,'-c');
    %     plot(dispvec,wmRvec0,'-r'); plot(dispvec,gmRvec0,'-b'); 
    %     aaa=findall(gcf,'type','line'); set(aaa,'LineWidth',2);
    %     plot(dispvec,zeros(size(dispvec)),'-k');
    %     legend('Ventr','Subcort','WM','GM','location','best');
    %     title('Rvec0'); xlabel('Disp.'); ylabel('Prob. dens.');
    %     fout = [INFO.plotdir filesep 'norm_Rvec_tissues.jpg'];
    %     print(f,fout,'-djpeg');
    %     
        f=figure; f.Color='w'; f.InvertHardcopy='off'; f.Position=[670 659 492 441];
        hold on; 
        plot(dispvec,csfRstd,'-g'); plot(dispvec,scRstd,'-c');
        plot(dispvec,wmRstd,'-r'); plot(dispvec,gmRstd,'-b'); 
        aaa=findall(gcf,'type','line'); set(aaa,'LineWidth',2);
        plot(dispvec,zeros(size(dispvec)),'-k');
        legend('Ventr','Subcort','WM','GM','location','best');
        title('Rstd'); xlabel('Disp.'); ylabel('Gen. aniso.');
        fout = [INFO.plotdir filesep 'no-norm_Rstd_tissues.jpg'];
        print(f,fout,'-djpeg');

        f=figure; f.Color='w'; f.InvertHardcopy='off'; f.Position=[670 659 492 441];
        hold on; 
        plot(dispvec,csfKurt,'-g'); plot(dispvec,scKurt,'-c');
        plot(dispvec,wmKurt,'-r'); plot(dispvec,gmKurt,'-b'); 
        aaa=findall(gcf,'type','line'); set(aaa,'LineWidth',2);
        plot(dispvec,zeros(size(dispvec)),'-k');
        yyl=ylim;
        ylim([2 yyl(2)]);
        legend('Ventr','Subcort','WM','GM','location','best');
        title('Kurt'); xlabel('Disp.'); ylabel('Kurtosis');
        fout = [INFO.plotdir filesep 'no-norm_Kurt_tissues.jpg'];
        print(f,fout,'-djpeg');
    end
    
    if 0 == 1
        
        normresults = load([INFO.outdir filesep 'gdsi_recon_normPDF.mat'],'DATAOUT','slices','Rvec0');
        
        wmRvec2 = zeros(1,100);
        gmRvec2 = zeros(1,100);
        csfRvec2 = zeros(1,100);
        scRvec2 = zeros(1,100);

        wmRvec02 = zeros(1,100);
        gmRvec02 = zeros(1,100);
        csfRvec02 = zeros(1,100);
        scRvec02 = zeros(1,100);

        wmRstd2 = zeros(1,100);
        gmRstd2 = zeros(1,100);
        csfRstd2 = zeros(1,100);
        scRstd2 = zeros(1,100);
        wmKurt2 = zeros(1,100);
        gmKurt2 = zeros(1,100);
        csfKurt2 = zeros(1,100);
        scKurt2 = zeros(1,100);

        for rr=1:100
            curr = normresults.DATAOUT.Rvec(:,:,92,rr);
            wmRvec2(rr) = mean(curr(normresults.slices.wmmask==1));
            gmRvec2(rr) = mean(curr(normresults.slices.gmmask==1));
            csfRvec2(rr) = mean(curr(normresults.slices.csfmask==1));
            scRvec2(rr) = mean(curr(normresults.slices.subcortmask==1));
            curr = normresults.DATAOUT.Rstd(:,:,92,rr);
            wmRstd2(rr) = mean(curr(normresults.slices.wmmask==1));
            gmRstd2(rr) = mean(curr(normresults.slices.gmmask==1));
            csfRstd2(rr) = mean(curr(normresults.slices.csfmask==1));
            scRstd2(rr) = mean(curr(normresults.slices.subcortmask==1));
            curr = normresults.DATAOUT.KurtProfile(:,:,92,rr);
            wmKurt2(rr) = mean(curr(normresults.slices.wmmask==1));
            gmKurt2(rr) = mean(curr(normresults.slices.gmmask==1));
            csfKurt2(rr) = mean(curr(normresults.slices.csfmask==1));
            scKurt2(rr) = mean(curr(normresults.slices.subcortmask==1));
            curr = normresults.Rvec0(:,:,92,rr);
            wmRvec02(rr) = mean(curr(normresults.slices.wmmask==1));
            gmRvec02(rr) = mean(curr(normresults.slices.gmmask==1));
            csfRvec02(rr) = mean(curr(normresults.slices.csfmask==1));
            scRvec02(rr) = mean(curr(normresults.slices.subcortmask==1));
        end




        %%%%%%%%%%%%%5
        f=figure; f.Color='w'; f.InvertHardcopy='off'; f.Position=[670 659 492 441];
        hold on;   
        plot(dispvec,csfRvec2,'-g'); plot(dispvec,scRvec2,'-c');
        plot(dispvec,wmRvec2,'-r'); plot(dispvec,gmRvec2,'-b');
        aaa=findall(gcf,'type','line'); set(aaa,'LineWidth',2);
        plot(dispvec,zeros(size(dispvec)),'-k');
        legend('Ventr','Subcort','WM','GM','location','best');
        title('Rvec'); xlabel('Disp.'); ylabel('Prob. dens.');
        fout = [INFO.plotdir filesep 'norm_Rvec_tissues.jpg'];
        print(f,fout,'-djpeg');

        f=figure; f.Color='w'; f.InvertHardcopy='off'; f.Position=[670 659 492 441];
        hold on;  
        plot(dispvec,csfRvec02,'-g'); plot(dispvec,scRvec02,'-c');
        plot(dispvec,wmRvec02,'-r'); plot(dispvec,gmRvec02,'-b'); 
        aaa=findall(gcf,'type','line'); set(aaa,'LineWidth',2);
        plot(dispvec,zeros(size(dispvec)),'-k');
        legend('Ventr','Subcort','WM','GM','location','best');
        title('Rvec0'); xlabel('Disp.'); ylabel('Prob. dens.');
        fout = [INFO.plotdir filesep 'norm_Rvec0_tissues.jpg'];
        print(f,fout,'-djpeg');
    %     
        f=figure; f.Color='w'; f.InvertHardcopy='off'; f.Position=[670 659 492 441];
        hold on; 
        plot(dispvec,csfRstd2,'-g'); plot(dispvec,scRstd2,'-c');
        plot(dispvec,wmRstd2,'-r'); plot(dispvec,gmRstd2,'-b'); 
        aaa=findall(gcf,'type','line'); set(aaa,'LineWidth',2);
        plot(dispvec,zeros(size(dispvec)),'-k');
        legend('Ventr','Subcort','WM','GM','location','best');
        title('Rstd'); xlabel('Disp.'); ylabel('Gen. aniso.');
        fout = [INFO.plotdir filesep 'norm_Rstd_tissues.jpg'];
        print(f,fout,'-djpeg');

        f=figure; f.Color='w'; f.InvertHardcopy='off'; f.Position=[670 659 492 441];
        hold on; 
        plot(dispvec,csfKurt,'-g'); plot(dispvec,scKurt,'-c');
        plot(dispvec,wmKurt,'-r'); plot(dispvec,gmKurt,'-b'); 
        aaa=findall(gcf,'type','line'); set(aaa,'LineWidth',2);
        plot(dispvec,zeros(size(dispvec)),'-k');
        yyl=ylim;
        ylim([2 yyl(2)]);
        legend('Ventr','Subcort','WM','GM','location','best');
        title('Kurt'); xlabel('Disp.'); ylabel('Kurtosis');
        fout = [INFO.plotdir filesep 'norm_Kurt_tissues.jpg'];
        print(f,fout,'-djpeg'); 
    end
end



if 0 == 1
    curr=DATAOUT.Rvec.D19;
    currslice = squeeze(curr(50,1:75,:,25));
    slices19.R25 = currslice;
    currslice = squeeze(curr(50,1:75,:,35));
    slices19.R35 = currslice;
    currslice = squeeze(curr(50,1:75,:,50));
    slices19.R50 = currslice;
    currslice = squeeze(curr(50,1:75,:,75));
    slices19.R75 = currslice;

    curr=DATAOUT.Rvec.D49;
    currslice = squeeze(curr(50,1:75,:,25));
    slices49.R25 = currslice;
    currslice = squeeze(curr(50,1:75,:,35));
    slices49.R35 = currslice;
    currslice = squeeze(curr(50,1:75,:,50));
    slices49.R50 = currslice;
    currslice = squeeze(curr(50,1:75,:,75));
    slices49.R75 = currslice;

    curr=DATAOUT.Rstd.D19;
    currslice = squeeze(curr(50,1:75,:,25));
    slices19.Rstd25 = currslice;
    currslice = squeeze(curr(50,1:75,:,50));
    slices19.Rstd50 = currslice;
    currslice = squeeze(curr(50,1:75,:,75));
    slices19.Rstd75 = currslice;

    curr=DATAOUT.Rstd.D49;
    currslice = squeeze(curr(50,1:75,:,25));
    slices49.Rstd25 = currslice;
    currslice = squeeze(curr(50,1:75,:,50));
    slices49.Rstd50 = currslice;
    currslice = squeeze(curr(50,1:75,:,75));
    slices49.Rstd75 = currslice;


    curr=DATAOUT.Rvec.D19;
    currslice = squeeze(curr(50,1:75,:,:));
    slices19.Rvec = currslice;
    curr=DATAOUT.Rstd.D19;
    currslice = squeeze(curr(50,1:75,:,:));
    slices19.Rstd = currslice;

    curr=DATAOUT.Rvec.D49;
    currslice = squeeze(curr(50,1:75,:,:));
    slices49.Rvec = currslice;
    curr=DATAOUT.Rstd.D49;
    currslice = squeeze(curr(50,1:75,:,:));
    slices49.Rstd = currslice;



    [positions] = subplot_pos(8, 3.5, 0.2, 0.2, ...
    0.1, 0.25, 2, 1, 0.1, 0);


    f=figure; f.Color='w';
    f.Units='inches'; f.Position=[4.5 6 8 3.5]; f.Units='pixels';
    axes('position',positions{1});
    imshow(slices19.MeanEAP,[0 2.8]);  colorbar;
    title('Mean EAP | \Delta=19 ms');
    axes('position',positions{2});
    imshow(slices49.MeanEAP,[0 0.6]);  colorbar;
    title('Mean EAP | \Delta=49 ms');

    f=figure; f.Color='w';
    f.Units='inches'; f.Position=[4.5 6 8 3.5]; f.Units='pixels';
    axes('position',positions{1});
    imshow(slices19.MeanGA,[0 0.5]);  colorbar;
    title('Mean GA | \Delta=19 ms');
    axes('position',positions{2});
    imshow(slices49.MeanGA,[0 0.15]);  colorbar;
    title('Mean GA | \Delta=49 ms');

    f=figure; f.Color='w';
    f.Units='inches'; f.Position=[4.5 6 8 3.5]; f.Units='pixels';
    axes('position',positions{1});
    imshow(slices19.MeanPA,[0 3]); colorbar;
    title('Mean PA | \Delta=19 ms');
    axes('position',positions{2});
    imshow(slices49.MeanPA,[0 1]);  colorbar;
    title('Mean PA | \Delta=49 ms');

    f=figure; f.Color='w';
    f.Units='inches'; f.Position=[4.5 6 8 3.5]; f.Units='pixels';
    axes('position',positions{1});
    imshow(slices19.RTOP,[0 8.5]); colorbar;
    title('RTOP | \Delta=19 ms');
    axes('position',positions{2});
    imshow(slices49.RTOP,[0 3]);  colorbar;
    title('RTOP | \Delta=49 ms');

    f=figure; f.Color='w';
    f.Units='inches'; f.Position=[4.5 6 8 3.5]; f.Units='pixels';
    axes('position',positions{1});
    imshow(slices19.FWHM,[4 5.5]); colorbar;
    title('HWHM | \Delta=19 ms');
    axes('position',positions{2});
    imshow(slices49.FWHM,[4 5.5]);  colorbar;
    title('HWHM | \Delta=49 ms');

    f=figure; f.Color='w';
    f.Units='inches'; f.Position=[4.5 6 8 3.5]; f.Units='pixels';
    axes('position',positions{1});
    imshow(slices19.FWHMstd,[0 0.7]); colorbar;
    title('HWHM std | \Delta=19 ms');
    axes('position',positions{2});
    imshow(slices49.FWHMstd,[0 0.7]);  colorbar;
    title('HWHM std | \Delta=49 ms');

    f=figure; f.Color='w';
    f.Units='inches'; f.Position=[4.5 6 8 3.5]; f.Units='pixels';
    axes('position',positions{1});
    imshow(slices19.MSD,[2000 4000]); colorbar;
    title('MSD | \Delta=19 ms');
    axes('position',positions{2});
    imshow(slices49.FWHMstd,[0 0.7]);  colorbar;
    title('HWHM std | \Delta=49 ms');


    f=figure; f.Color='w';
    f.Units='inches'; f.Position=[4.5 6 8 3.5]; f.Units='pixels';
    axes('position',positions{1});
    imshow(slices19.R25,[0 6]); colorbar;
    title('p_{0.25} | \Delta=19 ms');
    axes('position',positions{2});
    imshow(slices49.R25,[0 0.6]);  colorbar;
    title('p_{0.25} | \Delta=49 ms');

    f=figure; f.Color='w';
    f.Units='inches'; f.Position=[4.5 6 8 3.5]; f.Units='pixels';
    axes('position',positions{1});
    imshow(slices19.R50,[0 1]); colorbar;
    title('p_{0.50} | \Delta=19 ms');
    axes('position',positions{2});
    imshow(slices49.R50,[0 0.2]);  colorbar;
    title('p_{0.50} | \Delta=49 ms');

    f=figure; f.Color='w';
    f.Units='inches'; f.Position=[4.5 6 8 3.5]; f.Units='pixels';
    axes('position',positions{1});
    imshow(slices19.R75,[0 0.25]); colorbar;
    title('p_{0.75} | \Delta=19 ms');
    axes('position',positions{2});
    imshow(slices49.R75,[0 0.03]);  colorbar;
    title('p_{0.75} | \Delta=49 ms');

    f=figure; f.Color='w';
    f.Units='inches'; f.Position=[4.5 6 8 3.5]; f.Units='pixels';
    axes('position',positions{1});
    imshow(slices19.R35,[0 3]); colorbar;
    title('p_{0.35} | \Delta=19 ms');
    axes('position',positions{2});
    imshow(slices49.R35,[0 0.25]);  colorbar;
    title('p_{0.35} | \Delta=49 ms');



    f=figure; f.Color='w';
    f.Units='inches'; f.Position=[4.5 6 8 3.5]; f.Units='pixels';
    axes('position',positions{1});
    imshow(slices19.Rstd25,[0 0.5]); colorbar;
    title('a_{0.25} | \Delta=19 ms');
    axes('position',positions{2});
    imshow(slices49.Rstd25,[0 0.3]);  colorbar;
    title('a_{0.25} | \Delta=49 ms');

    f=figure; f.Color='w';
    f.Units='inches'; f.Position=[4.5 6 8 3.5]; f.Units='pixels';
    axes('position',positions{1});
    imshow(slices19.Rstd50,[0 0.9]); colorbar;
    title('a_{0.50} | \Delta=19 ms');
    axes('position',positions{2});
    imshow(slices49.Rstd50,[0 0.12]);  colorbar;
    title('a_{0.50} | \Delta=49 ms');

    f=figure; f.Color='w';
    f.Units='inches'; f.Position=[4.5 6 8 3.5]; f.Units='pixels';
    axes('position',positions{1});
    imshow(slices19.Rstd75,[0 0.5]); colorbar;
    title('a_{0.75} | \Delta=19 ms');
    axes('position',positions{2});
    imshow(slices49.Rstd75,[0 0.1]);  colorbar;
    title('a_{0.75} | \Delta=49 ms');


    fsdir = '/autofs/space/hemera_001/users/rjjones/rsi/susie_invivo/Tract_HC_008/diff/preproc/mri';
    faparc = [fsdir filesep 'aparc_p_aseg_diff.nii.gz'];
    vol=load_untouch_nii(faparc);
    aseg = double(vol.img);
    asegslice = squeeze(aseg(50,:,:));

    fwmmask = [fsdir filesep 'wm_bin.nii.gz'];
    vol=load_untouch_nii(fwmmask);
    wmmask=double(vol.img);
    wmmaskslice = squeeze(wmmask(50,:,:));

    asegslice = asegslice(1:75,:);
    wmmaskslice = wmmaskslice(1:75,:);

    gmmaskslice = asegslice;
    gmmaskslice(wmmaskslice==1)=0;
    gmmaskslice(gmmaskslice>0)=1;


    MeanEAP_wm.D19 = zeros(100,1);
    MeanEAP_wm.D49 = zeros(100,1);
    MeanGAp_wm.D19 = zeros(100,1);
    MeanGAp_wm.D49 = zeros(100,1);
    MeanEAP_gm.D19 = zeros(100,1);
    MeanEAP_gm.D49 = zeros(100,1);
    MeanGAp_gm.D19 = zeros(100,1);
    MeanGAp_gm.D49 = zeros(100,1);

    for ii=1:100
        temp = slices19.Rvec(:,:,ii);
        MeanEAP_wm.D19(ii) = mean(temp(wmmaskslice==1));
        MeanEAP_gm.D19(ii) = mean(temp(gmmaskslice==1));
        temp = slices49.Rvec(:,:,ii);
        MeanEAP_wm.D49(ii) = mean(temp(wmmaskslice==1));
        MeanEAP_gm.D49(ii) = mean(temp(gmmaskslice==1));

        temp = slices19.Rstd(:,:,ii);
        MeanGAp_wm.D19(ii) = mean(temp(wmmaskslice==1));
        MeanGAp_gm.D19(ii) = mean(temp(gmmaskslice==1));
        temp = slices49.Rstd(:,:,ii);
        MeanGAp_wm.D49(ii) = mean(temp(wmmaskslice==1));
        MeanGAp_gm.D49(ii) = mean(temp(gmmaskslice==1));

    end


    figure; hold on;
    plot(MeanEAP_wm.D19,'-b','LineWidth',2.5);
    plot(MeanEAP_gm.D19,'-g','LineWidth',2.5);
    legend('wm','gm'); 
    title('Mean EAP pro. | D19 | TractHC008 slice x=50');

    figure; hold on;
    plot(MeanEAP_wm.D49,'-b','LineWidth',2.5);
    plot(MeanEAP_gm.D49,'-g','LineWidth',2.5);
    legend('wm','gm'); 
    title('Mean EAP pro. | D49 | TractHC008 slice x=50');


    figure; hold on;
    plot(MeanGAp_wm.D19,'-b','LineWidth',2.5);
    plot(MeanGAp_gm.D19,'-g','LineWidth',2.5);
    legend('wm','gm'); 
    title('Mean GA pro. | D19 | TractHC008 slice x=50');

    figure; hold on;
    plot(MeanGAp_wm.D49,'-b','LineWidth',2.5);
    plot(MeanGAp_gm.D49,'-g','LineWidth',2.5);
    legend('wm','gm'); 
    title('Mean GA pro. | D49 | TractHC008 slice x=50');

end


% % Save volumes in .nii format
fafile='/autofs/space/nyx_002/users/rjones/Bay8/CS_DSI_10_01_2021/preproc/dtifit/dtifit_FA.nii.gz';
vol = load_untouch_nii(fafile);
% vol.hdr.dime.datatype = 16;
% vol.hdr.dime.bitpix = 32;


filelist = {'RTOP','MeanEAP','MeanGA','MeanKurt','GaussDiv',...
            'RTOP0','StdEAP','StdEAP0','MeanEAP0','sintheta',...
            'MD','MSD','MFD',...
            'HWHM','HWHMstd','HWHMmin','HWHMmax','HWHMminmax'};
for f=1:length(filelist)
    F = filelist{f};
    eval(['curr=DATAOUT.' F ';']);
    outfile = [INFO.outdir filesep 'nii/' F '.nii.gz'];
    write2nifti( vol, curr, outfile );
end
filelist = {'Rvec','Rstd','KurtProfile'};
for f=1:length(filelist)
    F = filelist{f};
    eval(['curr=DATAOUT.' F ';']);
    outfile = [INFO.outdir filesep 'nii/' F '.nii.gz'];
    write2nifti( vol, curr, outfile );
end


% % INFO.ODFoutdir = [INFO.outdir filesep 'gdsi_odf'];
% % if ~exist(INFO.ODFoutdir,'dir'), mkdir(INFO.ODFoutdir); end
% % 
% % INFO.fvecs = [INFO.procdir filesep 'dtifit/dtifit_V1.nii.gz'];
% if exist(INFO.fvecs,'file')
%     vecvol = load_untouch_nii(INFO.fvecs);
% else
%     disp('COULD NOT FIND VECS VOL!!');
%     vecvol = vol;
%     vecvol.hdr.dime.dim([1 5]) = [4 3];
% end
% 
% fields = fieldnames(DATAODF);
% for f=1:length(fields)
%     F = fields{f};
%     eval(['curr=DATAODF.' F ';']);
%     outfile = [INFO.ODFoutdir filesep F '.nii.gz'];
%     if size(curr,4)==1 
%         write2nifti( vol, curr, outfile );
%     elseif size(curr,4)==3
%         write2nifti( vecvol, curr, outfile );
%     end
% end

% if 0 == 1
% %                     dirmaxs = max(Rmatrix,[],1);
% %                     Rmatrix_n = Rmatrix./repmat(dirmaxs,nr,1);
% %                     tmp = Rmatrix_n<0;
% 
%     tmp2=zeros(nb_pdf_dirs,1);
% 
%     for tt=1:nb_pdf_dirs
% %                         currtmp=tmp(:,tt);
% %                         a=find(currtmp,1,'first');
% 
%         currtmp=Rmatrix(:,tt);
%         a=find(currtmp<0,1,'first');
%         if ~isempty(a)
%             tmp2(tt)=a;
%         else
%             mintmp=find(currtmp==min(currtmp));
%             tmp2(tt)=mintmp(1);
%         end
%     end
%     tmp2(tmp2==0)=1;
%     tmp2b = reshape(displacement_microns(tmp2),nb_pdf_dirs,1);
% 
%     % For a voxel, make pdf contour glyphs
%     pdf_actor = diff_disp_sphere;
%     pdf_actor.vertices = diff_disp_sphere.vertices .* repmat(tmp2b, [1, 3]); % scale radial distance
%     pdf_actor.facevertexcdata = tmp2b; % change color data to represent pdf values
%     if min(pdf_actor.faces(:))==0, pdf_actor.faces=pdf_actor.faces+1; end
% 
%     % display pdf contour
%     ff=figure('color','w','position',[288 415 855 460]);
%     ff.InvertHardcopy='off';
% %                 subplot(1,2,1);
%     h = patch(pdf_actor);
%     view(180,0); %view(90, 90);
%     lighting gouraud; shading faceted; camlight
%     set(h, 'EdgeColor', 'none');
%     % colormap hsv; %colorbar;
%     caxis([min(tmp2b(:)), max(tmp2b(:))]);
%     axis equal, axis off, axis tight
% %                 title(['pdf contour ' sprintf('%.03f',r0) 'xMDD']);
%     ax1=gca;
% end
                
                

% disp(' --- Finished reconstructing voxels --- ');
% ftoc=toc(ftic)
% 
% 
% disp(' - Saving output .mat file...');
% 
% if ~exist(INFO.outdir,'dir')
%     fprintf('Making outdir:\n %s\n',INFO.outdir);
%     mkdir(INFO.outdir);
% end
% 
% foutmat = [INFO.outdir filesep 'gdsi_recon_v2.mat'];
% fprintf('Saving recon to file:\n %s\n',foutmat);
% save(foutmat, 'DATAOUT', 'INFO', 'ACQ', 'pdf_dirs','nr','displacement_microns', 'r_step_size', 'ftoc', '-v7.3');
% 
% 
% fprintf('\n-------------\n - Finished GDSI recon - \n-------------\n');


% end