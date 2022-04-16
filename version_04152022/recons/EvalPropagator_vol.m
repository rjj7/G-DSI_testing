function [ DATAOUT ] = EvalPropagator_vol( DATAIN, DATAOUT, opts ) 
% 
% [ DATAOUT ] = EvalPropagator_vol( DATAIN, params, opts ) 
%
% DATAIN has:   dwi, mask, index_b0, qvol, pdf_dirs, isshelled, nr, ndirs
% DATAOUT has:  empty matrices to store measures/stats from recon
% opts has:     StoreRmatrix, SaveDATAOUT
% 
% 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%    RJ  02 17 2022   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

verbose = true;

if verbose
    disp(' --- Starting GDSI recon --- ');
end

totalvoxels = nnz(DATAIN.mask);
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
                b0s = data(DATAIN.index_b0);       % b0 volumes
                dwis = data(~DATAIN.index_b0);     % dwi volumes
                data = [mean(b0s); dwis];   % average b0 volumes
                data = data / mean(b0s);    % normalize diff signal
                
                % metrics on q-space signal 
                DATAOUT.qiv(ix,iy,iz) = inv(sum(data));
%                 DATAOUT.qmsd(ix,iy,iz) = mean(data.*(mag_qvec.^2));
                DATAOUT.qvar(ix,iy,iz) = var(data);
                DATAOUT.qentropy(ix,iy,iz) = entropy(data);

                %%% Evaluate EAP               
                if DATAIN.isshelled
                    data = diag(DATAIN.qvol) * data;
                end
                Rmatrix = DATAIN.Fmtx * data;
                
%                 [ Rmatrix ] = compute_1d_eap_mtx( data, Fmtx );

                Rmatrix = reshape(Rmatrix, DATAIN.nr, DATAIN.ndirs);
                
                if opts.StoreRmatrix
                    DATAOUT.Rmatrix(ix,iy,iz,:,:) = Rmatrix;
                end
                
                Rvec        = mean(Rmatrix,2);
                Rstd        = std(Rmatrix,0,2);
                Rrms        = rms(Rmatrix,2);
                Rsigma      = std(Rmatrix(:));   %std of eap
                
                
                DATAOUT.RTOP(ix,iy,iz)             = Rvec(1); %Rmatrix(1,1);      %rtop of eap
                
                DATAOUT.Rvec(ix,iy,iz,:)           = Rvec;      
                DATAOUT.Rstd(ix,iy,iz,:)           = Rstd;       
                DATAOUT.Rrms(ix,iy,iz,:)           = Rrms;
                
                DATAOUT.MeanGA(ix,iy,iz)           = mean(Rstd); 
                DATAOUT.EntropyEAP(ix,iy,iz)       = entropy(Rmatrix(:));
                
%                 EAPEntropyProfile = zeros(nr,1);
%                 for rr=1:nr
%                     Rm = Rmatrix(rr,:);
%                     Rentropy = entropy(Rm);
%                     EAPEntropyProfile(rr) = Rentropy;
%                 end

                
                %%%%%%%% (NON)GAUSSIAN METRICS  %%%%%%%%
                % Divergence of Mean 1D EAP from Gaussian
                gaussmtx = exp(-DATAIN.dispvec2'*Rsigma);        %make gaussian
                gaussmtx = (gaussmtx/sum(gaussmtx))*sum(Rvec(:));
                                
                if min(Rvec)<0
                    rtmp = Rvec-min(Rvec);
                else
                    rtmp = Rvec;
                end
%                 if min(Rmatrix)<0
%                     mtmp = Rmatrix-min(Rmatrix(:));
%                 else
%                     mtmp = Rmatrix;
%                 end
                
                linearmtx = linspace(max(rtmp),0,DATAIN.nr)';
                
                isomtx = mean(rtmp)*ones(size(rtmp));
%                 isomtx = repmat(rtmp,1,size(Rmatrix,2));
                
%                 (KLDiv(mtmp(1:50,:),isomtx(1:50,:)));   %compute KL divergence
                
                tmp = (KLDiv(rtmp(1:50)',gaussmtx(1:50)'));
                tmp(isinf(tmp))=0;
                tmp(isnan(tmp))=0;
                DATAOUT.KLDg(ix,iy,iz)   = tmp;   %compute KL divergence
                
                tmp = (KLDiv(rtmp(1:50)',isomtx(1:50)'));
                tmp(isinf(tmp))=0;
                tmp(isnan(tmp))=0;
                DATAOUT.KLDi(ix,iy,iz)   = tmp;   %compute KL divergence
                
                DATAOUT.NLrms(ix,iy,iz) = sqrt(mean((rtmp-linearmtx).^2));
                
                % Calculate non-Gaussianity (NG) [Ning 2015]
                costheta = dot(Rvec,gaussmtx)/sqrt(dot(Rvec,Rvec)*dot(gaussmtx,gaussmtx));
                DATAOUT.NG(ix,iy,iz) = sqrt(1-costheta^2); %really sintheta
                
                %%%%%%%% KURTOSIS %%%%%%%%
                % Calculate kurtosis at each distance r
%                 DATAOUT.RotKurt(ix,iy,iz,:) = kurtosis(Rmatrix');
                
                % Calculate average eap kurtosis
                DATAOUT.MeanKurt(ix,iy,iz) = kurtosis(Rmatrix(:));
                
%                 % Mean squared displacement
%                 DATAOUT.MSD(ix,iy,iz) = sum(dispvec2' .* Rvec);
%                 DATAOUT.MFD(ix,iy,iz) = sum(dispvec4' .* Rvec);
                
                % HWHM radius (only half width, not full!)
                [ hwhm ] = compute_ralpha( Rmatrix, 0.5, DATAIN.dispvec );
                DATAOUT.HWHM(ix,iy,iz) = hwhm.mean;
                DATAOUT.HWHMstd(ix,iy,iz) = hwhm.std;
                DATAOUT.HWHMmin(ix,iy,iz) = hwhm.min;
                DATAOUT.HWHMmax(ix,iy,iz) = hwhm.max;
                
                tmp_mvec = DATAIN.pdf_dirs((hwhm.vec==hwhm.min),:);
                if size(tmp_mvec,1)>1
                    tmp_mvec = tmp_mvec(1,:);
                end
                DATAOUT.HWHMmin_vec(ix,iy,iz,:) = tmp_mvec;
                
                tmp_mvec = DATAIN.pdf_dirs((hwhm.vec==hwhm.max),:);
                if size(tmp_mvec,1)>1
                    tmp_mvec = tmp_mvec(1,:);
                end
                DATAOUT.HWHMmax_vec(ix,iy,iz,:) = tmp_mvec;
                
%                 if opts.ClipRmatrix
%                     % Clip after first zero crossing
%                     [ nRmatrix ] = clip_negative_Rmatrix( Rmatrix );
% %                     [ rRmatrix ] = clip_ringing_Rmatrix( Rmatrix );
% %                     [ nrRmatrix ] = clip_negative_Rmatrix( rRmatrix );
%                     nRvec = mean(nRmatrix, 2);
%                     nRstd = std( nRmatrix, 0, 2);
% 
%                     if opts.MakePlots.rgbprofiles
%                         [ f1 ] = plot_all_profiles_voxel( dispvec, Rmatrix, ...
%                             abs(diff_disp_sphere.vertices)' );
%                          hold on; plot(dispvec,zeros(size(dispvec)),'-k');
%                          
%                         pause(1);
%                         
%                         [ f1 ] = plot_all_profiles_voxel( dispvec, nRmatrix, ...
%                             abs(diff_disp_sphere.vertices)' ); 
%                          hold on; plot(dispvec,zeros(size(dispvec)),'-k');
%                     end
%                 end
%                 
% 
%                 
%                 if ReconODFs
%                     %%%%%%%% ODF %%%%%%%%
%                     % Recon dODF
%                     Omatrix = Rmatrix.*odfpar.Wmtx';
%                     if Use_DSI_Studio_dirs
%                         Omatrix = repmat(Omatrix,1,2);
%                     end                
%                     ODF = sum( Omatrix( odfpar.idx_start:odfpar.idx_end, :) );
%                     ODF_points = diff_disp_sphere.vertices .* repmat(ODF', [1, 3]);
%                     
%                     [ f1 ] = plot_all_profiles_voxel( dispvec, Omatrix, ...
%                             abs(diff_disp_sphere.vertices)' ); 
%                          hold on; plot(dispvec,zeros(size(dispvec)),'-k');
%                          
%                          
%                          
%                     nOmatrix = nRmatrix.*odfpar.Wmtx';
%                     if Use_DSI_Studio_dirs
%                         nOmatrix = repmat(nOmatrix,1,2);
%                     end                
%                     nODF = sum( nOmatrix ); %( odfpar.idx_start:odfpar.idx_end, :) );
%                     nODF_points = diff_disp_sphere.vertices .* repmat(nODF', [1, 3]);
%                     
%                     [ f1 ] = plot_all_profiles_voxel( dispvec, nOmatrix, ...
%                             abs(diff_disp_sphere.vertices)' ); 
%                          hold on; plot(dispvec,zeros(size(dispvec)),'-k');
%                          
%                     [pks2, allpks2] = gdsi_find_peak(nODF', diff_disp_sphere, odfpar.sep);     
%                          
% 
%                     % find ODF peaks
%                     [pks, allpks] = gdsi_find_peak(ODF', diff_disp_sphere, odfpar.sep);
%                     npeaks.D19 = length(pks.D19.pkvals);
%                     totalpeaks.D19 = npeaks.D19;
%                     if npeaks.D19>3                    
%                         npeaks.D19=3;
%                     end
%                     
%                     %%% Compute RTAP, RTPP on major fiber axis
%                     pdd = pks.pkvecs(1,:);
%                     [ rtap ] = compute_RTAP( data, v, pdd, qvec, rs, dispvec2 );
%                     [ rtpp ] = compute_RTPP( data, v, a_plane, pdd, qvec, rs, dispvec2 );
% 
%                     %%% Compute fixel metrics
%                     [fixels] = compute_fixel_metrics_v2( pks, diff_disp_sphere.vertices, ...
%                         Rmatrix, data, dispvec, zdps); 
% 
% 
%                     if MakePlots.odf
%                         %%%%% plot ODF
%     %                     plotODF = ODF;
%                         plotODF = ODF/max(ODF);
%         % %                 plotODF = (ODF-min(ODF))/(max(ODF)-min(ODF));
%                         plotPeaks = pks.pkvecs;
% 
%                         [ f ] = gdsi_plot_odf( diff_disp_sphere, plotODF, plotPeaks );
%                         
%                         
%                         
%                         
%                         plotODF = nODF/max(nODF);
%         % %                 plotODF = (ODF-min(ODF))/(max(ODF)-min(ODF));
%                         plotPeaks = pks2.pkvecs;
% 
%                         [ f ] = gdsi_plot_odf( diff_disp_sphere, plotODF, plotPeaks );
%                         
%                         
%                     end
%                 
%                     %%% STORE RESULTS
%                     % ODF GFA, RMS, STD, npeaks, iso, mass, max
%                     DATAODF.GFA(ix,iy,iz) = std(nODFhalf)/rms(nODFhalf);
%                     DATAODF.RMS(ix,iy,iz) = rms(nODFhalf);
%                     DATAODF.STD(ix,iy,iz) = std(nODFhalf);
%                     DATAODF.npeaks(ix,iy,iz) = n_totalpeaks;
%                     % ODF iso (min value)
%                     DATAODF.iso(ix,iy,iz) = min(nODF);
%                     % ODF mass (sum of all odf values)
%                     DATAODF.mass(ix,iy,iz) = sum(nODF);
% 
%                     %%% Fixels
%                     for nn=1:n_npeaks
%                         eval(['curr = n_fixels.v' num2str(nn) ';']);
%                         currind = curr.diridx;
%                         currO = nOmatrix(odfpar.idx_start:odfpar.idx_end, currind);
%                         currfa = sum(currO);
%     %                     currnfa = currfa/sum(currO);
% 
%                         eval(['DATAODF.ind' num2str(nn) '(ix,iy,iz) = curr.diridx;']);
%                         eval(['DATAODF.vec' num2str(nn) '(ix,iy,iz,:) = curr.dir;']);
% 
%                         eval(['DATAODF.f' num2str(nn) '(ix,iy,iz) = curr.f;']);
%                         eval(['DATAODF.fa' num2str(nn) '(ix,iy,iz) = curr.fa;']);
% 
%                         eval(['DATAODF.Rvec' num2str(nn) '(ix,iy,iz,:) = curr.Rvec;']);
%                         eval(['DATAODF.rtap' num2str(nn) '(ix,iy,iz) = curr.Rmean;']);
%                         eval(['DATAODF.Rstd' num2str(nn) '(ix,iy,iz) = curr.Rstd;']);
%                         eval(['DATAODF.Rmsd' num2str(nn) '(ix,iy,iz) = curr.Rmsd;']);            
%                         eval(['DATAODF.kurt' num2str(nn) '(ix,iy,iz) = curr.Kurtosis;']);                    
%                         eval(['DATAODF.ip' num2str(nn) '(ix,iy,iz) = curr.inflpnt;']);
%                         eval(['DATAODF.r_50_' num2str(nn) '(ix,iy,iz) = curr.r_50;']);
%                         eval(['DATAODF.r_25_' num2str(nn) '(ix,iy,iz) = curr.r_25;']);
% 
%                         eval(['DATAODF.rtpp.Rvec' num2str(nn) '(ix,iy,iz,:) = curr.rtpp.Rvec;']);                    
%                         eval(['DATAODF.rtpp.Rstd' num2str(nn) '(ix,iy,iz,:) = curr.rtpp.Rstd;']);                    
%                         eval(['DATAODF.rtpp.rtpp' num2str(nn) '(ix,iy,iz,:) = curr.rtpp.RTPP;']);                    
%                         eval(['DATAODF.rtpp.rtpv' num2str(nn) '(ix,iy,iz,:) = curr.rtpp.RTPV;']);                    
%                         eval(['DATAODF.rtpp.msd' num2str(nn) '(ix,iy,iz,:) = curr.rtpp.msd;']);                    
%                         eval(['DATAODF.rtpp.kurt' num2str(nn) '(ix,iy,iz,:) = curr.rtpp.kurtosis;']);                    
% 
%                     end
%                 end
                   
                cnt=cnt+1;
            end
        end
    end
end

disp(' --- Finished reconstructing voxels --- ');
ftoc=toc(ftic);
fprintf('%.1f seconds | %d voxels\n',ftoc,totalvoxels);

DATAOUT.recon_vox_count = cnt;
DATAOUT.recon_elap_time = ftoc;

if opts.SaveDATAOUT
    disp(' - Saving output .mat file...');
    
    if ~exist(DATAIN.outdir,'dir'), mkdir(DATAIN.outdir); end
    foutmat = [DATAIN.outdir filesep 'gdsi_vol_recon_prof+scal.mat'];

    fprintf('Saving recon to file:\n %s\n',foutmat);
    save(foutmat, 'DATAOUT', '-v7.3');
    fprintf('\n-------------\n - Finished saving GDSI recon - \n-------------\n');
end



