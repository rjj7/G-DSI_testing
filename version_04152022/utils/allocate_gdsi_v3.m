function [ DATAOUT, DATAODF ] = allocate_gdsi_v4( voldims, nr )
%
% [ DATAOUT, DATAODF ] = allocate_gdsi_v4( voldims, nr, nr_odf )
% 

if nargin<2
    disp('USAGE: [ DATAOUT, DATAODF ] = allocate_gdsi_v3( voldims, nr )');
    return
end


% Define structs
DATAOUT = [];  DATAODF = [];


% Mean EAP r-profile (spherically averaged)
DATAOUT.Rvec        = zeros([voldims nr]);
% Generalized anisotropy r-profile (spherically averaged) 
DATAOUT.Rstd        = zeros([voldims nr]);

% The mean of the GA profile
DATAOUT.MeanGA              = zeros(voldims);

% the mean + std dev of entire Pr
DATAOUT.StdEAP              = zeros(voldims);
DATAOUT.MeanEAP             = zeros(voldims);

% the mean + std dev before PDF normalization
DATAOUT.RmatrixSum0         = zeros(voldims);

% the mean + std dev before PDF normalization
DATAOUT.StdEAP0              = zeros(voldims);
DATAOUT.MeanEAP0             = zeros(voldims);


%Mean kurtosis of EAP 
DATAOUT.MeanKurt            = zeros(voldims);
%Kurtosis profile in r-space
DATAOUT.KurtProfile         = zeros([voldims nr]);

%Divergence of EAP from gaussian
DATAOUT.GaussDiv            = zeros(voldims);

%Non gaussianity (Ning et al 2015 IEEE)
DATAOUT.costheta            = zeros(voldims);
DATAOUT.sintheta            = zeros(voldims);
DATAOUT.NG                  = zeros(voldims);

%mean displacement, mean sq displacement, mean fourth-order displacement
DATAOUT.MD             = zeros(voldims);
DATAOUT.MSD            = zeros(voldims);
% DATAOUT.MFD            = zeros(voldims);

%FWHM
DATAOUT.HWHM = zeros(voldims);
DATAOUT.HWHMstd = zeros(voldims);

%RTOP
DATAOUT.RTOP = zeros(voldims);
DATAOUT.RTOP0 = zeros(voldims);
% DATAOUT.RTAP = zeros(voldims);
% DATAOUT.RTAV = zeros(voldims);
% DATAOUT.RTAP_vec = zeros([voldims nr]);
% DATAOUT.RTAV_vec = zeros([voldims nr]);
% DATAOUT.RTAP_msd = zeros(voldims);
% DATAOUT.RTPP = zeros(voldims);


% %%%%% ODFs
% 
% % ODF GFA
% DATAODF.GFA = zeros(voldims);
% DATAODF.RMS = zeros(voldims);
% DATAODF.STD = zeros(voldims);
% % ODF iso (min value)
% DATAODF.iso = zeros(voldims);
% % ODF max value
% DATAODF.max = zeros(voldims);
% % ODF mass (sum of all odf values)
% DATAODF.mass = zeros(voldims);
% %NbPeaks
% DATAODF.npeaks = zeros(voldims);
% 
% 
% %%% ODF - FIXEL METRICS
% % peak index in odf_vertices
% DATAODF.ind1 = zeros(voldims);
% DATAODF.ind2 = zeros(voldims);
% DATAODF.ind3 = zeros(voldims);
% % peak vectors/direction
% DATAODF.vec1 = zeros([voldims 3]);
% DATAODF.vec2 = zeros([voldims 3]);
% DATAODF.vec3 = zeros([voldims 3]);
% 
% % F (peak vol. frac. norm s.t. F(peak1)==1
% DATAODF.f1 = zeros(voldims);
% DATAODF.f2 = zeros(voldims);
% DATAODF.f3 = zeros(voldims);
% % FA (OFD peak vol. frac.)
% DATAODF.fa1 = zeros(voldims);
% DATAODF.fa2 = zeros(voldims);
% DATAODF.fa3 = zeros(voldims);
% % EAP FA (RTAP, peak vol. frac.)
% DATAODF.rtap1 = zeros(voldims);
% DATAODF.rtap2 = zeros(voldims);
% DATAODF.rtap3 = zeros(voldims);
% 
% % mean squared displacement of 1d eap profile
% DATAODF.msd1 = zeros(voldims);
% DATAODF.msd2 = zeros(voldims);
% DATAODF.msd3 = zeros(voldims);
% 
% % eap 1d profile
% DATAODF.Rvec1 = zeros([voldims nr]);
% DATAODF.Rvec2 = zeros([voldims nr]);
% DATAODF.Rvec3 = zeros([voldims nr]);
% % std of eap 1d profile
% DATAODF.Rstd1 = zeros([voldims ]);
% DATAODF.Rstd2 = zeros([voldims ]);
% DATAODF.Rstd3 = zeros([voldims ]);
% % sum of eap 1d profile
% DATAODF.Rsum1 = zeros([voldims ]);
% DATAODF.Rsum2 = zeros([voldims ]);
% DATAODF.Rsum3 = zeros([voldims ]);
% % fwhm of eap 1d profile
% DATAODF.fwhm1 = zeros(voldims);
% DATAODF.fwhm2 = zeros(voldims);
% DATAODF.fwhm3 = zeros(voldims);
% % mean kurtosis of eap 1d profile
% DATAODF.kurt1 = zeros(voldims);
% DATAODF.kurt2 = zeros(voldims);
% DATAODF.kurt3 = zeros(voldims);
% % inflection point of eap 1d profile
% DATAODF.ip1 = zeros(voldims);
% DATAODF.ip2 = zeros(voldims);
% DATAODF.ip3 = zeros(voldims);
% 
% 
% % r_alpha maps (disp dist where EAP prob dens = r_alpha)
% DATAODF.r_50_1 = zeros(voldims);
% DATAODF.r_50_2 = zeros(voldims);
% DATAODF.r_50_3 = zeros(voldims);
% DATAODF.r_25_1 = zeros(voldims);
% DATAODF.r_25_2 = zeros(voldims);
% DATAODF.r_25_3 = zeros(voldims);
% 
% 
% % RTPP fixel 1D vector/profile
% DATAODF.rtpp.Rvec1 = zeros([voldims nr]);
% DATAODF.rtpp.Rvec2 = zeros([voldims nr]);
% DATAODF.rtpp.Rvec3 = zeros([voldims nr]);
% 
% % Fixel RTP 1D anisotropy vector/profile
% DATAODF.rtpp.Rstd1 = zeros([voldims nr]);
% DATAODF.rtpp.Rstd2 = zeros([voldims nr]);
% DATAODF.rtpp.Rstd3 = zeros([voldims nr]);
% 
% % Fixel RTPP
% DATAODF.rtpp.rtpp1 = zeros(voldims);
% DATAODF.rtpp.rtpp2 = zeros(voldims);
% DATAODF.rtpp.rtpp3 = zeros(voldims);
% 
% % Fixel RTP variance (std dev of avg RTPP vector
% DATAODF.rtpp.rtpv1 = zeros(voldims);
% DATAODF.rtpp.rtpv2 = zeros(voldims);
% DATAODF.rtpp.rtpv3 = zeros(voldims);
% 
% % Fixel RTPP MSD
% DATAODF.rtpp.msd1 = zeros(voldims);
% DATAODF.rtpp.msd2 = zeros(voldims);
% DATAODF.rtpp.msd3 = zeros(voldims);
% 
% % Fixel RTPP kurtosis
% DATAODF.rtpp.kurt1 = zeros(voldims);
% DATAODF.rtpp.kurt2 = zeros(voldims);
% DATAODF.rtpp.kurt3 = zeros(voldims);


end