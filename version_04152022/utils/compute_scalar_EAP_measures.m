function [ DATAOUT ] = compute_scalar_EAP_measures( DATAOUT, mask, outdims, gparams )


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Compute metrics from results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Calculate measures from the spherically averaged profiles
%%%  saved in DATAOUT struct. 
%%%  - mean of "Std. Dev. of EAP" profile
%%%  - 1st deriv of mean EAP profile
%%%  - 2nd deriv of mean EAP profile
%%%  - min of 1st deriv.
%%%  - inflection point of 2nd deriv
%%%  - "rtop" of mean EAP profile, its 1st deriv and 2nd deriv
%%%  - kurtosis of mean EAP profile
%%%  - KL divergence from isotropic of mean EAP profile

% DATAOUT = rmfield(DATAOUT,'scalarmaps');

% 1st derivative of mean 1D EAP
DATAOUT.scalarmaps.diff1_avg_disp_sph = diff(DATAOUT.avg_disp_sph,1,4)/gparams.r_step_size;

% 2nd derivative of mean 1D EAP
DATAOUT.scalarmaps.diff2_avg_disp_sph = diff(DATAOUT.avg_disp_sph,2,4)/gparams.r_step_size;

% Allocate
DATAOUT.scalarmaps.diff1_min = zeros(outdims(1),outdims(2),outdims(3));
DATAOUT.scalarmaps.std_disp_sph_mean = zeros(outdims(1),outdims(2),outdims(3));
DATAOUT.scalarmaps.inflection_point = zeros(outdims(1),outdims(2),outdims(3));
DATAOUT.scalarmaps.rtop1 = zeros(outdims(1),outdims(2),outdims(3));
DATAOUT.scalarmaps.rtop2 = zeros(outdims(1),outdims(2),outdims(3));
DATAOUT.scalarmaps.rtop3 = zeros(outdims(1),outdims(2),outdims(3));
DATAOUT.scalarmaps.kurtosis = zeros(outdims(1),outdims(2),outdims(3));
DATAOUT.scalarmaps.gausskldiv = zeros(outdims(1),outdims(2),outdims(3));
DATAOUT.scalarmaps.isokldiv = zeros(outdims(1),outdims(2),outdims(3));

% Create gaussian pdf (not sure about this?)
xgauss = linspace(0,3,gparams.nr);
ygauss = normpdf(xgauss,0,1);
ygauss0 = ygauss/max(ygauss);

%Find nonzero voxels, analyze them iteratively
nzinds = find(mask);
for nzid=1:length(nzinds)
    nzind = nzinds(nzid);
    %get subscript coordinate from index
    [nx,ny,nz] = ind2sub(size(mask),nzind);
    
    %calculate the inflection point of the 2nd deriv profile
    temp = reshape(DATAOUT.scalarmaps.diff2_avg_disp_sph(nx,ny,nz,:),[],gparams.nr-2);
    tempind = find(temp>0,1,'first');
    DATAOUT.scalarmaps.inflection_point(nx,ny,nz) = tempind*gparams.r_step_size;
    
    %calculate the RTOP of mean EAP, its 1st deriv, and its 2nd deriv.
    DATAOUT.scalarmaps.rtop1(nx,ny,nz) = DATAOUT.avg_disp_sph(nx,ny,nz,1);
    DATAOUT.scalarmaps.rtop2(nx,ny,nz) = DATAOUT.scalarmaps.diff1_avg_disp_sph(nx,ny,nz,1);
    DATAOUT.scalarmaps.rtop3(nx,ny,nz) = DATAOUT.scalarmaps.diff2_avg_disp_sph(nx,ny,nz,1);
    
    %calculate the min of the 1st derivative of mean EAP profile
    temp = reshape(DATAOUT.scalarmaps.diff1_avg_disp_sph(nx,ny,nz,:),[],gparams.nr-1);
    DATAOUT.scalarmaps.diff1_min(nx,ny,nz) = min(temp);
    
    %calculate kurtosis of mean EAP profile
    temp = reshape(DATAOUT.avg_disp_sph(nx,ny,nz,:),[],gparams.nr);
    DATAOUT.scalarmaps.kurtosis(nx,ny,nz) = kurtosis(temp);
    
    %calculate KL divergence from isotropic pdf of mean EAP profile
    sfact = max(temp);
    ygauss_ = ygauss0*sfact;
    DATAOUT.scalarmaps.gausskldiv(nx,ny,nz) = real(KLDiv(ygauss_,temp));
    
    %calculate KL divergence from isotropic pdf of mean EAP profile
    sfact = max(temp);
    yiso_ = ones(size(temp))*sfact;
    DATAOUT.scalarmaps.isokldiv(nx,ny,nz) = real(KLDiv(yiso_,temp));
    
    %calculate mean of the std dev eap profile
    temp = reshape(DATAOUT.std_disp_sph(nx,ny,nz,:),[],gparams.nr);
    DATAOUT.scalarmaps.std_disp_sph_mean(nx,ny,nz) = mean(temp);
    
end


end