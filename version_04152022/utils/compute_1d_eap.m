function [ Rmean, Rstd, Rmatrix, Rmean_nneg ] = compute_1d_eap( data, Fmtx, do_nneg )
%
% [ Rmean, Rstd, Rmatrix, Rmean_nneg ] = compute_1d_eap( data, Fmtx, do_nneg )
% 
% Compute GDSI metrics on 1D EAP displacement profiles averaged across all
% the (EAP) points on the sphere. 
% For each direction on the sphere, split the vector into nr equally spaced
% points, and evaluate the EAP at these points. Do this for each direction,
% and average across all directions.
%
%  Rmean = spherically averaged probability density at each distance
%  Rstd  = spherically averaged std. dev. of prob. densities along each dir
%  Rmatrix = matrix with all EAP values at all distances + for all dirs
%  Rmean_nneg = Rmean but clipping PDF after first zero-crossing
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%nb of pdf dirs on sphere
nb_dirs = size(Fmtx,3);

%nb of radial displacement distances
nb_r = size(Fmtx,1);

%Pre-allocate storage
Rvec_ = zeros(nb_r,1); %
Rvec_std_ = 0;         % Add to these, then average after looping all dirs
if do_nneg
    Rvec_nneg_ = zeros(nb_r,1);
end
Rmatrix = zeros(nb_r,nb_dirs); %to store all EAP prob. dens. 

% Loop pdf_dirs, eval pdf at radial points on each dir
for dirind=1:nb_dirs
    F = Fmtx(:,:,dirind);   %get F for current dir
    pdf_1d = F * data;      %calculate 1D PDF
    Rvec_ = Rvec_ + pdf_1d;   %add 1D PDF
    Rvec_std_ = Rvec_std_ + std(pdf_1d); %add 1D PDF std dev
    Rmatrix(:,dirind) = pdf_1d; %store 1D PDF
    if do_nneg
        pdf_nneg = pdf_1d;      %to calc non-neg PDF
        ind_negative = find(pdf_nneg < 0);
        if ~isempty(ind_negative)   %to zero after 1st 0-cross
            pdf_nneg(ind_negative(1):end) = 0;
        end
        Rvec_nneg_ = Rvec_nneg_ + pdf_nneg; %store nneg PDF
    end
end


%%%% Average across dirs (spherical average)
%sph average of mean prob. dens.
Rmean = Rvec_/nb_dirs;
%sph average of std deviation of prob dens
Rstd    = Rvec_std_/nb_dirs;
%sph average of nonnegative mean prob. dens.
if do_nneg
    Rmean_nneg = Rvec_nneg_/nb_dirs;
else
    Rmean_nneg = [];
end


end



