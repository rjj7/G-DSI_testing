function [ Rmatrix ] = compute_1d_eap_mtx( data, Fmtx )
%
% [ Rmatrix ] = compute_1d_eap_mtx( data, Fmtx )
% 
%  Rmatrix = matrix with all EAP values at all distances + for all dirs
%  data    = input dwi data from one voxel (Assumes DSI grid data)
%  Fmtx    = 3D F matrix made with precompute_Fmtx()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%nb of pdf dirs on sphere
nb_dirs = size(Fmtx,3);

%nb of radial displacement distances
nb_r = size(Fmtx,1);

%Pre-allocate storage
Rmatrix = zeros(nb_r,nb_dirs); %to store all EAP prob. dens. 

% Loop pdf_dirs, eval pdf at radial points on each dir
for dirind=1:nb_dirs
    F = Fmtx(:,:,dirind);   %get F for current dir
    pdf_1d = F * data;      %calculate 1D PDF
    Rmatrix(:,dirind) = pdf_1d; %store 1D PDF
end

end



