function [ F ] = precompute_Fmtx( qvec, pdf_dirs, rs )
% 
% [ F ] = precompute_Fmtx( qvec, pdf_dirs, rs )
% 
%  Compute F matrix for GDSI recon. Dimension is [R*D,N], where R is nb of
%  (radial) displacements, D is nb of (angular) directions, and N is nb of
%  qspace samples.
% 

% From: Q. Tian et al. (NeuroImage, 2019), Equation 5

if nnz(size(pdf_dirs)==3)==0
    error('Invalid PDF directions [req: x, y, z]');
end
if size(pdf_dirs,1)==3
    pdf_dirs = pdf_dirs';
end

R = length(rs);
N = length(qvec);
D = size(pdf_dirs,1);
F = zeros(R*D,N);

for ind=1:D
    pstart = (ind-1)*R + 1;
    pend = pstart + R - 1;
    pdf_dir = pdf_dirs(ind,:);
    rvec = repmat(pdf_dir, [R, 1]) .* repmat(rs, [1, 3]);
    f = cos(rvec * qvec') / N;
    F(pstart:pend,:) = f;
end

disp(' ');
