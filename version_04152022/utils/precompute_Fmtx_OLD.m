function [ Fmtx ] = precompute_Fmtx( qvec, pdf_dirs, rs )

% From: Q. Tian et al. (NeuroImage, 2019), Equation 5

nr = length(rs);

if nnz(size(pdf_dirs)==3)==0
    error('Invalid directions');
end

if size(pdf_dirs,1)==3
    pdf_dirs = pdf_dirs';
end

ndirs = size(pdf_dirs,1);
Fmtx = zeros(nr,length(qvec),ndirs);

for pdfind=1:ndirs
    pdf_dir = pdf_dirs(pdfind,:);
    rvec = repmat(pdf_dir, [nr, 1]) .* repmat(rs, [1, 3]);
    F = cos(rvec * qvec') / length(qvec);
    Fmtx(:,:,pdfind) = F;
end
 

end
