function [ out ] = compute_1d_gdsi_metrics( data, Fmtx, r_step_size_micron, do_nneg )
%
% [ out ] = compute_1d_gdsi_metrics( data, Fmtx, r_step_size_micron[, do_nneg] )
% 
% Compute GDSI metrics on 1D EAP displacement profiles averaged across all
% the (EAP) points on the sphere. 
% For each direction on the sphere, split the vector into nr equally spaced
% points, and evaluate the EAP at these points. Do this for each direction,
% and average across all directions.
%
% USAGE: 
%   data                = DSI DWI data from one voxel
%   Fmtx                = precomputed F matrix containing qvec*rvec'
%   r_step_size_micron  = distance b/w points in radial vector
%   do_nneg [opt]       = compute non-negative pdf metrics (def. false)
%                           ( use true to compute them)
% 
% OUTPUTS:
%   out.Rvec   = metrics on 1D EAP profiles (sph. avg. at each r)
%   out.dRvec  = metrics on 1st deriv. of EAP profiles
%   out.ddRvec = metrics on 2nd deriv. of EAP profiles
%
% 
% In the code: (might have changed....) 
%   - "Rvec" is a 1xnr vector, where nr is the number of displacement
%       distances. It is the directionally-averaged EAP at each displacement
%       distance. 
%   - "Rvec_auc" the sum of Rvec. 
%   - "Rvec_nneg" is Rvec but clipped (to 0) after first zero-crossing
%   - "Rvec_std" is the average standard deviation of the 1D EAP
%       Calculate std(pdf_1d) at each direction, then average. 



% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%% ANALYZE RADIAL POINTS ALONG EACH DIRECTION
% 

nb_dirs = size(Fmtx,3);
nb_r = size(Fmtx,1);

Rvec_ = zeros(nb_r,1);
Rvec_std_ = 0;
Rvec_nneg = zeros(nb_r,1);

Rmatrix = zeros(nb_r,nb_dirs);

% Loop pdf_dirs, eval pdf at radial points on each dir
for dirind=1:nb_dirs
    F = Fmtx(:,:,dirind);   %get F for current dir
    pdf_1d = F * data;      %calculate 1D PDF
    Rvec_ = Rvec_ + pdf_1d;   %store 1D PDF
    Rvec_std_ = Rvec_std_ + std(pdf_1d);
    Rmatrix(:,dirind) = pdf_1d;
    if do_nneg
        pdf_nneg = pdf_1d;      %to calc non-neg PDF
        ind_negative = find(pdf_nneg < 0);
        if ~isempty(ind_negative)   %to zero after 1st 0-cross
            pdf_nneg(ind_negative(1):end) = 0;
        end
        Rvec_nneg = Rvec_nneg + pdf_nneg; %store nneg PDF
    end
end


%%%% Average 1D profile across dirs
%mean
Rvec.mean = Rvec_/nb_dirs;
%std deviation
Rvec.std    = Rvec_std_/nb_dirs;
%nonneg mean
if do_nneg
    Rvec.mean_nneg = Rvec_nneg/nb_dirs;
end


%auc
Rvec.auc    = sum(Rvec.mean);
%rtop
Rvec.rtop = Rvec.mean(1); % RTOP ( P(r)=0 )
% Save step size
Rvec.r_step_size_micron = r_step_size_micron;


%%%% Calculate derivatives of P(r)
%1st derivative mean
dRvec.mean = diff(Rvec.mean)/r_step_size_micron;
%1st deriv. auc
dRvec.auc = sum(dRvec.mean);
%1st deriv. "rtop"
dRvec.rtop = dRvec.mean(1);
%min of 1st deriv.
dRvec.min = min(dRvec.mean);
%std of 1st deriv.
dRvec.std = std(dRvec.mean);

%2nd derivative mean
ddRvec.mean = diff(dRvec.mean)/r_step_size_micron;
%2nd deriv. auc
ddRvec.auc = sum(ddRvec.mean);
%2nd deriv. "rtop"
ddRvec.rtop = ddRvec.mean(1);
%max of 2nd deriv.
ddRvec.max = max(ddRvec.mean);
%std of 2nd deriv.
ddRvec.std = std(ddRvec.mean);

%%%% Find zero-crossings
%1D EAP profile
zind = find(Rvec.mean<0); 
if ~isempty(zind)
    Rvec.zc = zind(1)*r_step_size_micron;
else
    Rvec.zc = 0;
end

%1st derivative
dzind = find(dRvec.mean>0); 
if ~isempty(dzind)
    dRvec.zc = dzind(1)*r_step_size_micron;
else
    dRvec.zc = 0;
end

%2nd derivative
ddzind = find(ddRvec.mean>0); 
if ~isempty(ddzind)
    ddRvec.zc = ddzind(1)*r_step_size_micron;
else
    ddRvec.zc = 0;
end

% Save step size
Rvec.r_step_size_micron = r_step_size_micron;
dRvec.r_step_size_micron = r_step_size_micron;
ddRvec.r_step_size_micron = r_step_size_micron;

% Save to output structure
out.Rvec = Rvec;
out.dRvec = dRvec;
out.ddRvec = ddRvec;
out.Rmatrix = Rmatrix;

end



