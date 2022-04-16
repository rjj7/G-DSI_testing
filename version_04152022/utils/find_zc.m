function [ zc ] = find_zc( Rvec )
% 
% [ Rvec ] = find_zc( Rvec )
% 
% Find (first) zero crossing of 1D displacement profile. 
% 
% INPUTS:
%   Rvec    = 1D vector to analyze.
% 
% OUTPUTS:
%   zc      = index of first zero crossing 
%             [really index of first point whose sign is -1*sign(Rvec(1))]
% 
% RJ 2021

%1D EAP profile
zind = find(Rvec<0); 
if ~isempty(zind)
    zc = zind(1);
else
    zc = NaN;
end

end