function [ dRvec, ddRvec ] = differentiate_Rvec( Rvec, stepsize )
% 
% [ dRvec, ddRvec ] = differentiate_Rvec( Rvec )
% 
% Differentiate vector of displacement probabilities of P(r). 
% Calculates 1st and 2nd derivatives Rvec.
% Rvec is 1D EAP profile, evaluated at nr distances along one EAP dir.
% Requires Rvec format to be as in output of compute_1d_gdsi_metrics().
% 
% INPUTS:
%   Rvec    = struct with 1D EAP displacement profile, calculated from
%             compute_1d_gdsi_metrics(). 
% 
% OUTPUTS:
%   dRvec   = struct with 1st derivative of Rvec metrics
%   ddRvec  = struct with 2nd derivatives of Rvec metrics
%
% RJ 2021

if nargin<1
    error('args:TooFew',...
        ' Must specify 1 input arguments! Run help for usage');
end


%%%% Calculate derivatives of P(r)
%1st derivative mean
dRvec = diff(Rvec)/stepsize;

%2nd derivative mean
ddRvec = diff(dRvec)/stepsize;

    
    
end