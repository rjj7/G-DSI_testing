function [ Rstats ] = fetch_metrics_Rvec( Rvec, opts )
% 
% [ dRvec, ddRvec ] = differentiate_Rvec( Rvec )
% 
% Differentiate vector of displacement probabilities of P(r). 
% Calculates 1st and 2nd derivatives Rvec.
% Rvec is 1D EAP profile, evaluated at nr distances along one EAP dir.
% Requires Rvec format to be as in output of compute_1d_gdsi_metrics().
% 
% INPUTS:
%   Rvec    = vec with 1D EAP displacement profile. 
% 
% OUTPUTS:
%   Rstats  = struct with stats
%      .min, .max, .mean, .std, .auc, .rtop
%
% RJ 2021

if nargin<1
    error('args:TooFew',...
        ' Must specify 1 input arguments! Run help for usage');
end

% optional "opts" 
if nargin<2
    %if not specified, compute all stats
    opts.doall = true;
    %otherwise, just use specified options
else
    if ~isfield(opts,'doall')
        opts.doall = false;
    end
end


%%%% Calculate stats
Rstats.opts = opts;

%mean
if opts.doall || opts.mean
    Rstats.mean = mean(Rvec);
end

%std dev
if opts.doall || opts.std
    Rstats.std = std(Rvec);
end

%sum (area under curve) (integrate)
if opts.doall || opts.auc
    Rstats.auc = sum(Rvec);
end

%rtop (first value, zero-displacement-prob)
if opts.doall || opts.rtop
    Rstats.rtop = Rvec(1);
end

%global minimum
if opts.doall || opts.min
    Rstats.min = min(Rvec);
end

%global maximum
if opts.doall || opts.max
    Rstats.max = max(Rvec);
end



end

    
    
    