function [ Rstats ] = fetch_metrics_Rtx( Rmatrix, opts )
% 
% [ dRvec, ddRvec ] = fetch_metrics_Rtx( Rmatrix, opts )
% 
% Get stats from Rmtx, which has EAP values at each distance for all
% directions on sphere. i.e. 2D matrix of size nb_r x nb_dirs
% 
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

%To calculate stats by averaging across directions
if opts.doall || opts.bydirs
    Rstats.bydirs.mean = mean(Rmtx,2);
    Rstats.bydirs.std = std(Rmtx,0,2);
end


%To calculate stats by averaging across distances
if opts.doall || opts.bydistances
    
    
end





%mean
if opts.doall || opts.mean
    Rstats.mean = mean(Rmatrix);
end

%std dev
if opts.doall || opts.std
    Rstats.std = std(Rmatrix);
end

%sum (area under curve) (integrate)
if opts.doall || opts.auc
    Rstats.auc = sum(Rmatrix);
end

%rtop (first value, zero-displacement-prob)
if opts.doall || opts.rtop
    Rstats.rtop = Rmatrix(1);
end

%global minimum
if opts.doall || opts.min
    Rstats.min = min(Rmatrix);
end

%global maximum
if opts.doall || opts.max
    Rstats.max = max(Rmatrix);
end




end

    
    
    