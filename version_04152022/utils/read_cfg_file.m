function [ DATAIN ] = read_cfg_file( fconfig )

if ~exist(fconfig,'file')
    error('fileSpecification:doesNotExist','\nConfig file does not exist');
end

fid = importdata(fconfig);

DATAIN = [];

for row=1:length(fid)
    data = strsplit(fid{row});
    param = data{1};
    value = data{2};
    switch param
        case 'fdwi'
            DATAIN.fdwi = value;
        case 'fmask'
            DATAIN.fmask = value;
        case 'fbvec'
            DATAIN.fbvec = value;
        case 'fbval'
            DATAIN.fbval = value;
        case 'outdir'
            DATAIN.outdir = value;
        case 'acqspecs'
            DATAIN.fspecs = value;
        case 'procdir'
            DATAIN.procdir = value;
        case 'aparcaseg'
            DATAIN.faparcaseg = value;
        case 'fsdir'
            DATAIN.FSdir = value;
        case 'isshelled'
            DATAIN.shelledflag = value;
        otherwise
            fprintf(' -Could not parse config option %s...\n',param);
    end  
end

if ~isfield(DATAIN,'fdwi')
    error('fieldSpecMissing:dwi','\nNo fdwi specified\n');
end
if ~isfield(DATAIN,'fmask')
    error('fieldSpecMissing:mask','\nNo fmask specified\n');
end
if ~isfield(DATAIN,'fbvec')
    error('fieldSpecMissing:bvec','\nNo fbvec specified\n');
end
if ~isfield(DATAIN,'fbval')
    error('fieldSpecMissing:bval','\nNo fbval specified\n');
end
if ~isfield(DATAIN,'outdir')
    error('fieldSpecMissing:outdir','\nNo outdir specified\n');
end

end