function [ ACQ ] = read_acqspecs_file( fspecs )

if ~exist(fspecs,'file')
    error('fileSpecification:doesNotExist','\nSpecs file does not exist');
end

fid = readtable(fspecs);
fid = table2cell(fid);

ACQ = [];
ACQ.units = [];

for row=1:length(fid)
    param = fid{row,1};
    value = fid{row,2};
    switch param
        case {'gmax','Gmax'}
            ACQ.Gmax = value/1000;
            ACQ.units.Gmax = 'T/m';
        case {'Delta','bigdelta'}
            ACQ.big_delta = value/1000;
            ACQ.units.big_delta = 's';
        case {'delta','smalldelta'}
            ACQ.small_delta = value/1000;
            ACQ.units.small_delta = 's';
        case 'qmax'
            ACQ.qmax = value/1000;
            ACQ.units.qmax = 'um^-1';
        case 'deltaq'
            ACQ.deltaq = value/1000;
            ACQ.units.deltaq = 'um^-1';
        case 'bmax'
            ACQ.bmax = value;
            ACQ.units.bmax = 'mm^2/s';
        case 'qradius'
            ACQ.qradius = value;
            ACQ.units.qradius = 'radius of q-space Cartesian lattice';
        case {'te','TE'}
            ACQ.TE = value;
            ACQ.units.TE = 'ms';
        case {'tr','TR'}
            ACQ.TR = value;
            ACQ.units.TR = 'ms';
        case 'echospacing'
            ACQ.echo_spacing = value;
            ACQ.units.echo_spacing = 'ms';
        otherwise
            fprintf(' -Could not parse config option %s...\n',param);
    end  
end

end