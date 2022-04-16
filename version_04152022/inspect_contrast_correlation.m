
voldir = '/autofs/space/nyx_002/users/rjones/Bay8/CS_DSI_10_01_2021/preproc/nufft/4sh/gdsi/nii';

flist = {'RTOP','qvar','qmsd','qiv','qentropy','NLrms','NG','MeanKurt',...
    'MeanGA','KLDi','KLDg','HWHM','HWHMstd','HWHMmin','HWHMmax',...
    'EntropyEAP','T1w','T2w'};

for fl=1:length(flist)
    fname = flist{fl};
    fprintf(' - on %s \n',fname);
    vol = load_untouch_nii([voldir filesep fname '.nii.gz']);
    eval([fname ' = double(vol.img);']);
end

fprintf(' - on brain_mask \n');
fmask = '/autofs/space/nyx_002/users/rjones/Bay8/CS_DSI_10_01_2021/preproc/brain_mask.nii.gz';
vol = load_untouch_nii(fmask);
mask = double(vol.img);


ContrastCorrMtx = zeros(length(flist));

for fl=1:length(flist)
    fname = flist{fl};
    eval(['curr = ' fname ';']);
    for ff=1:length(flist)
        fname2 = flist{ff};
        eval(['mov = ' fname2 ';']);
        ContrastCorrMtx(fl,ff) = corr(curr(mask==1),mov(mask==1));
    end
end
