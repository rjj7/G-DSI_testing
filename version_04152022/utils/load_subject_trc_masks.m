function [ trc ] = load_subject_trc_masks( trcdir )

trc.directory = trcdir;

tfile= [trcdir filesep 'fmajor.nii'];
[trc.fmajor]=load_trc_data(tfile);

tfile= [trcdir filesep 'fminor.nii'];
[trc.fminor]=load_trc_data(tfile);

tfile= [trcdir filesep 'lh.atr.nii'];
[trc.fmajor]=load_trc_data(tfile);

tfile= [trcdir filesep 'lh.cab.nii'];
[trc.cab.lh]=load_trc_data(tfile);
tfile= [trcdir filesep 'rh.cab.nii'];
[trc.cab.rh]=load_trc_data(tfile);
trc.cab.all = trc.cab.lh | trc.cab.rh;

tfile= [trcdir filesep 'lh.ccg.nii'];
[trc.ccg.lh]=load_trc_data(tfile);
tfile= [trcdir filesep 'rh.ccg.nii'];
[trc.ccg.rh]=load_trc_data(tfile);
trc.ccg.all = trc.ccg.lh | trc.ccg.rh;

tfile= [trcdir filesep 'lh.cst.nii'];
[trc.cst.lh]=load_trc_data(tfile);
tfile= [trcdir filesep 'rh.cst.nii'];
[trc.cst.rh]=load_trc_data(tfile);
trc.cst.all = trc.cst.lh | trc.cst.rh;

tfile= [trcdir filesep 'lh.ilf.nii'];
[trc.ilf.lh]=load_trc_data(tfile);
tfile= [trcdir filesep 'rh.ilf.nii'];
[trc.ilf.rh]=load_trc_data(tfile);
trc.ilf.all = trc.ilf.lh | trc.ilf.rh;

tfile= [trcdir filesep 'lh.slfp.nii'];
[trc.slfp.lh]=load_trc_data(tfile);
tfile= [trcdir filesep 'rh.slfp.nii'];
[trc.slfp.rh]=load_trc_data(tfile);
trc.slfp.all = trc.slfp.lh | trc.slfp.rh;

tfile= [trcdir filesep 'lh.slft.nii'];
[trc.slft.lh]=load_trc_data(tfile);
tfile= [trcdir filesep 'rh.slft.nii'];
[trc.slft.rh]=load_trc_data(tfile);
trc.slft.all = trc.slft.lh | trc.slft.rh;

tfile= [trcdir filesep 'lh.unc.nii'];
[trc.unc.lh]=load_trc_data(tfile);
tfile= [trcdir filesep 'rh.unc.nii'];
[trc.unc.rh]=load_trc_data(tfile);
trc.unc.all = trc.unc.lh | trc.unc.rh;

end


function [b]=load_trc_data(a)
    if exist(a,'file')
        b=load_untouch_nii(a);
        b=(b.img);
        b=b>0;
    else
        [i,j,k]=fileparts(a);
        if strcmpi(k,'.nii')
            newa = [i filesep j '.nii.gz'];
        elseif strcmp(k,'.gz')
            newa = [i filesep j];
        end
        if exist(newa,'file')
            b=load_untouch_nii(a);
            b=b.img;
            b=b>0;
        else
            disp('Could not find file!');
            disp(a);
            return
        end
    end
end
    

