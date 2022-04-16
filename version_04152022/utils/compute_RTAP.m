function [ rtap ] = compute_RTAP( data, v, pdd, qvec, rs, dispvec2 )

[ R ] = get_rot_mtx( v, pdd' );
rtap_axis = R\v; %inv(R) * v;
[ rtap_Fmtx ] = precompute_Fmtx( qvec, rtap_axis', rs );
[ rtap.Rvec ] = compute_1d_eap_mtx( data, rtap_Fmtx );
rtap.RTAP = mean(rtap.Rvec);
rtap.RTAV = std(rtap.Rvec);
if nnz(size(dispvec2)~=size(rtap.Rvec))>0
    dispvec2=dispvec2';
end
rtap.msd  = mean(rtap.Rvec.*dispvec2);
rtap.kurtosis = kurtosis(rtap.Rvec);


% 
