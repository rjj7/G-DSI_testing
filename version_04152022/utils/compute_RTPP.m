function [ rtpp ] = compute_RTPP( data, v, vperp_plane, pdd, qvec, rs, dispvec2 )

[ R ] = get_rot_mtx( v, pdd' );
rtpp_plane = R\vperp_plane; %inv(R) * vperp_plane;
[ rtpp_Fmtx ] = precompute_Fmtx( qvec, rtpp_plane', rs );
[ rtpp.Rmatrix ] = compute_1d_eap_mtx( data, rtpp_Fmtx );
rtpp.Rvec = mean(rtpp.Rmatrix,2);
rtpp.Rstd = std(rtpp.Rmatrix,0,2);
rtpp.RTPP = mean(rtpp.Rvec);
rtpp.RTPV = mean(rtpp.Rstd);
if nnz(size(dispvec2)~=size(rtpp.Rvec))>0
    dispvec2=dispvec2';
end
rtpp.msd  = mean(rtpp.Rvec.*dispvec2);
rtpp.kurtosis = kurtosis(rtpp.Rvec);


% ( data, v, pdd, qvec, rs, dispvec2 )
% 
% [ R ] = get_rot_mtx( v, pdd' );
% rtpp_axis = inv(R) * v;
% [ rtpp_Fmtx ] = precompute_Fmtx( qvec, rtpp_axis', rs );
% [ rtpp.Rvec ] = compute_1d_eap_mtx( data, rtpp_Fmtx );
% rtpp.RTPP = mean(rtpp.Rvec);
% rtpp.RTPV = std(rtpp.Rvec);
% if nnz(size(dispvec2)~=size(rtpp.Rvec))>0
%     dispvec2=dispvec2';
% end
% rtpp.msd  = mean(rtpp.Rvec.*dispvec2);
% rtpp.kurtosis = kurtosis(rtpp.Rvec);