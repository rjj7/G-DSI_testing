function [ R ] = get_rot_mtx( ref, mov )
% 
% [ R ] = get_rot_mtx( mov, ref )
% 
% Returns rotation matrix R that rotates vector mov to vector ref.
%  
%  Ex.: 
%    Let a=[1;0;0] and b=[0.99;-0.13;0.06].
%     R = get_rot_mtx(b,a);
%     R = [ 0.9902   -0.1251    0.0628
%           0.1251    0.9921    0.0039
%           -0.0628    0.0039    0.9980];
%     R * mov = [1.0000; 0.0000; -0.0000] == ref
% 

%%%%%% EXAMPLE:
% ref = [1; 0; 0];
% mov = [0.9902; -0.1251; 0.0628];

if nargin<2
    disp('[ R ] = get_rot_mtx( mov, ref )');
    return
end

if size(ref,2)>size(ref,1), ref=ref'; end
if size(mov,2)>size(mov,1), mov=mov'; end

v = cross(mov, ref)';
% s = norm(v);
c = dot(mov, ref);
vx = [0,-v(3),v(2); v(3),0,-v(1); -v(2),v(1),0];
% vx2 = vx^2 * ((1-c)/s^2);
vx2 = vx^2 * (1/(1+c));

R = eye(3) + vx + vx2;

