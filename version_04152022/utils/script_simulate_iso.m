function [S] = script_simulate_iso(S0,b, grad)
% function [S, D] = create_signal_multi_tensor (ang, f, Eigenvalues, b, grad, S0, SNR, add_noise)
% -Normalizing the gradient vector and then transforming the b-value.
% -This part is only for non-unitary gradient vectors
% Transform_b_value_and_grad = repmat(sqrt(diag(grad*grad')+eps), [1 3]);
% grad = grad./Transform_b_value_and_grad;
% b = b.*(Transform_b_value_and_grad(:,1)).^2;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % S0        : signal at b=0
% % SNR       : snr
% % lambda1   : diffusivity, fiber 1
% % lambda2   : diffusivity, fiber 2
% % f1        : volume fraction, fiber 1
% % f2        : volume fraction, fiber 2 
% % anglesF   : 
% %         [ang1 ang2]: -> ang1 is the rotation in degrees from x to y axis (rotation using the z-axis);
% %                         ang2 is the rotation in degrees from the x-y plane to the z axis (rotation using the y axis)
% %         anglesF(i,:) defines the orientation of the fiber with volume fraction fi. 
% 
% adc_range = [1e-4,1e-3];
% 
% % Gmax = 480e-3; %480e-3; %T/m
% % %adc = 1.7e-4; %mm^2/s
% % adc =  9.65e-5; %5e-5; %mm^2/s
% % fw_adc = 2.5e-3;  %mm^2/s
% % smdel = 15e-3; %s
% % bigdel = 21.116035e-3; %19e-3; %s
% % gmr = 42.57747892e6;  %T^-1*s^-1
% % 
% % mdd = sqrt(6*adc*(bigdel-(smdel/3)))*1000;  %um
% 
% % --- Experimental parameters
% % diffusivities
% lambda1 = 1*10^-4;% mm^2/s
% lambda2 = 0.1*10^-4;% mm^2/s
% % S0 is the image at b-value = 0
% 
% S0 = 1;
% % Signal to noise ratio
% SNR = 30;
% display(['The signal to noise ratio in this experiment is SNR = ' num2str(SNR) '.']);
% 
% % --- Intravoxel properties
% % Using a mixture of two fibers (more fibers can be used)
% % volume fraction for each fiber (the sum of this elements is 1 by definition)
% f1 = 0.5; f2 = 0.5;
% % angular orientation corresponding to each fiber (main eigenvector for each diffusion tensor)
% anglesF = [0  0;
%            90 0];
% % [ang1 ang2]: -> ang1 is the rotation in degrees from x to y axis (rotation using the z-axis);
% %                 ang2 is the rotation in degrees from the x-y plane to the z axis (rotation using the y axis)
% % anglesF(i,:) defines the orientation of the fiber with volume fraction fi. 
% 
% %  --- Computation of the diffusion signal using the above experimental parameters and intravoxel properties
% % *** Signal with realistic rician noise
% add_rician_noise = 1;
% S = create_signal_multi_tensor(anglesF, [f1  f2],...
%     [lambda1, lambda2, lambda2], bvalues, grad, S0, SNR, add_rician_noise);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f = [1]; %[1/3 1/3 1/3];
ang = [0,0]; %[0,0; 0,90; 90,90];
adc = 2.5e-3;
Eigenvalues = [adc adc adc];
% b = DATAIN.bvals;
% grad = DATAIN.bvecs;
% S0 = 1;
SNR = 1;
add_noise = 0;


A = diag(Eigenvalues);

S = 0;
Nfibers = length(f);
f = f/sum(f);

for i = 1:Nfibers
    phi(i) = ang(i, 1);
    theta(i) = ang(i, 2);
    R = RotMatrix(phi(i),theta(i));
    D = R*A*R';
    S = S + f(i)*exp(-b.*diag(grad*D*grad'));
end
S = S0*S;

sigma = S0/SNR;

standar_deviation = sigma.*(ones(length(grad),1));
med = zeros(length(grad),1);

er1 = normrnd(med, standar_deviation);
er2 = normrnd(med, standar_deviation);

if add_noise == 1
    S = sqrt((S + er1).^2 + er2.^2);
end
disp(' ');

return

% --- private funtions -----------
function R = RotMatrix(phi,theta)

c = pi/180;
phi = phi*c;
theta = theta*c;

Rz = [ cos(phi)  -sin(phi)  0
       sin(phi)   cos(phi)  0
           0         0      1];


Ry = [cos(theta)   0   sin(theta)
          0        1         0
     -sin(theta)   0   cos(theta)];

R =  Rz*Ry;
return