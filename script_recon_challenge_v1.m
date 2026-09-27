%%-------------------------------------------------------------------------
% QSM DEMO
% showcasing QSM Recon Challenge 2016
% data available at
%                       _____________________
%                       ||    qsm.rocks    ||
%
% includes input phase, ground truth QSM, evaluation metrics, simple dipole
% inversion methods for benchmarking
% 
% goals:
%   - test ability of existing algorithms to recover susceptibility from in
%   vivo phase data
%   - serve as a common dataset for testing and benchmarking future
%   algorithms
% 
% 
% Challenge toolbox did not include the phase processing pipeline, only raw
% and tissue phase images
% 
%
% in this demo:
%   - bridge this gap by starting from raw, wrapped phase and go through 
%   all the processing steps that went into the Challenge 
%   - a simple, stand-alone toolbox that could be a starting point for QSM
%   recon in clinical and research studies
% 
% software available at
%                       ________________________________
%                       ||    martinos.org/~berkin    ||
%%-------------------------------------------------------------------------


%%-------------------------------------------------------------------------
%% load data
% complex image
% header info
% brain mask from BET
% COSMOS reference
%%-------------------------------------------------------------------------

set(0,'DefaultFigureWindowStyle','docked')

addpath(genpath(pwd))


load Img;           % complex data
load prot;          % header information
load Mask_bet;      % brain mask from BET
load chi_cosmos;    % COSMOS ground truth QSM


spatial_res = [255/240, 255/240, 180/168];
N = size(Img);

TE = prot.alTE * 1e-6;           % sec
B0 = prot.flNominalB0;           % Tesla
gyro = 2*pi*42.58;
 

imagesc3d2(Img, N/2, 1, [90,90,-90], [0,5e-2], [], 'Magnitude'), set(gcf, 'color', [1,1,1]*0.15)
imagesc3d2(Img .* Mask_bet, N/2, 2, [90,90,-90], [0,5e-2], [], 'Magnitude * Bet Mask'), set(gcf, 'color', [1,1,1]*0.15)



%%-------------------------------------------------------------------------
%% Laplacian unwrapping: 
% STI Suite (v1.42)
%%-------------------------------------------------------------------------


Phase_wrapped = angle(Img) .* Mask_bet;


tic
    Phase_unwrapped = LaplacianPhaseUnwrap(Phase_wrapped, spatial_res, [0,0,0]) .* Mask_bet;
toc
    

imagesc3d2(Phase_wrapped - 10 * (~Mask_bet), N/2, 3, [90,90,-90], [-pi,pi], [], 'Wrapped Phase')
set(gcf, 'color', [1,1,1]*0.15), colorbar, set(gca, 'fontsize', 32)

imagesc3d2(Phase_unwrapped - 10 * (~Mask_bet), N/2, 4, [90,90,-90], [-pi,pi]*2, [], 'Unwrapped Phase')
set(gcf, 'color', [1,1,1]*0.15), colorbar, set(gca, 'fontsize', 32)



%%-------------------------------------------------------------------------
%% remove background phase: 
% Laplacian Boundary Value (LBV) in MEDI toolbox
%%-------------------------------------------------------------------------


tol = 1e-2;         % default
depth = -1;         % default
peel = 5;           % no of voxels to peel from the mask boundary
     
 
tic
    Phase_lbv = LBV(Phase_unwrapped, Mask_bet, N, spatial_res, tol, depth, peel);
toc    


% scale by TE, B0, and gyromagnetic ratio to convert to ppm
Nfm_lbv = Phase_lbv / (TE*B0*gyro);

% tissue mask from LBV
Mask_lbv = Phase_lbv~=0;

% estimated background field map
Nfm_back = (Phase_unwrapped - Phase_lbv) .* Mask_lbv / (TE*B0*gyro);


imagesc3d2(Nfm_lbv - 10 * (~Mask_lbv), N/2, 5, [90,90,-90], [-0.05,0.05], [], 'LBV: normalized field map')
set(gcf, 'color', [1,1,1]*0.15), colorbar, set(gca, 'fontsize', 32)

imagesc3d2(Nfm_back - 10 * (~Mask_lbv), N/2, 6, [90,90,-90], [-0.2,0.2], [], 'background field map')
set(gcf, 'color', [1,1,1]*0.15), colorbar, set(gca, 'fontsize', 32)


%%-------------------------------------------------------------------------
%% 3d polynomial fit to remove transmit phase from LBV
% Roemer coil combination fails to remove B1+ phase and includes some
% contribution from B1- of the body coil 
%%-------------------------------------------------------------------------

Order = 4;      % degree of 3d polyfit


tic
    I_fitted = polyfit3D_NthOrder(Nfm_lbv, Mask_lbv, Order);
toc


% subtract fitted phase
Nfm_final = Nfm_lbv - I_fitted .* Mask_lbv;


imagesc3d2(Nfm_final - 10 * (~Mask_lbv), N/2, 7, [90,90,-90], [-0.05,0.05], [], 'Remove Tx phase: normalized field map')
set(gcf, 'color', [1,1,1]*0.15), colorbar, set(gca, 'fontsize', 32)


imagesc3d2(I_fitted - 10 * (~Mask_lbv), N/2, 8, [90,90,-90], [-0.05,0.05], [], 'Removed component')
set(gcf, 'color', [1,1,1]*0.15), colorbar, set(gca, 'fontsize', 32)



%%-------------------------------------------------------------------------
%% create dipole kernel
%%-------------------------------------------------------------------------

[ky,kx,kz] = meshgrid(-N(1)/2:N(1)/2-1, -N(2)/2:N(2)/2-1, -N(3)/2:N(3)/2-1);

kx = (kx / max(abs(kx(:)))) / spatial_res(1);
ky = (ky / max(abs(ky(:)))) / spatial_res(2);
kz = (kz / max(abs(kz(:)))) / spatial_res(3);

k2 = kx.^2 + ky.^2 + kz.^2;


% transverse plane: identity rotation
Rotation_Mtx = eye(3);      


kernel = fftshift( 1/3 - (kx * Rotation_Mtx(3,1) + ky * Rotation_Mtx(3,2) + kz * Rotation_Mtx(3,3)).^2 ./ (k2 + eps) );    

    
kernel_disp = fftshift(abs(kernel));
figure(9), imagesc(squeeze(kernel_disp(1+end/2,:,:)), [0,2/3]), axis image off, colormap gray, title('dipole kernel', 'fontsize', 32, 'color', 'w')
set(gcf, 'color', [1,1,1]*0.2), colorbar, set(gca, 'fontsize', 32)

 

%%-------------------------------------------------------------------------
%% QSM: TKD recon
% some noise amplification and dipole artifacts present
% avoids spatial smoothing
%%-------------------------------------------------------------------------


thre_tkd = 0.17;        % TKD threshold parameter


% compute inverse kernel by truncating small values in the inversion
kernel_inv = zeros(N);
kernel_inv( abs(kernel) > thre_tkd ) = 1 ./ kernel(abs(kernel) > thre_tkd);


chi_tkd = real( ifftn( fftn(Nfm_final) .* kernel_inv ) ) .* Mask_lbv; 

rmse_tkd = 100 * norm2( Mask_lbv .* (chi_tkd - chi_cosmos) ) / norm2( Mask_lbv .* chi_cosmos )


imagesc3d2(chi_tkd - 10 * (~Mask_lbv), N/2, 10, [90,90,-90], [-0.10,0.14], [], ['QSM: TKD   rmse=', num2str(rmse_tkd), '%'])
set(gcf, 'color', [1,1,1]*0.2), colorbar, set(gca, 'fontsize', 32)



%%-------------------------------------------------------------------------
%% QSM: closed-form L2 recon
% dipole inversion with smooth image gradients
% dipole artifacts mitigated at the cost of over smoothing
%%-------------------------------------------------------------------------


[k1, k2, k3] = ndgrid(0:N(1)-1,0:N(2)-1,0:N(3)-1);

E1 = 1 - exp(2i .* pi .* k1 / N(1));
E2 = 1 - exp(2i .* pi .* k2 / N(2));
E3 = 1 - exp(2i .* pi .* k3 / N(3));

EtE = abs(E1).^2 + abs(E2).^2 + abs(E3).^2;

DtD = abs(kernel).^2;


% gradient regularization parameter
reg_param = 6e-2;       

chi_L2 = real( ifftn(conj(kernel) .* fftn(Nfm_final) ./ (DtD + reg_param * EtE))) .* Mask_lbv;

rmse_L2 = 100 * norm2( Mask_lbv .* (chi_L2 - chi_cosmos) ) / norm2( Mask_lbv .* chi_cosmos )


imagesc3d2(chi_L2 - 10 * (~Mask_lbv), N/2, 11, [90,90,-90], [-0.10,0.14], [], ['QSM: CF-L2   rmse=', num2str(rmse_L2), '%'])
set(gcf, 'color', [1,1,1]*0.2), colorbar, set(gca, 'fontsize', 32)

  

%%-------------------------------------------------------------------------
%% QSM: FANSI   
% fast nonlinear solver using total generalized variation 
% C. Milovic et al #3669
% Monday @ 16:15
% computer 103
%   fast version of nonlinear-MEDI recon using variable splitting
%   good trade-off between dipole artifact mitigation and smoothing 
%   uses magnitude image for noise weighting
%   two parameters that need to be tuned
%%-------------------------------------------------------------------------
 

magn_use = abs(Img) .* Mask_lbv;
magn_use = magn_use / max(abs(magn_use(:)));


params = [];
params.K = kernel;
params.N = N;
params.input = Nfm_final;
params.weight = magn_use; 

params.maxOuterIter = 40;

params.mu1 = 5e-5;          % gradient consistency
params.alpha1 = 3e-5;       % gradient L1 penalty (TGV constraint)



outnl = nlTGV(params); 


rmse_fansi = 100 * norm2( Mask_lbv .* (chi_cosmos - outnl.x) ) / norm2( Mask_lbv .* chi_cosmos )


imagesc3d2(outnl.x - 10 * (~Mask_lbv), N/2, 12, [90,90,-90], [-0.10,0.14], [], ['QSM: FANSI #3669   rmse=', num2str(rmse_fansi), '%'])
set(gcf, 'color', [1,1,1]*0.2), colorbar, set(gca, 'fontsize', 32)


imagesc3d2(chi_cosmos - 10 * (~Mask_lbv), N/2, 13, [90,90,-90], [-0.10,0.14], [], 'QSM: COSMOS')
set(gcf, 'color', [1,1,1]*0.2), colorbar, set(gca, 'fontsize', 32)


 
