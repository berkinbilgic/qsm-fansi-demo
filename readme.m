%%
% script_recon_challenge_v1.m : recon script that replicates the phase
% processing pipeline in QSM Recon Challenge 2016:
%
% http://www.neuroimaging.at/qsm2016/pages/welcome.php
%
%
% Data acquisition: 3D-GRE with 15-fold accelerated Wave-CAIPI acquisition
%   Rapid multi-orientation quantitative susceptibility mapping, B. Bilgic 
%   et al NIMG'16
%
%
% Brain masking: Brain Extraction Tool, FSL
%   Fast robust automated brain extraction, S. Smith HBM'02
%
% 
% Phase unwrapping: Laplacian Unwrapping, as implemented in STI Suite v1.42
% http://people.duke.edu/~cl160/
%   Quantitative susceptibility mapping of human brain reflects spatial 
%   variation in tissue composition, W. Li NIMG'11
%
%
% Background phase removal: Laplacian Boundary Value, as implemented in
% MEDI toolbox
% http://weill.cornell.edu/mri/pages/qsm.html
%   Background field removal by solving the Laplacian boundary value 
%   problem, D. Zhou et al NMR in Biomed'14  
%
%
% QSM recon 1: Thresholded K-space Division
%   Susceptibility mapping in the human brain using threshold?based k?space 
%   division, S. Wharton et al MRM'10
% and
%   Magnetic Susceptibility Mapping of Brain Tissue In VivoUsing MRI Phase
%   Data, K. Shmueli et al MRM'09
%
%
% QSM recon 2: Closed-form L2 regularization
%   Fast image reconstruction with L2 regularization, B. Bilgic et al JMRI'
%   14
%
%
% QSM recon 3: nonlinear dipole inversion with TGV regularization
%   A Fast Algorithm for Nonlinear QSM Reconstruction, C. Milovic et al
%   ISMRM'17 p3669


