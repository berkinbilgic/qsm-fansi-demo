# QSM Software Demo: 2016 Reconstruction Challenge with FANSI

ISMRM 2017 demo that replicates the phase processing pipeline of the QSM Reconstruction Challenge 2016 on 3D-GRE data acquired with 15-fold accelerated Wave-CAIPI, then compares dipole inversions: thresholded k-space division (TKD), closed-form L2 regularization, and nonlinear dipole inversion with TGV regularization using FANSI (Fast Nonlinear Susceptibility Inversion).

## Requirements: MEDI toolbox and STI Suite

Two third-party toolboxes used by the pipeline are **not included**. Download them and place them inside this folder; the script adds all subfolders to the path.

- STI Suite v1.42 (Laplacian phase unwrapping), by Wei Li and Chunlei Liu: https://people.eecs.berkeley.edu/~chunlei.liu/software.html
- MEDI toolbox (Laplacian boundary value background removal), Cornell MRI Research Lab: https://pre.weill.cornell.edu/mri/pages/qsm.html

`FANSI_v2/` is included. It was developed by Carlos Milovic in collaboration with Berkin Bilgic, Bo Zhao, Julio Acosta-Cabronero and Cristian Tejos; see `FANSI_v2/README.txt` for the publications to cite.

## Usage

Run `script_recon_challenge_v1.m` from the repository folder in MATLAB.

## Reference

B Bilgic. QSM Software Demo: 2016 Reconstruction Challenge. *ISMRM*, 2017. [[PDF]](https://www.martinos.org/~berkin/1645_Bilgic.pdf)

C Milovic, B Bilgic, B Zhao, J Acosta-Cabronero, C Tejos. Fast Nonlinear QSM Method with Variational Regularizations. *ISMRM*, 2017.

Contact: berkin AT nmr.mgh.harvard.edu
