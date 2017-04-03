%% Simple example on Saturn
%
% We present here the NLSPCA algorithm to denoise Poisson corrupted
% images.
% 
% Authors: \(\textbf{J. Salmon} \), \(\textbf{Z.  Harmany} \),
% \( \textbf{C-A. Deledalle} \) and \(\textbf{R. Willett} \) 
% Copyright (C) 2012 NLSPCA project.  
% See The GNU Public License (GPL)


%% Initialization
clear all
close all

addpath('functions')
addpath('tools')

%% Loading 'Saturn' image 
% ima_ori=double(imread('./data/saturn.tif'));
% ima_ori= ima_ori((1:256)+70,(1:256));
% Loading 'Saturn' image 
ima_ori=double(imread('./data/Ridges.png'));
%% Noisy image generation 
peak=0.1;
sd=2;
rng(sd)
Q = max(max(ima_ori)) /peak;
ima_lambda = ima_ori / Q;
ima_lambda(ima_lambda == 0) = min(min(ima_lambda(ima_lambda > 0)));
ima_nse_poiss = knuth_poissrnd(ima_lambda);
[m,n]=size(ima_nse_poiss);
%% Parameters:numbers

param.Patch_width=20; 
param.nb_axis=4; 
param.nb_clusters=14;
param.bandwith_smooth=2;
param.sub_factor=2;
param.big_cluster1=1;% special case for the biggest cluster 1st pass
param.big_cluster2=1;% special case for the biggest cluster 2nd pass                   
param.double_iteration=0;%1 or 2 pass of the whole algorithm
param.eps_stop=5e-3; %loop stoping criterion
param.epsilon_cond=1e-4;
param.nb_iterations=15;%number of iteration (Gordon)
param.cste=70;
param.func_tau=@(X) lasso_tau(X{1},X{2},param.cste);
param.bin=1; % subsampling factor, 0 is without it

%% Denoising with NLSPCA

tic
ima_fil=NLSPCA(ima_nse_poiss,ima_nse_poiss,param);
toc




%% Result display
figure('Position',[100 100   1200 600])
ax(1) = subplot(1, 3, 1);
plotimage(Q * ima_nse_poiss);
title(sprintf('Noisy PSNR = %f',psnr(Q*ima_nse_poiss, Q*ima_lambda, 255)));
ax(2) = subplot(1, 3, 2);
plotimage(Q * ima_lambda);
title('Original');
ax(3) = subplot(1, 3, 3);
plotimage(Q * ima_fil);
title(sprintf('NLSPCA, PSNR = %f',psnr(Q*ima_fil, Q*ima_lambda, 255)));
linkaxes(ax);

