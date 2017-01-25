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
ima_ori=double(imread('./data/saturn.tif'));
GT_Original_image_dir = 'C:\Users\csjunxu\Desktop\CVPR2017\cc_Results\Real_ccnoise_denoised_part\';
GT_fpath = fullfile(GT_Original_image_dir, '*mean.png');
TT_Original_image_dir = 'C:\Users\csjunxu\Desktop\CVPR2017\cc_Results\Real_ccnoise_denoised_part\';
TT_fpath = fullfile(TT_Original_image_dir, '*real.png');
GT_im_dir  = dir(GT_fpath);
TT_im_dir  = dir(TT_fpath);
im_num = length(TT_im_dir);

% ima_ori= ima_ori((1:256)+70,(1:256));
% Loading 'Saturn' image
% ima_ori=double(imread('./data/Ridges.png'));
% %% Noisy image generation
peak=120; % 120 1 5 10 30 60
% rand('seed',0)
% Q = max(max(ima_ori)) /peak;
% ima_lambda = ima_ori / Q;
% ima_lambda(ima_lambda == 0) = min(min(ima_lambda(ima_lambda > 0)));
% ima_nse_poiss = knuth_poissrnd(ima_lambda);
% [m,n]=size(ima_nse_poiss);
%% Parameters:numbers
method = 'NLPCA';
%% write image directory
write_sRGB_dir = ['C:/Users/csjunxu/Desktop/CVPR2017/cc_Results/' method '/'];
if ~isdir(write_sRGB_dir)
    mkdir(write_sRGB_dir)
end

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

PSNR = [];
SSIM = [];
CCPSNR = [];
CCSSIM = [];
for i = 1 : im_num
    IMin = im2double(imread(fullfile(TT_Original_image_dir,TT_im_dir(i).name) ));
    IM_GT = im2double(imread(fullfile(GT_Original_image_dir, GT_im_dir(i).name)));
    S = regexp(TT_im_dir(i).name, '\.', 'split');
    IMname = S{1};
    [h,w,ch] = size(IMin);
    fprintf('%s: \n',TT_im_dir(i).name);
    CCPSNR = [CCPSNR csnr( IMin*255,IM_GT*255, 0, 0 )];
    CCSSIM = [CCSSIM cal_ssim( IMin*255, IM_GT*255, 0, 0 )];
    fprintf('The initial PSNR = %2.4f, SSIM = %2.4f. \n', CCPSNR(end), CCSSIM(end));
    IMout = zeros(size(IMin));
    for cc = 1:ch
        %% denoising
        ima_fil=NLSPCA(IMin(:,:,cc),IMin(:,:,cc),param);
        IMout(:,:,cc) = ima_fil;
    end
    %% output
    PSNR = [PSNR csnr( IMout*255, IM_GT*255, 0, 0 )];
    SSIM = [SSIM cal_ssim( IMout*255, IM_GT*255, 0, 0 )];
    fprintf('The final PSNR = %2.4f, SSIM = %2.4f. \n', PSNR(end), SSIM(end));
    %% output
    imwrite(IMout, [write_sRGB_dir method '_RID_' IMname '.png']);
end
mPSNR = mean(PSNR);
mSSIM = mean(SSIM);
mCCPSNR = mean(CCPSNR);
mCCSSIM = mean(CCSSIM);
save(['C:/Users/csjunxu/Desktop/CVPR2017/cc_Results/BM3DPoisson_RID' num2str(im_num) '.mat'],'PSNR','mPSNR','SSIM','mSSIM','CCPSNR','mCCPSNR','CCSSIM','mCCSSIM');


% %% Denoising with NLSPCA
%
% tic
% ima_fil=NLSPCA(ima_nse_poiss,ima_nse_poiss,param);
% toc




% %% Result display
% figure('Position',[100 100   1200 600])
% ax(1) = subplot(1, 3, 1);
% plotimage(Q * ima_nse_poiss);
% title(sprintf('Noisy PSNR = %f',csnr(Q*ima_nse_poiss, Q*ima_lambda, 0,0)));
% ax(2) = subplot(1, 3, 2);
% plotimage(Q * ima_lambda);
% title('Original');
% ax(3) = subplot(1, 3, 3);
% plotimage(Q * ima_fil);
% title(sprintf('NLSPCA, PSNR = %f',csnr(Q*ima_fil, Q*ima_lambda, 0,0)));
% linkaxes(ax);

