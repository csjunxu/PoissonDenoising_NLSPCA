function [ima_fil,ima_int,IDX_fil,IDX_int]=NLSPCA(...
                ima_nse_poiss,ima_nse_to_cluster,param)
% NLSPCA computes the denoising  of a noisy image. It is the main fucntion
% in the entiere NLSPCA process.
%                 
%   INPUT:
%     ima_nse_poiss          : noisy image 
%     ima_nse_to_cluster     : image used for the clustering part
%     param 	             : structure containing all the parameters 
%
%   OUTPUT:
%     ima_fil                : final denoised estimate 
%     ima_fil_int            : intermediate denoised estimate (when 
%                              double_iteration=0), otherwise = ima_fil
%     IDX_fil                : final clustering of the patches 
%     IDX_int                : intermediate clustering of the patches (when 
%                              double_iteration=0), otherwise = IDX_fil
%
% [ima_fil,ima_int,IDX_fil,IDX_int]=denoise_poisson_kmeans_poisson_PCA_l1(...
%                ima_nse_poiss,ima_nse_to_cluster,param)
%   Produce the denoised image ima_fil from the noisy input ima_nse_poiss
%   using the NLSPCA approach. Intermediate result is given when the
%   process is performed with two iterations. Moreover the indexes IDX_fil
%   ,IDX_int are provided mainly for visualiszing  the clustering 
%   performance.
%
%   Copyright (C) 2012 NLSPCA project
%   Joseph Salmon, Zachary Harmany, Charles-Alban Deledalle, Rebecca Willett 
%
%   See The GNU Public License (GPL)
%
%---------------------------------------------------------------------
%
%   This file is part of NLSPCA.
%
%   NLSPCA is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as
%   published by the Free Software Foundation, either version 3 of
%   the License, or (at your option) any later version.
%
%   NLSPCA is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public
%   License along with NLSPCA.  If not, see
%   <http://www.gnu.org/licenses/>.
%   Joseph Salmon, Zachary Harmany, Charles-Alban Deledalle, Rebecca Willett
%
%   See The GNU Public License (GPL)











[m,n]=size(ima_nse_poiss);
func_thresholding = @(ima_ppca) no_thresholding(ima_ppca);
func_clustering=@(X) clustering_poisson_litekmeans(X{1},param.Patch_width,param.nb_clusters,m,n,X{2});

%reconstruction used from the coefficient of the PCA to the patches
%could be reconstruction_gaussian_ppca for gaussian noise

func_recontruction=@(X)reconstruction_eppca(X);

%type of PCA calculation used: 
func_factorization=@(X)newton_SPIRAL(X{1},param.nb_axis,param.nb_iterations,X{2},X{3},param.eps_stop,param.epsilon_cond,param.func_tau);


func_denoising_patches=@(X)...
         poisson_NL_PCA(X{1},X{2},X{3},func_factorization);


     
binRadius=[param.bin,param.bin];
binSize=1+2*binRadius;
ima_nse_binned=conv2(ima_nse_poiss,ones(binSize),'same');
ima_nse_to_cluster_binned=conv2(ima_nse_to_cluster,ones(binSize),'same');


samples1=binRadius(1)+1:binSize(1):size(ima_nse_poiss,1)-binRadius(1);
samples2=binRadius(2)+1:binSize(2):size(ima_nse_poiss,2)-binRadius(2);

ind_sample1=1:size(ima_nse_poiss,1);
ind_sample2=1:size(ima_nse_poiss,2);

ima_nse_downsized=ima_nse_binned(samples1,samples2);
ima_nse_to_cluster_downsized=ima_nse_to_cluster_binned(samples1,samples2);





[ima_fil1,ima_int,IDX_fil,IDX_int]=NL_PCA(ima_nse_downsized,ima_nse_to_cluster_downsized,param.Patch_width,param.nb_axis,...
                          param.nb_clusters,func_thresholding,func_recontruction,...
                         func_denoising_patches,func_clustering,...
                         param.double_iteration,param.bandwith_smooth,param.sub_factor,...
                         param.big_cluster1,param.big_cluster2);

                     
ima_fil=ima_fil1/prod(binSize);

ima_fil=interp1(samples1,ima_fil,ind_sample1,'linear','extrap');
ima_fil=interp1(samples2,ima_fil',ind_sample2,'linear','extrap')';

                     
                     



