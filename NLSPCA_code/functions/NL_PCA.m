function [ima_fil,ima_fil_int,IDX_fil,IDX_int]=NL_PCA(ima_nse_poiss,...
		ima_nse_to_cluster,Patch_width,nb_axis,...
                nb_clusters,func_thresholding,func_recontruction,...
                func_denoising_patches,func_clustering,double_iteration,...
		bandwith_smooth,sub_factor,big_cluster1,big_cluster2)

% NL_PCA computes the denoising  of a noisy image. It is the main fucntion
% in the entire NLSPCA process.
%                 
%   INPUT:
%     ima_nse_poiss          : noisy image 
%     ima_nse_to_cluster     : image used fotr the clustering
%     Patch_width            : width of the square patches to consider
%     nb_axis		     : number of axis/components choosen 
%     nb_clusters            : number of clusters choosen
%     func_thresholding      : only no thresholding is yet supported
%     func_recontruction     : build the representaion of the patches
%                              thanks to the model used (gausiann or 
%                              poisson are handle)
%                              (current is the litekmeans)
%     func_denoising_patches : choice of the denoising perform 
%                                (gausiann or poisson are handle)
%     func_clustering	     : method to cluster the patches 
%     double_iteration       : either do the whole process one time (0)
%                              or use the denoised image to improve the 
%                              clustering (1)
%     bandwith_smooth        : smoothing parameter if subsampling is used 
%	             	       for clustering
%     sub_factor             : the subsampling factor used for the clustering 
%     big_cluster1           : special case for the biggest cluster at iter 1
%     big_cluster2           : special case for the biggest cluster at iter 2
%
%   OUTPUT:
%     ima_fil                : final denoised estimate 
%     ima_fil_int            : intermediate denoised estimate (when 
%                              double_iteration=0), otherwise = ima_fil
%     IDX_fil                : final clustering of the patches 
%     IDX_int                : intermediate clustering of the patches (when 
%                              double_iteration=0), otherwise = IDX_fil
%
% [ima_fil,ima_fil_int,IDX_fil,IDX_int]=NL_PCA(ima_nse_poiss,...
%		ima_nse_to_cluster,Patch_width,nb_axis,...
%               nb_clusters,func_thresholding,func_recontruction,...
%               func_denoising_patches,func_clustering,double_iteration,...
%		bandwith_smooth,sub_factor,big_cluster1,big_cluster2)
%
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
	
ima_patchs_sub_sampled = spatial_patchization_sub_sample(ima_nse_to_cluster, Patch_width,bandwith_smooth,sub_factor);
ima_patchs = spatial_patchization(ima_nse_poiss, Patch_width);
[~,~,w2]=size(ima_patchs);
ima_patchs_vect=reshape(ima_patchs,[(m-Patch_width+1)*(n-Patch_width+1),w2]); 
clear ima_patchs
label=ceil(nb_clusters*rand(1,(m-Patch_width+1)*(n-Patch_width+1)));
[~,IDX_int]=func_clustering({ima_patchs_sub_sampled,label});
nb_clusters=max(IDX_int(:));


normalization=normalization_UWA(m-Patch_width+1,n-Patch_width+1,Patch_width);


if big_cluster1==1    
    ima_fil_int=denoise_clusters_black(ima_patchs_vect,func_denoising_patches,func_thresholding,func_recontruction,IDX_int,Patch_width,nb_axis,nb_clusters,m,n,normalization);
else
    ima_fil_int=denoise_clusters(ima_patchs_vect,func_denoising_patches,func_thresholding,func_recontruction,IDX_int,Patch_width,nb_axis,nb_clusters,m,n,normalization);
end

label=ceil(nb_clusters*rand(1,(m-Patch_width+1)*(n-Patch_width+1)));

%Double iteration if wanted
if double_iteration==1
    ima_patchs_sub_sampled = spatial_patchization_sub_sample(ima_fil_int, Patch_width,bandwith_smooth,sub_factor);
    [~,IDX_fil]=func_clustering({ima_patchs_sub_sampled;label});
    nb_clusters=max(IDX_fil(:));
    clear ima_patchs_sub_sampled    
    if big_cluster2==1  
        ima_fil=denoise_clusters_black(ima_patchs_vect,func_denoising_patches,func_thresholding,func_recontruction,IDX_fil,Patch_width,nb_axis,nb_clusters,m,n,normalization);
    else        
        ima_fil=denoise_clusters(ima_patchs_vect,func_denoising_patches,func_thresholding,func_recontruction,IDX_fil,Patch_width,nb_axis,nb_clusters,m,n,normalization);
    end
else
    ima_fil=ima_fil_int;
    IDX_fil=IDX_int;    
end
clear ima_patchs_vect


