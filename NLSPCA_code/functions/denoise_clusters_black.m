function ima_fil=denoise_clusters_black(ima_patchs_vect,func_denoising_patches,func_thresholding,...
func_recontruction,IDX,Patch_width,nb_axis,nb_clusters,m,n,normalization)

% denoise_clusters_black computes the denoising operation cluster by cluster
%                 
%   INPUT:
%     ima_patchs_vect        : collection of patches (matrix)
%     func_denoising_patches : choice of the denoising perform 
%                                (gausiann or poisson are handle)
%     Patch_width            : width of the square patches to consider
%     func_thresholding      : only no thresholding is yet supported
%     func_recontruction     : build the representaion of the patches
%                              thanks to the model used (gausiann or 
%                              poisson are handle)
%     IDX                    : indexes of the cluster
%     nb_clusters            : number of clusters choosen
%     nb_axis		     : number of axis/components choosen 
%     Patch_width            : width of the square patches to consider
%     m,n                    : size of the original image
%   OUTPUT:
%     ima_fil                : final denoised image
%     final_estimate         : collection of denoised patches
%   
%   ima_fil=denoise_clusters_black(ima_patchs_vect,func_denoising_patches,func_thresholding,...
%           func_recontruction,IDX,Patch_width,nb_axis,nb_clusters,m,n,normalization)
%   Perform the denoising of the collection of patches and the reprojection
%   to get an estimation of the original image (ima_fil).
%
%%   Copyright (C) 2012 NLSPCA project
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

final_estimate=zeros(size(ima_patchs_vect));
scores_patches=zeros(size(ima_patchs_vect,1),1);
startv = randn( nb_axis,Patch_width*Patch_width);


indexes=cell(nb_clusters,1);
ima_ppca_cluster=cell(nb_clusters,1);
ima_patchs_fil_cluster=cell(nb_clusters,1);

size_cluster=zeros(nb_clusters,1);
best=1;
max_number=0;

for k=1:nb_clusters
    
    indexes{k}=find(IDX==k);
    scores_patches(indexes{k})=length(indexes{k});
    size_cluster(k)=size(scores_patches(indexes{k}),1);
    if length(indexes{k})>max_number
        best=k;
        max_number=length(indexes{k});   
    end    
end
clear scores_patches

 value_mean=repmat(mean(ima_patchs_vect(indexes{best},:),2),[1,Patch_width*Patch_width]); 
 final_estimate(indexes{best},:)=value_mean;
%final_estimate(indexes{best},:)=mean(ima_patchs_vect(indexes{best},:),2)*ones(1,Patch_width*Patch_width);
index=1:nb_clusters;
index(best)=[];
        
for j=1:(nb_clusters-1)

        k=index(j);

%        fprintf('cluster size: %3d\n', size_cluster(k));

%         ima_ppca_cluster{k} = func_denoising_patches({ima_patchs_vect(indexes{k},:);zeros(size_cluster(k), nb_axis);startv});   
        ima_ppca_cluster{k} = func_denoising_patches({ima_patchs_vect(indexes{k},:);zeros(size_cluster(k), nb_axis);startv} );   
        ima_ppca_cluster{k} = func_thresholding(ima_ppca_cluster{k});           
        ima_patchs_fil_cluster{k} = func_recontruction(ima_ppca_cluster{k});            
        final_estimate(indexes{k},:) = ima_patchs_fil_cluster{k};
        ima_ppca_cluster{k} = [];
        ima_patchs_fil_cluster{k} = [];
                
end

ima_fil = reprojection_UWA(reshape(final_estimate,[(m-Patch_width+1),(n-Patch_width+1),Patch_width^2]));
ima_fil = ima_fil./normalization;



clear final_estimate

