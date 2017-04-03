function [ima_patchs_vect,IDX]=clustering_poisson_litekmeans(ima_patchs,Patch_width,nb_clusters,m,n,label)
% clustering_poisson_litekmeans computes the clusters of patches extracted 
%                 
%   INPUT:
%     ima_patchs          : collection of patches 
%     Patch_width         : width of the square patches to consider
%     nb_clusters	  : number of clusters choosen
%     m,n                 : size of the original image
%     label               : initialization vector 
%   OUTPUT:
%     ima_patchs_vect     : reshape the 3D collection of patches to 2D
%     IDX                 : indexes of the patches clusters
%   
% [ima_patchs_vect,IDX]=clustering_poisson_litekmeans(ima_patchs,Patch_width,nb_clusters,m,n,label)
%   tranform the whole collections of  3D patches of size 
%   Patch_width x Patch_width to a 2D (matrix). The indexes give the
%   membership of each patch to a cluster (among nb_clusters choices)
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
%   Joseph Salmon, Zachary Harmany, Charles-Alban Deledalle, Rebecca Willet
%
%   See The GNU Public License (GPL)

[m2,n2,w2]=size(ima_patchs);
ima_patchs_vect=reshape(ima_patchs,[(m2)*(n2),w2]); 
IDX = litekmeans_poisson(ima_patchs_vect',nb_clusters,label);

