function ima_patchs = spatial_patchization_sub_sample(im_ini,w,bandwith_smooth,sub_factor)
% spatial_patchization computes the 3D  collection of patches from the
% noisy images. Patches have size Patch_width x Patch_width, but are 
% coded as vectors. A sub-sampling factor, using also a smoothing parameter,
% is provided.                 
%
%   INPUT:
%     im_ini                 : noisy image
%     Patch_width            : patch width 
%     bandwith_smooth        : bandwith for blurring before  subsampling 
%     sub_factor             : subsampling factor
%
%   OUTPUT:
%     ima_patchs             : 3D collection of patches (each patch is a 
%                              vector)
%    
%     ima_normalization      : "image" coding the number of patches a pixel
%                               pixel belong to.
%
%
%[ima_patchs, ima_normalization] = ...
%                                  spatial_patchization(im_ini,Patch_width)
%
% computes the 3D  collection of patches from the noisy images. Patches 
% have size Patch_width x Patch_width, but are coded as vectors. 
% A sub-sampling factor ,sub_factor, using also a smoothing parameter,
%  bandwith_smooth, is provided.          
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

[m,n] = size(im_ini);

w1=w-1;
m1=m-w1;
n1=n-w1;

ima_patchs=zeros(m1,n1,ceil(w/sub_factor)^2);
ima_filtered= conv2(im_ini,ones(bandwith_smooth)/bandwith_smooth^2,'same');
delta=(0:sub_factor:w1)-1;

for j = 1:n1
    
    yrange = mod(delta+j,n)+1;
    
    
    for i = 1:m1

        xrange = mod(delta+i,m)+1;  


        B = ima_filtered(xrange, yrange);
        ima_patchs(i,j,:) = B(:);

    end
end
