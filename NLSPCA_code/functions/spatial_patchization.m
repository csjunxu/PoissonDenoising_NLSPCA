function ima_patchs = spatial_patchization(im_ini,w)
% spatial_patchization computes the 3D  collection of patches from the
% noisy images. Patches have size Patch_width x Patch_width, but are 
% coded as vectors.                 
%
%   INPUT:
%     im_ini                 : noisy image
%  	  Patch_width            : patch width 
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

ima_patchs=zeros(m1,n1,w^2);
delta=(0:w1)-1;
for j = 1:n1
    
    yrange = mod(delta+j,n)+1;
    
    for i = 1:m1

        xrange = mod(delta+i,m)+1;  

        
        B = im_ini(xrange, yrange);
        ima_patchs(i,j,:) = B(:);
    end
end
