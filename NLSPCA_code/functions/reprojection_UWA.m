function ima = reprojection_UWA(ima_patchs)
%   reprojection_UWA compute the reprojection step that from a collections 
%   of patches representing an image, provide an image estimate. The 
%   various pixels estimates (due to overlapping) are uniformly average.
%                 
%   INPUT:
%     ima_patchs             : denoised collections of patches 
%     xg,yg                  : potential subpart of the image
%
%   OUTPUT:
%    ima                     : denoised images by uniformly averaging the
%                              patches containing a pixel of interest to
%                              proposed a denoised pixel estimator
%
%   ima = reprojection_UWA(ima_patchs,xg,yg)
%   
%   Compute the reprojection step that from a collections 
%   of patches ima_patchs representing an image, provide an image estimate.
%   The various pixels estimates (due to overlapping) are uniformly
%   average.  Potentially restrict the process to a smaller window,
%   centered around xg,yg
%   
%   Copyright (C) 2012 NLSPCA project
%   Joseph Salmon,  Zachary Harmany, Charles-Alban Deledalle, Rebecca Willett
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
%   Joseph Salmon,  Zachary Harmany, Charles-Alban Deledalle, Rebecca Willett
%
%   See The GNU Public License (GPL)

    [M, N, P] = size(ima_patchs);

    p = sqrt(P);
    p1=p-1;
    M1=M+p1;
    N1=N+p1;
    ima = zeros(M+p-1, N+p-1);

    delta=(-1)+(1:p)-1;
    
for j = 1:N
    
    yrange = mod(delta+j,N1)+1;
        
    for i = 1:M

        xrange = mod(delta+i,M1)+1;        

        ima(xrange, yrange) = ima(xrange, yrange) + ...
                                reshape((ima_patchs(i,j,:)), p, p);

    end
end

