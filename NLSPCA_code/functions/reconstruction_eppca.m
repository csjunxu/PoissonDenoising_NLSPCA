function patchs = reconstruction_eppca(patchs_pca)

% reconstruction_eppca computes from the coefficients in the PCA domain the 
% representation of the patches in the original patch-space, when using a
% Poisson model for the data generation process.
%                 
%   INPUT:
%    patchs_pca      : structure containing the patches coefficients in
%                      the patch domain (patchs_pca.axis) and the
%                      dictionary elements/compenents (patchs_pca.dicos)
%                              
%   OUTPUT:
%     patchs         : collection of patches in the original
%                      patch-space. 
%
% patchs = reconstruction_eppca(patchs_pca)
%
%   Produce the structure the collection of denoising patching starting 
%   from the information coming from the PCA (Poisson) information.
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
coefs = patchs_pca.coefs;
[nb_axis,P]=size(patchs_pca.dicos{1}.axis);

    if size(coefs, 3) > 1
        isimage = 1;
        [M,N,R] = size(coefs);
        MN = M*N;
        coefs=reshape(coefs,[M*N, R]);

    else
        if nb_axis==1
        
            [M,N]=size(coefs);   
            MN=M*N;
            coefs=reshape(coefs,[M*N, 1]);
        else
            isimage = 0;
            [MN,R] = size(coefs);        

    
        end
    end

    L = length(patchs_pca.dicos);
    patchs = zeros(MN, P);

    
    for i = 1:L
        

         patchs =exp(coefs*patchs_pca.dicos{i}.axis);

         
    end
   
    clear coefs;


    if nb_axis==1
        
        patchs = reshape(patchs,M,N,P);

    else
        if isimage
            patchs = reshape(patchs, M, N, P);


        end
    end
