function patchs_ppca = poisson_NL_PCA(patchs,startu,startv,func_factorization)

% poisson_NL_PCA computes the denoising  of the noisy patches using a 
% Poisson model generating the data. Equivalent to a PCA decomposition. 
%                 
%   INPUT:
%     patchs                 : noisy patches 
%     nb_axis                : number of axis/components choosen 
%     nb_iterations          : number of iterations to perform the Newton
%                              step
%     startu                 : initialization of the coefficient in the
%                              matrix factorization / pca representation
%     startv                 : initialization of the axis/components in the
%                              matrix factorization / pca representation
%     func_factorization     : factorization used for representating the data
%                              (PCA, Poisson PCA, etc.) default is newton_SPIRAL
%   OUTPUT:
%     patchs_ppca            : structure representing the dictionnary/pca 
%                              representaion (coefficients, axis, etc...)
%
%
% patchs_ppca = poisson_NL_PCA(patchs,nb_axis,nb_iterations,...
%                        startu,startv,func_factorization)
%   Produce the structure of PCA representing all the patches through a
%   matrix factorization/PCA decomposition adapted to the data generated 
%   with a Poisson distribution.
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


    if size(patchs, 3) > 1
        isimage = 1;
        [M,N,P] = size(patchs);
        MN = M*N;
        patchs = reshape(patchs, MN, P);
    else
        isimage = 0;
        [MN,P] = size(patchs);
    end
    
     [u v ] =func_factorization({patchs,startu,startv});        
    

    patchs_ppca.coefs = u;
    patchs_ppca.dicos{1}.axis = v;

