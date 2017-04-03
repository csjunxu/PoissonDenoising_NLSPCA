function ima_ppca = no_thresholding(ima_ppca)
% denoise_clusters computes the denoising operation cluster by cluster
%                 
%   INPUT:
%     ima_ppca        : collection of patches (matrix)

%   OUTPUT:
%     ima_fil                : final denoised image
%   
%  ima_ppca = no_thresholding(ima_ppca)
%   Performs pretty much... nothing  on the PCA coeffecients. More is
%   coming
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
%   Joseph Salmon,Zachary Harmany, Charles-Alban Deledalle, Rebecca Willett
%
%   See The GNU Public License (GPL)

ima_ppca.coefs = ima_ppca.coefs;
