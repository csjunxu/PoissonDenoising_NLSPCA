function label = litekmeans_poisson(X, k,label)
% litekmeans_poisson computes a clustering of the   noisy patches using a 
% Poisson model generating the data. Equivalent to a Kmeans for Poisson
% generated data. 
%                 
%   INPUT:
%     X                      : elements to cluster 
%     K                      : number of clusters used 
%     label                  : initialization of the clusters
%   OUTPUT:
%     label                  : indexes of the elements: each value refers 
%                              to the correponding cluster of the each point  
%
%
%   label = litekmeans_poisson(X, k,label)
%   Produce the clustering of the inpute in k clusters, iniatliazing the 
%   clusters with label
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


n = size(X,2);
last = 0;

if (nargin < 3)
    label = ceil(k*rand(1,n));  % random initialization    
end

while any(label ~= last)
    %[~,~,label] = unique(label);   % remove empty clusters
    E = sparse(1:n,label,1,n,k,n);  % transform label into indicator matrix
    center = X*(E*spdiags(1./sum(E,1)',0,k,k));    % compute center of each cluster
    last = label;
%    [~,label] = max(bsxfun(@minus,center'*X,0.5*sum(center.^2,1)')); % assign samples to the nearest centers
    [~,label] = max(bsxfun(@minus,log(center)'*X,sum(center,1)')); % assign samples to the nearest centers
end
