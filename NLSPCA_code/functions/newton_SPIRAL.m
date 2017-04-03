function [u, v, lastexpy] = newton_SPIRAL(x, l, iters, startu, startv,epsilon_stop,epsilon_cond,func_tau)

% newp1svd_SPIRAL computes the Newton algorithm to factorize the matrix of 
% of the noisy patches in the Poisson domain, using an ell_1 regularization
% handle adapting the SPIRAL algorithm.
%                 
%   INPUT:
%     x                      : noisy patches 
%     l                      : number of axis/components choosen 
%     iters                  : number of iterations to perform the Newton
%                              step
%     startu                 : initialization of the coefficient in the
%                              matrix factorization / pca representation
%     startv                 : initialization of the axis/components in the
%                              matrix factorization / pca representation
%     eps_stop               : relative stopping  criterion for the Newton
%                              step
%     epsilon_cond           : small conditioning number to be sure the
%                              Hessian is well conditionned when doing 
%                              the Newton step
%     func_tau               : function for coding the tuning parameter for
%                              ell_1 constraint
%
%   OUTPUT:
%     u                      : coefficient of the patches in the PCA-domain
%    
%     v                      : axis/component/dictionary element to 
%                              represent the patches in the PCA-domain
%
%     lastexpy               : representation of the patches using only the
%                              nb_axis elements
%
% [u, v, lastexpy] = newton_SPIRAL(x, l, iters, startu, ...
%	       	                   startv,epsilon_stop,epsilon_cond,func_tau)
%
%   Compute the Newton algorithm to factorize the matrix of noisy patches,
%   using nb_axis, a maximum amount of nb_iterations for the main loop and
%   a stoping criterion eps_stop, for controlling the relative change. 
%   startu and starv are initial guess for the factorization, and 
%   func_tau is the trade-off on the ell_1 constraint.
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

[n, d] = size(x);

tau=func_tau({d,n});

u = startu;
v = startv;

norma=sqrt(sum(v(1,:).^2,2));
v(1,:)=norma*ones(1,d)/sqrt(d);


cond_mat=epsilon_cond*eye(l);
product_uv=u*v;
lastexpy=exp(product_uv/2);
curnorm = sum(sum(lastexpy.^2)) ;

L_tot=curnorm-sum(sum(x.*product_uv));


xinit = zeros(l,1);
miniter = 2;
maxiter = 20;


for iter = 1:iters
            
    vprime=v';
    for i = 1:n               
        u(i,:) = SPIRALTAPBregman(x(i,:)',vprime,v,xinit,tau,miniter,maxiter);
    end

  lastexpy2 = exp(u*v/2);
  lasterru= u'*(x-lastexpy2.^2);
  
  for j = 1:d    

     int_v=bsxfun(@times,u,lastexpy2(:,j));
     covar = (int_v'*int_v);       
     covar=covar+cond_mat;         
     v(:,j) =v(:,j)+ covar \ (lasterru(:,j));
    
  end
 
    
    product_uv=u*v;
    lastexpy2=exp(product_uv/2);


    lastexp_square=lastexpy2.^2;
    curnorm = sum(sum(lastexp_square)) ;


    L_tot2=curnorm-sum(sum(x.*product_uv));
 
    ratio=abs((L_tot-L_tot2)/L_tot);
    %verbose or not
%     fprintf('%3d:  %g \n', iter,ratio);

    if ( (ratio < epsilon_stop))
    break;
    end

    L_tot=(L_tot2);
    lastexpy=lastexpy2;


end
