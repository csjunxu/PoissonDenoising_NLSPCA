function [x] = SPIRALTAPBregman(y,A,Aprime,xinit,tau,miniter,maxiter)
% SPIRALTAPBregman reconstructs a sparse image or vector from photon-limited
% (Poisson distributed) measurements.
%                 
%   INPUT:
%     y                         : Noisy photon count data
%     A          		: Forward system operator
%     Aprime 			: Adjoint system operator 
%     xinit 			: Initialization
%     tau             		: l1 regularization parameter
%     miniter 			: minimum number of iterations to perform
%     maxiter 			: maximum number of iterations to perform
%
%   OUTPUT:
%     x 			: Reconstructed image
% [x] = SPIRALTAPBregman(y,A,Aprime,xinit,tau,miniter,maxiter)
%
% This function solves the optimization problem
%
%   minimize - log p( y | A x ) + tau ||x||_1                           
%      x
%   subject to x >= 0                                                    
%
% where p( y | A x ) is the Poisson likelihood, tau a regularization parameter.
% If tau is a vector (of the same size as x), then this function will also solve 
%                                                                        
%   minimize - log p( y | A x ) + || diag(tau) x||_1 
%      x
%   subject to x >= 0                                                    
%
% where diag(tau) is a diagonal matrix with entries specified by tau. For more
% algorithmic details see our journal paper:
%
% Z. T. Harmany, R. F. Marcia, and R. M. Willett, “This is SPIRAL-TAP: Sparse
% Poisson Intensity Reconstruction ALgorithms—Theory and Practice,” IEEE 
% Transactions on Image Processing, vol. 21, pp. 1084–1096, Mar. 2012.
%
%   Copyright (C) 2012 NLSPCA project, 2009-2012 SPIRALTAP project
%   Joseph Salmon, Zachary T. Harmany, Charles-Alban Deledalle,
%   Roummel F. Marcia, Rebecca M. Willett
%   Corresponding author: Zachary T. Harmany (zth@duke.edu)                   
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
   

%====================================
%= Set default / initial parameters =
%====================================
% Problem specification options
%AT = [];
%simplebound = 1;
%lowerbound = 0;
truth = [];
%init = 0;
% Treatment of alpha (step size) parameter
alpha = 1;
% alphamin = 1e-30;
% alphamax = 1e+30;

alphamin = 1e-30;
alphamax = 1e+30;


%alphamethod = 1;
monotone = 1;
acceptalphamax = 1e+30;
acceptmult = 2;
acceptdecrease = 0.1;
acceptpast = 10;
% Algorithm termination options
%stopcriterion = 3;
tolerance = (1e-6)^2;%%modification no more sqrt

% Output options
verbose = 0;
saveobjective = 0;
savereconerror = 0;
savecputime = 0;
reconerrortype = 0;

%saveiteratepath = 0; % Can be memory intensive, don't save by default
%savesteppath = 0; 
% Advanced options
%warnings = 1;
% Set non-optional initial parameters
converged = 0;
iter = 0;
% Reserve variable names
termval = [];
objective = [];
reconerror = [];
cputime = [];
alphapath = [];


  saveobjective = 1;

termval = zeros(maxiter+1,1);

x = xinit;
Ax = A*x;
xprevious = x;
Axprevious = Ax;
dx = x - xprevious;
grad = Aprime*(exp(Ax) - y); 


[termval(iter+1), ~] = checkTermination;
 if saveobjective;   objective(iter+1) = computeObjective;    end

%=============================
%= Display beginning message =
%=============================
if (verbose > 0); thetime = fix(clock);
  for ii = 1:60; fprintf('='); end; fprintf('\n');
  fprintf([ '= Beginning SPIRAL Canonical            ',...
            '@ %2d:%02d %02d/%02d/%4d =\n',...
            '=   Miniter:  %-5d               ',...
            'Maxiter:  %-5d          =\n'],...
    thetime(4),thetime(5),thetime(2),thetime(3),thetime(1),miniter,maxiter);  
  for ii = 1:60; fprintf('='); end; fprintf('\n');
end
iter = iter+1;


%=============================
%= Begin main algorithm loop =
%=============================
while (iter <= miniter) || ((iter <= maxiter) && ~converged)
  
      if monotone % Check acceptance criterion
        past = (max(iter-1-acceptpast,0):iter-1) + 1;
        maxpastobjective = max(objective(past));
        accept = 0;
        while (accept == 0)
          % Compute a candidate next iterate
          dx = xprevious;
          step = xprevious - grad./alpha;
%jo
%           x = computeSubproblemSolution;
          quotient=-tau./alpha; 
          x = max( 0, step+ quotient) - max( 0, -step + quotient); 
          dx = x - dx;
          %Adx = Axprevious;
          Ax = A*x;
          Adx = Ax - Axprevious;
          normsqdx = sum(dx(:).^2);
          
          % Compute the resulting objective 
          objective(iter+1) = computeObjective;
          if ( objective(iter+1) <= (maxpastobjective ...
                - acceptdecrease*alpha/2*normsqdx) ) ...
                || (alpha >= acceptalphamax);
            accept = 1;
          end
          %alphaaccept = alpha;  % Keep value for displaying
          alpha = acceptmult*alpha;
        end
%         if savealphapath; alphapath(iter+1) = alphaaccept; end
      else
        % Just take the Barzilai-Borwein step, no enforcing monotonicity.
        dx = xprevious;
        step = xprevious - grad./alpha;
%jo
%         x = computeSubproblemSolution;
        quotient=-tau./alpha; 
        x = max( 0, step+ quotient) - max( 0, -step + quotient); 
        dx = x - dx;
        Adx = Axprevious;
        Ax = A*x;
        Adx = Ax - Adx;
        normsqdx = sum(dx(:).^2);
%         if saveobjective; objective(iter+1) = computeObjective; end
%         if savealphapath; alphapath(iter+1) = alpha; end
      end
      
      
      gamma = Adx'*(Adx.*exp(Ax));
%jo       
%        if gamma == 0
%          alpha = alphamin;
%        else
        alpha = min(alphamax, max(gamma/normsqdx, alphamin));
%        end
      
%   end
  
  %=== Compute items for next iteration ===
  xprevious = x;
  Axprevious = Ax; 
  diff_expAx_y=exp(Ax) - y;
  grad = Aprime*(diff_expAx_y); 

  %=== Check convergence and calculate output quantities ===
  [termval(iter+1) converged] = checkTermination;
  % Note: Objective is saved above (since it is needed when monotone = 1)
%   if savereconerror;  reconerror(iter+1) = computeReconError; end
%   if savecputime;     cputime(iter+1) = toc;                  end
%   if saveiteratepath; iteratepath{(iter+1)} = x;              end
%   if savesteppath;    steppath{iter+1} = step;                end
  
  %=== Display Progress ===
  %jo
  %displayProgress;
    
	iter = iter + 1;
end


%==============================
%= Display completion message =
%==============================
if (verbose > 0); thetime = fix(clock);
  for ii = 1:60; fprintf('='); end; fprintf('\n');
  fprintf([ '= Completed SPIRAL Canonical           ',...
            ' @ %2d:%02d %02d/%02d/%4d =\n',...
            '=   Iterations performed:  %-5d ',...
            '                          =\n'],...
    thetime(4),thetime(5),thetime(2),thetime(3),thetime(1),iter-1);  
  for ii = 1:60; fprintf('='); end; fprintf('\n');
end


%===============================================================================
%= Helper Subfunctions =========================================================
%===============================================================================

%=========================
%= Objective Computation =
%=========================
function objective = computeObjective
    precompute = exp(Ax)-y.*Ax;
    objective = sum(precompute)+ tau*sum(abs(x));
end



%====================================
%= Reconstruction Error Computation =
%====================================
function reconerror = computeReconError
  errorvect = x - truth;
  switch reconerrortype
    % - Based on squared error -
    case 0 % Squared error
      reconerror = sum(errorvect(:).^2);
    case 1 % Squared error per pixel
      reconerror = sum(errorvect(:).^2)./numel(errorvect);
    case 2 % Relative squared error
      reconerror = sum(errorvect(:).^2)./sum(truth(:).^2);
    case 3 % Relative squared error, as a percent
      reconerror = 100*sum(errorvect(:).^2)./sum(truth(:).^2);
      
    % - Based on l2 norm -
    case 4 % Root squared error (l2 norm)
      reconerror = sqrt(sum(errorvect(:).^2));
    case 5 % Relative root squared error (l2 norm)
      reconerror = sqrt(sum(errorvect(:).^2))./sqrt(sum(truth(:).^2));
    case 6 % Relative root squared error (l2 norm), as a percent
      reconerror = 100*sqrt(sum(errorvect(:).^2))./sqrt(sum(truth(:).^2));
    
    % - Based on l1 norm -
    case 7 % Absolute error (l1 norm)
      reconerror = sum(abs(errorvect(:)));
    case 8 % Absolute error (l1 norm) per pixel
      reconerror = sum(abs(errorvect(:)))./numel(errorvect);
    case 9 % Relative absolute error (l1 norm)
      reconerror = sum(abs(errorvect(:)))./sum(abs(truth(:)));
    case 10 % Relative absolute error (l1 norm), as a percent
      reconerror = 100*sum(abs(errorvect(:)))./sum(abs(truth(:)));
      
    % - Based on PSNR -
    case 11 % PSNR, using maximum true intensity
      reconerror = sum(errorvect(:).^2)./numel(errorvect);
      reconerror = max(truth(:)).^2./reconerror;
      reconerror = 10*log10(reconerror);
    case 12 % PSNR, using dynamic range
      reconerror = sum(errorvect(:).^2)./numel(errorvect);
      reconerror = (max(truth(:)) - min(truth(:))).^2./reconerror;
      reconerror = 10*log10(reconerror);
  end
end

%====================================
%= Termination Criteria Computation =
%====================================
function [termval converged] = checkTermination
      if iter == 0
        termval = NaN;
      else
        termval = (sum(dx(:).^2)./sum(x(:).^2));
      end
      converged = (termval <= tolerance);
    
    
end

%====================
%= Display Progress =
%====================
function displayProgress
  if ~mod(iter,verbose)
    fprintf('%3d|',iter)
    if savecputime
      fprintf(' t:%3.2fs',cputime(iter+1))
    end
    fprintf(' dx:%10.4e', sqrt(sum(dx(:))));
    if savealphapath
      fprintf(' alpha used:%10.4e',alphapath(iter+1));
    end
    fprintf(' alpha BB:%10.4e',alpha);
    if saveobjective
      fprintf(' obj:%11.4e',objective(iter+1));
      fprintf(' dobj:%11.4e',objective(iter+1) - objective(iter));
    end
    if savereconerror
      fprintf(' err:%10.4e', reconerror(iter+1));
    end
    fprintf(' term:%11.4e (target %10.4e)', termval(iter+1),tolerance);
    fprintf('\n');
  end
end
    
        
  
end
