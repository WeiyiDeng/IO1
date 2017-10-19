% GMMEST: Main procedure for the unrestricted GMM estimation 
%
% SYNTAX 
% [theta_final, J_test, probJ, S_final, moments, moments_grad, bandw,...
% var_theta, std_theta, conf_inter] = 
% gmmest(options, data, 'popmom', stval, Wstart, varargin);
% 
% INPUTS
%   options   : Options structure for fmincon (use OPTIMSET)
%               Set the options structure for fminunc, the unconstrained
%               minimization function used by MATLAB. To see the list of
%               possible options, enter optimset('fminunc').          
%   data      : A Txm matrix with the dataset that the moments are based on 
%               (m variables with T observations each). 		
%   popmom    : An Matlab function that calculates the moment conditions and their 
%               gradient. The function must be of the form 
%               [mom, gradmom] = popmom(theta,data,varargin)
%               where mom is the moments and gradmom their gradient.                 
%   stval     : A vector with an initial guess for the parameters. 	
%   Wstart    : The weighting matrix for the first step estimation. 	
%   [varargin]: Additional parameters passed to the popmom function
%               For example, in the case of a GIV estimation, 
%               you can pass the vector of the instruments.
%
% OUTPUT
%   theta_final    : A vector with the final GMM estimates.
%   J_test         : The J-test of overidentifying restrictions.
%   probJ          : The probability of the J-test.
%   S_final          The long-run covariance matrix of the moments, 
%                    evaluated at the final GMM estimates.
%   moments        : A matrix with the values of the moment conditions,
%                    evaluated at the final GMM estimates
%   moments_grad   : The gradient of the moment conditions, evaluated at 
%                    the final GMM estimates
%          bandw   : The bandwidth used to calculate S_final
%   var_theta      : The variance-covariance matrix of the estimates.
%   std_theta      : The standard error of each estimates, defined as 
%                    [Var(estimate)/T]^0.5.
%   conf_inter     : The 95% confidence interval of each estimate.
% -------------------------------------------------------------------------
% 
% SETABLE OPTIONS (use OPTSET):
%   center : A dummy variable, taking the values 0 or 1
%            Set to 0: The variance of the moments is calculated using the 
%                    uncentered moment conditions               
%            Set to 1: The variance of the moments is calculated using the 
%                    centered moment conditions      
%            ##The default value is 0## 
%   method : Covariance matrix estimation method. Set this option equal to            
%           'HACC_B', for HAC with Bartlett kernel
%           'HACC_P', for HAC with Parzen kernel1 
%           'SerUnc', for Serially uncorrelated 
%           ##The default method is 'SerUnc'## 
%   bandw  : The bandwidth used in THE HAC estimation. The bandwidth must 
%           be a non-negative integer.    
%           ## If the bandwidth is not given by the user, and a HACC estimator 
%           has been selected by the user, the "optimal" bandwidth is automatically 
%           calculated using Newey and Wests's Method of Bandwidth Selection
%   itergmm : Maximum number for the iterated GMM estimator.
%           ##(The default is 50).
%   tol     : Tolerance criterion for the iteration procedure.
%           ##(The default is 1e-006).


%<--------------THE GMM ESTIMATION PROCEDURE STARTS HERE-------------->
function [theta_final, J_test, probJ, S_final, final_moments,... 
          final_moments_grad, bandw, var_theta, std_theta,conf_inter] = gmmest(options, data, popmom, startval, We, varargin);

% ERROR CHECK
if nargin<5, error('The first four inputs (data, popmom, stval, W) must be provided by the user');end
[stvr,stvc] = size(startval);
if stvc~=1, error('The starting values must be a column vector');end
pmc=feval(popmom, startval,data, varargin{:});
[~,q]=size(pmc);
if stvr>q, error('The system is under-identified. You must supply at least as many moment conditions as parameters.');end

% OPTIONS STRUCTURE FOR GMM (DEFAULT VALUES) 
center    = optget('gmmest','center',0);
method    = optget('gmmest','method','SerUnc');
bandw     = optget('gmmest','bandw',[]);
itergmm   = optget('gmmest','itergmm',50);
tol       = optget('gmmest','tol',1e-006);

% First step estimator
theta = fminunc('gobj', startval, options, popmom, data, We, varargin{:}); 

% Iterative estimation starts here
for i=2:itergmm
    pmc = feval(popmom, theta,data, varargin{:}); % Calculate the pmc and their gradient
    S = longvar(pmc, center, method, bandw);      % Calculate the covariance matrix of the moments
    invS = S\eye(size(S,1)); % Inverse of S, computed with Gaussian elimination, ...
 %   options = optimset(options,'Display','off');
    thetanew = fminunc('gobj', theta,...
        options, popmom, data, invS, varargin{:});   
    if norm(abs(theta - thetanew)) < tol
        result = sprintf('The algorithm converged to a solution. The optimal estimator was achieved in iteration %2.0f .', i);
        disp(result);
        break
    end
    theta = thetanew;
end
if exist('result','var')==0 & itergmm~=1
    disp('The algorithm didn''t converged to a solution.' );
end

if itergmm == 1
    thetanew = theta;
    pmc = feval(popmom, thetanew, data, varargin{:});
    S = longvar(pmc, center, method, bandw);
    disp('One step GMM estimation: Completed');
    disp('J-test has a non-standard distribution; p-value is not calculated');
end

theta_final = thetanew;
if nargout>1
    [pmc,dpmc] = feval(popmom, theta_final,data, varargin{:}); 
    W_final = S\eye(size(S,1)); % Inverse of S, computed with Gaussian elimination, ...
	[S_final, bandw] = longvar(pmc, center, method, bandw); 
	[N, pc] = size(pmc); 
	df = pc - stvr;
	final_moments = pmc;
	final_moments_grad = dpmc;
    if itergmm==1 
        J_test = gobj(theta_final, popmom, data, We, varargin{:});
        [VAR,SD,CI] = varest(dpmc, S_final, theta_final, N, We,1);
        probJ = [];
    else
        [VAR,SD,CI] = varest(dpmc, S_final, theta_final, N,We,itergmm);
        J_test = gobj(theta_final, popmom, data, inv(S_final), varargin{:});
        probJ = 1-chi2cdf(J_test, df);
    end
	var_theta = VAR;
	std_theta = SD;
	conf_inter = CI;
end

% USER NOTIFICATIONS
if isempty(optget('gmmest', 'bandw')) & lower(optget('gmmest', 'method'))~='serunc'
    message = sprintf('The optimum bandwidth, has been set to %4.0f', bandw);
    disp(message);
end