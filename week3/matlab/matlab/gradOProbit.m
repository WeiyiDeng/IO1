function [ll]=gradOProbit(par,x,y)

%% GRADOPROBIT calculates the loglikelihood for an ordered probit model
%
% [grad] = GRADOPROBIT(par,x,y) computes the length(y) x length(par) matrix
% grad of score contributions of an ordered probit model with parameters
%
%       par(1) - level of first threshold
%       par(2:end-1) - distances between subsequent thresholds
%       par(end) - standard deviation of probit error term
%
% for data
%
%       x - length(y) x 1 vector of covariate values
%       y - length(y) x 1 vector of outcomes in {1,2,...,length(par)}

thresholds = [-inf cumsum(par(1:end-1)) inf]';
sigma = par(end);
p = normcdf((x-thresholds(y))/sigma) - normcdf((x-thresholds(y+1))/sigma);
ll = sum(log(p));

end