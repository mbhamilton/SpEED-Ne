function [lower,upper] = normal_dist_CI(estimate,values,sample_size,alpha)
% Construct normal distribution 95% confidence intervals from a vector of delete-one jackknife estimates assuming a normal distribution.
%
%   http://people.bu.edu/aimcinto/jackknife.pdf
%
%   Inputs:
%   estimate - point estimate of r^2
%   values - column vector of estimates of statistic determined using resampling such as delete-one jackknife
%   alpha - probability mass in tail of confidence intervals such that lower
%       CI = (alpha/2) and upper %CI = 1 - (alpha/2).
%
%   Outputs:
%   lower - lower CI = estimate - z(alpha/2)*SE(jack)
%   upper - upper CI = estimate + z(1-alpha/2)*SE(jack)
%
%**************
%   Copyright 2017 Matthew B Hamilton.
%   
%   This file is part of SpEED-Ne: Simulation & Estimation of gEnotypic Disequilibrium Ne.
%
%   SpEED-Ne is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
% 
%   SpEED-Ne is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
% 
%   A copy of the GNU General Public License is available to http://www.gnu.org/licenses/.
%**************


    % compute deviate value for std normal based on alpha
    mu = 0;
    sigma = 1;
    pd = makedist('Normal',mu,sigma); % generate standard normal distribution with mean zero and unit variance
    constant = icdf(pd,1-alpha/2); % get standard normal deviate
    
    average = mean(values);
    sqr_diffs = (values - average).^2;
    
    var_j = 1/sample_size * sum(sqr_diffs); % variance of jackknife pseudovalues
    
    se_j = sqrt(var_j); % standard error of jackknife pseudovalues
    
    lower = estimate - constant*se_j;
    
    upper = estimate + constant*se_j;
