function [lower,upper] = percentile_CI(values,sample_size,alpha)
% Construct percentile confidence intervals from a vector of estimates.
%
%   Inputs:
%   values - column vector of estimates of statistic determined using resampling such as delete-one jackknife
%   alpha - probability mass in tail of confidence intervals such that lower
%       CI = (alpha/2) and upper %CI = 1 - (alpha/2).
%
%   Outputs:
%   lower - lower CI = (alpha/2) percentile
%   upper - upper CI = 1-(alpha/2) percentile
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


    sorted_values = sort(values); % sort in ascending order

    %sample_size = size(values); % determine number of observations
    
    bottom_index = round((alpha/2)*sample_size);
    if bottom_index == 0
        bottom_index = 1;
    end

    top_index = round((1 - alpha/2)*sample_size);
    if top_index == sample_size
        top_index = sample_size - 1;
    end

    lower = sorted_values(bottom_index,1);
    upper = sorted_values(top_index,1);
