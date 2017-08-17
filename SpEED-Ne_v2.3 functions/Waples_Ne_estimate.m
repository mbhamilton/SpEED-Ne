function [W_Ne_estimate] = Waples_Ne_estimate(sample_size,r_squared_delta)

% Function to estimate Ne from r^2 based on Waples(2006) post hoc
% regression fits. Equations based on Table 1 in Waples and Do (2007). Includes
% multiplication by N/(N-1) to give unbiased estimate of Ne from r_squared delta.
%
% This function works with scalar or vectorized inputs. Vector computation
% requires that both sample_size and r_squared inputs have identical row by single column dimensions.
%
% version 1.0  21 July 2017
%
% Inputs:
% sample_size - sample size of individuals for r^2
% r_squared_delta - r^2 delta estimate
%
% Outputs:
% W_Ne_estimate - Ne estimated using Waples and Do (2007) equations
%
%
%**************
%   Copyright 2017 by Matthew B Hamilton. 
%
%   This file is part of SpEED-Ne: Simulation & Estimation of gamEtic Disequilibrium Ne.
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
%   A copy of the GNU General Public License is available at http://www.gnu.org/licenses/.
%**************

    % check to verify that har_mean_S_delta has an equal number of rows as r_squared_delta_AFW,r_squared_delta_AFT,r_squared_delta_AFW_UB and r_squared_delta_AFT_UB
    % for vectorized computation
    if isequal(size(sample_size), size(r_squared_delta)) ~= 1
        fprintf('Error: sample_size and r_squared_delta are unequal sizes!');
        return;
    end
    
    index = sample_size < 30; % get logical index, 1 if less than 30 zero otherwise
    
    W_Ne_estimate = index .* (0.308 + sqrt(0.308^2 + 2.08.*(r_squared_delta - (0.0018 + 0.907./sample_size + 4.44./(sample_size.^2))).^2) ) ./ (2.*(r_squared_delta - (0.0018 + 0.907./sample_size + 4.44./(sample_size.^2)))) +...
        ~index .* (1/3 + sqrt(1/9 + 2.76 .* (r_squared_delta - (1./sample_size + 3.19./(sample_size.^2))).^2) )./(2 .* (r_squared_delta - (1./sample_size + 3.19./(sample_size.^2))));

    
%     if sample_size < 30
%         E_r2_delta_sample = 0.0018 + 0.907./sample_size + 4.44./(sample_size.^2);
% 
%         r_2_prime = r_squared_delta - E_r2_delta_sample;
%         W_Ne_estimate = (0.308 + sqrt(0.308^2 + 2.08.*r_2_prime.^2) )./(2.*r_2_prime); % Waples Ne from r_squared_delta
% 
%     else % sample_size >= 30
%         E_r2_delta_sample = 1./sample_size + 3.19./(sample_size.^2);
% 
%         r_2_prime = r_squared_delta - E_r2_delta_sample;
%         W_Ne_estimate = (1/3 + sqrt(1/9 + 2.76.*r_2_prime.^2) )./(2.*r_2_prime); % Waples Ne from r_squared_delta AFW
% 
%     end
%     
%     disp(W_Ne_estimate_vec)
%     disp(W_Ne_estimate)
    
    
        
    
    
    
    
    
    
    
    
    
    
