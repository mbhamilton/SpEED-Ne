function [jack_pct_CI_values,jack_norm_CI_values,chi_CI_values,effective_n_c,effective_n_delta] = wd_confidence_intervals_v2(num_loci,r_squared_c_values,r_squared_delta_values,c_locus_pairs,delta_locus_pairs,...
    r_squared_c_W,r_squared_c_TH,r_squared_delta_W,r_squared_delta_TH,correction_factor_har_mean_S,median_r_sqr_c_AFW_perm,median_r_sqr_c_AFT_perm,alpha,print_output,make_graphs)

% Function to compute confidence intervales for estimates of r^2 and Ne
% based on methods decribed by Waples (2006) and Waples & Do (2008).
%
% version 2.0 20 Sept 2016
%
% Inputs:
%   num_loci - integer half the number of columns in the genotype data matrix
%   r_squared_c_values - c_pairs by 3 vector of r_squared_c for each locus pair. 
%       1st column is numerator for AFW estimate, 2nd column is numerator for AFT estimate, 3rd column is denominator. 
%
%   r_squared_delta_values - delta_pairs by 3 vector of r_squared_delta for each locus pair. 
%       1st column is numerator for AFW estimate, 2nd column is numerator for AFT estimate, 3rd column is denominator.
%
%   c_locus_pairs - sample size of pairs of loci for r_squared_c estimates
%   delta_locus_pairs  - sample size of pairs of loci for estimates of r_squared_delta
%   r_squared_c_W - multilocus estimate of r^2_{comp} with AFW
%   r_squared_c_TH - multilocus estimate of r^2_{comp} with AFT
%   r_squared_delta_W - multilocus estimate of r^2_{delta} with AFW
%   r_squared_delta_TH - multilocus estimate of r^2_{delta} with AFT   
%   correction_factor_har_mean_S - real number E(r^2) based on finite sample size of
%       individuals, harmonic mean S over locus pairs
%   alpha - probability mass in tail of confidence intervals such that lower
%       CI = (alpha/2) and upper %CI = 1 - (alpha/2).
%   print_output - boolean switch for fprintf output of values
%   make_graphs - boolean switch to output histograms of jackknife estimates in a four panel plot
%
% Outputs:
%   jack_pct_CI_values - an 1 by y vector of confidence interval values
%           values are percentile confidence intervals of Ne estimates based on jackknifing over r_squared estimates from allele pairs
%           1,1 - Ne from r_squared_c AFW lower bound CI
%           1,2 - Ne from r_squared_c AFW upper bound CI
%           1,3 - Ne from r_squared_c AFT lower bound CI
%           1,4 - Ne from r_squared_c AFT upper bound CI
%
%           1,5 - Ne from r_squared_delta AFW lower bound CI
%           1,6 - Ne from r_squared_delta AFW upper bound CI
%           1,7 - Ne from r_squared_delta AFT lower bound CI
%           1,8 - Ne from r_squared_delta AFT upper bound CI
%
%   jack_norm_CI_values - an 1 by y vector of confidence interval values
%           values are normal confidence intervals of Ne estimates based on jackknifing over r_squared estimates from allele pairs
%           1,1 - Ne from r_squared_c AFW lower bound CI
%           1,2 - Ne from r_squared_c AFW upper bound CI
%           1,3 - Ne from r_squared_c AFT lower bound CI
%           1,4 - Ne from r_squared_c AFT upper bound CI
%
%           1,5 - Ne from r_squared_delta AFW lower bound CI
%           1,6 - Ne from r_squared_delta AFW upper bound CI
%           1,7 - Ne from r_squared_delta AFT lower bound CI
%           1,8 - Ne from r_squared_delta AFT upper bound CI
%
%   chi_CI_values - an 1 by y vector of confidence interval values
%           values are confidence intervals of Ne estimates based on the chi-square distribution
%           r^2_{comp} AFW:
%           1,1 - Ne from r_squared_c AFW chi-squared n' lower bound CI
%           1,2 - Ne from r_squared_c AFW chi-squared n' upper bound CI
%           1,3 - Ne from r_squared_c AFW chi-squared locus pairs lower bound CI
%           1,4 - Ne from r_squared_c AFW chi-squared locus pairs upper bound CI
%           1,5 - Ne from r_squared_c AFW chi-squared number of loci lower bound CI
%           1,6 - Ne from r_squared_c AFW chi-squared number of loci upper bound CI
%
%           r^2_{comp} AFT:
%           1,7 - Ne from r_squared_c AFT chi-squared n' lower bound CI
%           1,8 - Ne from r_squared_c AFT chi-squared n' upper bound CI
%           1,9 - Ne from r_squared_c AFT chi-squared locus pairs lower bound CI
%           1,10 - Ne from r_squared_c AFT chi-squared locus pairs upper bound CI
%           1,11 - Ne from r_squared_c AFT chi-squared number of loci lower bound CI
%           1,12 - Ne from r_squared_c AFT chi-squared number of loci upper bound CI
%    
%           r^2_{delta} AFW:
%           1,13 - Ne from r_squared_delta AFW chi-squared n' lower bound CI
%           1,14 - Ne from r_squared_delta AFW chi-squared n' upper bound CI
%           1,15 - Ne from r_squared_delta AFW chi-squared locus pairs lower bound CI
%           1,16 - Ne from r_squared_delta AFW chi-squared locus pairs upper bound CI
%           1,17 - Ne from r_squared_delta AFW chi-squared number of loci lower bound CI
%           1,18 - Ne from r_squared_delta AFW chi-squared number of loci upper bound CI
%
%           r^2_{delta} AFW:
%           1,19 - Ne from r_squared_delta AFT chi-squared n' lower bound CI
%           1,20 - Ne from r_squared_delta AFT chi-squared n' upper bound CI
%           1,21 - Ne from r_squared_delta AFT chi-squared locus pairs lower bound CI
%           1,22 - Ne from r_squared_delta AFT chi-squared locus pairs upper bound CI
%           1,23 - Ne from r_squared_delta AFT chi-squared number of loci lower bound CI
%           1,24 - Ne from r_squared_delta AFT chi-squared number of loci upper bound CI
%
%   effective_n_c - a 1 by 2 vector of reals
%           1st column is AFW estimate n' value, 2nd column is AFT n' value
%   effective_n_delta - a 1 by 2 vector of reals
%           1st column is AFW estimate n' value, 2nd column is AFT n' value
%
%**************
%   Copyright 2017 Matthew B Hamilton.
%   
%   This file is part of SpEED-Ne.
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



% get parametric CIs based on Waples (2006) equation 12
chi_lower_c = chi2inv(alpha/2,c_locus_pairs);
chi_upper_c = chi2inv((1-alpha/2),c_locus_pairs);

c_W_lower_bound = (c_locus_pairs*r_squared_c_W)/chi_lower_c;
c_W_upper_bound = (c_locus_pairs*r_squared_c_W)/chi_upper_c;

c_TH_lower_bound = (c_locus_pairs*r_squared_c_TH)/chi_lower_c;
c_TH_upper_bound = (c_locus_pairs*r_squared_c_TH)/chi_upper_c;

chi_lower_delta = chi2inv(alpha/2,delta_locus_pairs);
chi_upper_delta = chi2inv((1-alpha/2),delta_locus_pairs);

delta_W_lower_bound = (delta_locus_pairs*r_squared_delta_W)/chi_lower_delta;
delta_W_upper_bound = (delta_locus_pairs*r_squared_delta_W)/chi_upper_delta;

delta_TH_lower_bound = (delta_locus_pairs*r_squared_delta_TH)/chi_lower_delta;
delta_TH_upper_bound = (delta_locus_pairs*r_squared_delta_TH)/chi_upper_delta;

if print_output
    fprintf('\n');
    fprintf('==========================================================================\n');
    fprintf('Parametric confidence intervals for r^2 estimates based on the\n');
    fprintf('chi-square distribution and number of pairs of loci in data set.\n\n');
    fprintf('\t\testimate\tlower bound\tupper bound\tlocus pairs\n');
    fprintf('--------------------------------------------------------------------------\n');

    fprintf('r^2_{comp} AFW\t%8.6g\t%8.6f\t%8.6f\t%8.6g\n\n', r_squared_c_W,c_W_upper_bound,c_W_lower_bound,c_locus_pairs);
    fprintf('r^2_{comp} AFT\t%8.6g\t%8.6f\t%8.6f\t%8.6g\n\n', r_squared_c_TH,c_TH_upper_bound,c_TH_lower_bound,c_locus_pairs);

    fprintf('r^2_{delta} AFW\t%8.6g\t%8.6f\t%8.6f\t%8.6g\n\n', r_squared_delta_W,delta_W_upper_bound,delta_W_lower_bound,delta_locus_pairs);
    fprintf('r^2_{delta} AFT\t%8.6g\t%8.6f\t%8.6f\t%8.6g\n', r_squared_delta_TH,delta_TH_upper_bound,delta_TH_lower_bound,delta_locus_pairs);

    fprintf('==========================================================================\n\n');
end % if print_output


%---------------------------
if num_loci <= 2 && print_output
    fprintf('**************\n');
    fprintf('Warning! Too few loci to compute Waples & Do "jackknife" chi-square CIs and effective sample size n''.\n');
    fprintf('**************\n');
end 

% Estimate CI with jackknife over locus pairs as implemented by Waples & Do (2007). Requires at least 3 loci!

% carry out jackknife for all pairwise r_squared_c AFW and AFT values
% delete one row from column vector of all estimates, make new estimate
j_r_squared_c = zeros(c_locus_pairs,2); % 1st column is AFW, 2nd column is AFT 
j_r_squared_c_denominator = zeros(c_locus_pairs,1);

for i=1:c_locus_pairs,
    del_one_squared_c_value = r_squared_c_values(~ismember(1:c_locus_pairs, [i]), :); % get subset of estimates with entry i deleted

    j_r_squared_c(i,1) = sum(del_one_squared_c_value(:,1)); % numerator for AFW estimates
    j_r_squared_c(i,2) = sum(del_one_squared_c_value(:,2)); % numerator for AFT estimates
    j_r_squared_c_denominator(i,1) = sum(del_one_squared_c_value(:,3));
end

j_r_squared_c(:,1) = j_r_squared_c(:,1)./j_r_squared_c_denominator; % get sub-estimate of r_squared_c_AFW
j_r_squared_c(:,2) = j_r_squared_c(:,2)./j_r_squared_c_denominator; % get sub-estimate of r_squared_c_AFT

[j_r_squared_c_W_lower,j_r_squared_c_W_upper] = percentile_CI(j_r_squared_c(:,1),c_locus_pairs,alpha); % get percentile CI for r_squared_c_AFW
[j_r_squared_c_TH_lower,j_r_squared_c_TH_upper] = percentile_CI(j_r_squared_c(:,2),c_locus_pairs,alpha); % get percentile CI for r_squared_c_AFT

[r_squared_c_AFW_norm_CI(1,1),r_squared_c_AFW_norm_CI(1,2)] = normal_dist_CI(r_squared_c_W,j_r_squared_c(:,1),c_locus_pairs,alpha); % get normal CI for r_squared_c_AFW
[r_squared_c_AFT_norm_CI(1,1),r_squared_c_AFT_norm_CI(1,2)] = normal_dist_CI(r_squared_c_TH,j_r_squared_c(:,2),c_locus_pairs,alpha); % get normal CI for r_squared_c_AFT


% carry out jackknife for all pairwise r_squared_delta AFW and AFT values
% delete one from column vector of all estimates, make new estimate
j_r_squared_delta = zeros(delta_locus_pairs,2); % 1st column is AFW, 2nd column is AFT 
j_r_squared_delta_denominator = zeros(delta_locus_pairs,1);

for i=1:c_locus_pairs,
    del_one_squared_delta_value = r_squared_delta_values(~ismember(1:delta_locus_pairs, [i]), :); % get subset of estimates with entry i deleted

    j_r_squared_delta(i,1) = sum(del_one_squared_delta_value(:,1)); % numerator for AFW estimates
    j_r_squared_delta(i,2) = sum(del_one_squared_delta_value(:,2)); % numerator for AFT estimates
    j_r_squared_delta_denominator(i,1) = sum(del_one_squared_delta_value(:,3));
end

j_r_squared_delta(:,1) = j_r_squared_delta(:,1)./j_r_squared_delta_denominator; % get sub-estimate of r_squared_c_AFW
j_r_squared_delta(:,2) = j_r_squared_delta(:,2)./j_r_squared_delta_denominator; % get sub-estimate of r_squared_c_AFT

[j_r_squared_delta_W_lower,j_r_squared_delta_W_upper] = percentile_CI(j_r_squared_delta(:,1),delta_locus_pairs,alpha); % get percentile CI for r_squared_delta_AFW
[j_r_squared_delta_TH_lower,j_r_squared_delta_TH_upper] = percentile_CI(j_r_squared_delta(:,2),delta_locus_pairs,alpha); % get percentile CI for r_squared_delta_AFT

[r_squared_delta_AFW_norm_CI(1,1),r_squared_delta_AFW_norm_CI(1,2)] = normal_dist_CI(r_squared_delta_W,j_r_squared_delta(:,1),delta_locus_pairs,alpha); % get normal CI for r_squared_delta_AFW
[r_squared_delta_AFT_norm_CI(1,1),r_squared_delta_AFT_norm_CI(1,2)] = normal_dist_CI(r_squared_delta_TH,j_r_squared_delta(:,2),delta_locus_pairs,alpha); % get normal CI for r_squared_delta_AFT


% make CIs for delete-one jackknife over all unique pairs of loci

% estimates for r^2_{comp} AFW
Ne_r_squared_c_AFW = 1/(3*(r_squared_c_W - median_r_sqr_c_AFW_perm));

% percentile
jack_pct_CI_values(1,1) = 1/(3*(j_r_squared_c_W_upper - median_r_sqr_c_AFW_perm)); % Ne from r_squared_c AFW lower bound CI
jack_pct_CI_values(1,2) = 1/(3*(j_r_squared_c_W_lower - median_r_sqr_c_AFW_perm)); % Ne from r_squared_c AFW upper bound CI

% normal distribution
jack_norm_CI_values(1,1) = 1/(3*(r_squared_c_AFW_norm_CI(1,2) - median_r_sqr_c_AFW_perm)); % Ne from r_squared_c AFW lower bound CI
jack_norm_CI_values(1,2) = 1/(3*(r_squared_c_AFW_norm_CI(1,1) - median_r_sqr_c_AFW_perm)); % Ne from r_squared_c AFW upper bound CI
%-----------------

% estimates for r^2_{comp} AFT
Ne_r_squared_c_AFT = 1/(3*(r_squared_c_TH - median_r_sqr_c_AFT_perm));

% percentile
jack_pct_CI_values(1,3) = 1/(3*(j_r_squared_c_TH_upper - median_r_sqr_c_AFT_perm)); % Ne from r_squared_c AFW lower bound CI
jack_pct_CI_values(1,4) = 1/(3*(j_r_squared_c_TH_lower - median_r_sqr_c_AFT_perm)); % Ne from r_squared_c AFW upper bound CI

% normal distribution
jack_norm_CI_values(1,3) = 1/(3*(r_squared_c_AFT_norm_CI(1,2) - median_r_sqr_c_AFT_perm)); % Ne from r_squared_c AFW lower bound CI
jack_norm_CI_values(1,4) = 1/(3*(r_squared_c_AFT_norm_CI(1,1) - median_r_sqr_c_AFT_perm)); % Ne from r_squared_c AFW upper bound CI
%-----------------

% estimates for r^2_{delta} AFW
Ne_r_squared_delta_AFW = 1/(3*(r_squared_c_W - correction_factor_har_mean_S));    

% percentile
jack_pct_CI_values(1,5) = 1/(3*(j_r_squared_delta_W_upper - correction_factor_har_mean_S)); % Ne from r_squared_c AFW lower bound CI
jack_pct_CI_values(1,6) = 1/(3*(j_r_squared_delta_W_lower - correction_factor_har_mean_S)); % Ne from r_squared_c AFW upper bound CI    

% normal distribution
jack_norm_CI_values(1,5) = 1/(3*(r_squared_delta_AFW_norm_CI(1,2) - correction_factor_har_mean_S)); % Ne from r_squared_c AFW lower bound CI
jack_norm_CI_values(1,6) = 1/(3*(r_squared_delta_AFW_norm_CI(1,1) - correction_factor_har_mean_S)); % Ne from r_squared_c AFW upper bound CI
%-----------------

% estimates for r^2_{delta} AFT
Ne_r_squared_delta_AFT = 1/(3*(r_squared_delta_TH - correction_factor_har_mean_S));

% percentile
jack_pct_CI_values(1,7) = 1/(3*(j_r_squared_delta_TH_upper - correction_factor_har_mean_S)); % Ne from r_squared_delta AFT lower bound CI
jack_pct_CI_values(1,8) = 1/(3*(j_r_squared_delta_TH_lower - correction_factor_har_mean_S)); % Ne from r_squared_delta AFT upper bound CI

% normal distribution
jack_norm_CI_values(1,7) = 1/(3*(r_squared_delta_AFT_norm_CI(1,2) - correction_factor_har_mean_S)); % Ne from r_squared_delta AFT lower bound CI
jack_norm_CI_values(1,8) = 1/(3*(r_squared_delta_AFT_norm_CI(1,1) - correction_factor_har_mean_S)); % Ne from r_squared_delta AFT upper bound CI
%-----------------


% make table for delete-one jackknife over all unique pairs of loci
if print_output

    fprintf('\n');
    fprintf('====================================================================================\n');
    fprintf('Confidence intervals from a delete-one jackknife over r^2 estimates over all \n');
    fprintf('unique pairs of loci.\n');
    fprintf('Ne = 1/(3*(r^2 - r^2 correction factor)).\n');
    fprintf('The r^2 correction factor is median r^2 permute for r^2_{c} and harmonic mean S \n');
    fprintf('for r^2_{delta}. Estimates are either allele frequency weighted (AFW) or allele \n');
    fprintf('frequency thresholded (AFT).\n\n');
    fprintf('\t\testimate\tlower bound\tupper bound\tlocus pairs\n');
    fprintf('------------------------------------------------------------------------------------\n');

    %-----------------------------
    % estimates for r^2_{comp} AFW
    fprintf('r^2_{comp} AFW:\n');

    % percentile
    fprintf('  percentile CI: %8.6g\t%8.6f\t%8.6f\t%8.6g\n',r_squared_c_W,j_r_squared_c_W_lower,j_r_squared_c_W_upper,c_locus_pairs);

    if jack_pct_CI_values(1,2) >= 0 
        fprintf('  Ne CI:\t%8.6g\t%8.6g\t%8.6g\n', Ne_r_squared_c_AFW, jack_pct_CI_values(1,1), jack_pct_CI_values(1,2));
    else
        fprintf('  Ne CI:\t%8.6g\t%8.6g\t%s\n', Ne_r_squared_c_AFW, jack_pct_CI_values(1,1), 'infinity');
    end
    fprintf('\n');

    % normal distribution
    fprintf('  normal CI:\t%8.6g\t%8.6f\t%8.6f\t%8.6g\n',r_squared_c_W,r_squared_c_AFW_norm_CI(1,1),r_squared_c_AFW_norm_CI(1,2),c_locus_pairs);

    if jack_norm_CI_values(1,2) >= 0 
        fprintf('  Ne CI:\t%8.6g\t%8.6g\t%8.6g\n', Ne_r_squared_c_AFW, jack_norm_CI_values(1,1), jack_norm_CI_values(1,2));
    else
        fprintf('  Ne CI:\t%8.6g\t%8.6g\t%s\n', Ne_r_squared_c_AFW, jack_norm_CI_values(1,1), 'infinity');
    end
    fprintf('\n');

    %-----------------------------
    % estimates for r^2_{comp} AFT
    fprintf('r^2_{comp} AFT:\n');

    % percentile
    fprintf('  percentile CI: %8.6g\t%8.6f\t%8.6f\t%8.6g\n',r_squared_c_TH,j_r_squared_c_TH_lower,j_r_squared_c_TH_upper,c_locus_pairs);

    if jack_pct_CI_values(1,4) >= 0 
        fprintf('  Ne CI:\t%8.6g\t%8.6g\t%8.6g\n', Ne_r_squared_c_AFT, jack_pct_CI_values(1,3), jack_pct_CI_values(1,4));
    else
        fprintf('  Ne CI:\t%8.6g\t%8.6g\t%s\n', Ne_r_squared_c_AFT, jack_pct_CI_values(1,3), 'infinity');
    end
    fprintf('\n');

    % normal distribution
    fprintf('  normal CI:\t%8.6g\t%8.6f\t%8.6f\t%8.6g\n',r_squared_c_TH,r_squared_c_AFT_norm_CI(1,1),r_squared_c_AFT_norm_CI(1,2),c_locus_pairs);

    if jack_norm_CI_values(1,4) >= 0 
        fprintf('  Ne CI:\t%8.6g\t%8.6g\t%8.6g\n', Ne_r_squared_c_AFT, jack_norm_CI_values(1,3), jack_norm_CI_values(1,4));
    else
        fprintf('  Ne CI:\t%8.6g\t%8.6g\t%s\n', Ne_r_squared_c_AFT, jack_norm_CI_values(1,3), 'infinity');
    end
    fprintf('\n');

    %-----------------------------
    % estimates for r^2_{delta} AFW
    fprintf('r^2_{delta} AFW:\n');

    % percentile
    fprintf('  percentile CI: %8.6g\t%8.6f\t%8.6f\t%8.6g\n',r_squared_delta_W,j_r_squared_delta_W_lower,j_r_squared_delta_W_upper,delta_locus_pairs);

    if jack_pct_CI_values(1,6) >= 0 
        fprintf('  Ne CI:\t%8.6g\t%8.6g\t%8.6g\n', Ne_r_squared_delta_AFW, jack_pct_CI_values(1,1), jack_pct_CI_values(1,2));
    else
        fprintf('  Ne CI:\t%8.6g\t%8.6g\t%s\n', Ne_r_squared_delta_AFW, jack_pct_CI_values(1,1), 'infinity');
    end
    fprintf('\n');

    % normal distribution
    fprintf('  normal CI:\t%8.6g\t%8.6f\t%8.6f\t%8.6g\n',r_squared_delta_W,r_squared_delta_AFW_norm_CI(1,1),r_squared_delta_AFW_norm_CI(1,2),delta_locus_pairs);

    if jack_norm_CI_values(1,6) >= 0 
        fprintf('  Ne CI:\t%8.6g\t%8.6g\t%8.6g\n', Ne_r_squared_delta_AFW, jack_norm_CI_values(1,1), jack_norm_CI_values(1,2));
    else
        fprintf('  Ne CI:\t%8.6g\t%8.6g\t%s\n', Ne_r_squared_delta_AFW, jack_norm_CI_values(1,1), 'infinity');
    end
    fprintf('\n');

    %-----------------------------
    % estimates for r^2_{delta} AFT
    fprintf('r^2_{delta} AFT:\n');

    % percentile
    fprintf('  percentile CI: %8.6g\t%8.6f\t%8.6f\t%8.6g\n',r_squared_delta_TH,j_r_squared_delta_TH_lower,j_r_squared_delta_TH_upper,delta_locus_pairs);

    if jack_pct_CI_values(1,8) >= 0 
        fprintf('  Ne CI:\t%8.6g\t%8.6g\t%8.6g\n', Ne_r_squared_delta_AFT, jack_pct_CI_values(1,7), jack_pct_CI_values(1,8));
    else
        fprintf('  Ne CI:\t%8.6g\t%8.6g\t%s\n', Ne_r_squared_delta_AFT, jack_pct_CI_values(1,7), 'infinity');
    end
    fprintf('\n');

    % normal distribution
    fprintf('  normal CI:\t%8.6g\t%8.6f\t%8.6f\t%8.6g\n',r_squared_delta_TH,r_squared_delta_AFT_norm_CI(1,1),r_squared_delta_AFT_norm_CI(1,2),delta_locus_pairs);

    if jack_norm_CI_values(1,8) >= 0 
        fprintf('  Ne CI:\t%8.6g\t%8.6g\t%8.6g\n', Ne_r_squared_delta_AFT, jack_norm_CI_values(1,7), jack_norm_CI_values(1,8));
    else
        fprintf('  Ne CI:\t%8.6g\t%8.6g\t%s\n', Ne_r_squared_delta_AFT, jack_norm_CI_values(1,7), 'infinity');
    end
    fprintf('\n');


    fprintf('------------------------------------------------------------------------------------\n');
    fprintf('Confidence intervals are %g%% to %g%%\n',100*alpha/2,100*(1-alpha/2));

    fprintf('====================================================================================\n');

end; % if print_output



if make_graphs
    % show distribution of jackknifed estimates to evalute distribution
    figure;
    hold on;
    subplot(2,2,1);
    hist(j_r_squared_c(:,1)); % 1st column is AFW, 2nd column is AFT
    xlabel('r_c^2 AFW')
    ylabel('Count')

    subplot(2,2,2);
    hist(j_r_squared_c(:,2)); % 1st column is AFW, 2nd column is AFT
    xlabel('r_c^2 AFT')
    ylabel('Count')

    subplot(2,2,3)
    hist(j_r_squared_delta(:,1)); % 1st column is AFW, 2nd column is AFT
    xlabel('r_{\Delta}^2 AFW')
    ylabel('Count')

    subplot(2,2,4)
    hist(j_r_squared_delta(:,2)); % 1st column is AFW, 2nd column is AFT
    xlabel('r_{\Delta}^2 AFT')
    ylabel('Count')

    suptitle('Distributions of delete-one jackknife r^2 sub-estimates from all pairs of alleles.');
    hold off;
end; % if print_graphs


%---------------------------
% Estimate variance from "jackknifed" estimates of j_r_squared_c and j_r_squared_delta vectors
% following Waples & Do (2007) p. 754-755.  This provides an estimate
% of effective number of independent comparisons or n'.
% All computations include both AFW and AFT r_squared estimates 
% (stored as 1st and 2nd columns in both j_r_squared_c and j_r_squared_delta). 

% compute n' for r_squared_comp AFW and AFT estimates
mean_r_c_j = mean(j_r_squared_c);
matrix_mean_r_c_j = repmat(mean_r_c_j,c_locus_pairs,1);
sqr_diffs = (j_r_squared_c - matrix_mean_r_c_j).^2;

var_j_c = (c_locus_pairs-1)/c_locus_pairs * sum(sqr_diffs);

theta_c(1,1) = var_j_c(1,1)/(r_squared_c_W.^2); % theta for AFW estimate
theta_c(1,2) = var_j_c(1,2)/(r_squared_c_TH.^2); % theta for AFT estimate

effective_n_c = 2./theta_c; % 1st column is AFW estimate value, 2nd column is AFT value

% compute n' for r_squared_delta AFW and AFT estimates
mean_r_delta_j = mean(j_r_squared_delta);
matrix_mean_r_delta_j = repmat(mean_r_delta_j,delta_locus_pairs,1);
sqr_diffs = (j_r_squared_delta - matrix_mean_r_delta_j).^2;

var_j_delta = (delta_locus_pairs-1)/delta_locus_pairs * sum(sqr_diffs);

theta_delta(1,1) = var_j_delta(1,1)/(r_squared_delta_W.^2); % theta for AFW estimate
theta_delta(1,2) = var_j_delta(1,2)/(r_squared_delta_TH.^2); % theta for AFT estimate

effective_n_delta = 2./theta_delta; % 1st column is AFW estimate value, 2nd column is AFT value


% compute estimates and make table for Chi-square confidence intervals and estimate of n'
% from delete-one jackknife variance estimate.
if print_output
    fprintf('\n');
    fprintf('==========================================================================\n');
    fprintf('Chi-square confidence intervals. \n');
    fprintf('Estimate of effective number of independent comparisons (n'') based on\n');
    fprintf('delete-one-locus-pair jackknike estimate of variance of r^2 (following Waples and Do 2007\n');
    fprintf('p. 754-755). Estimates are either allele frequency weighted (AFW) or allele\n');
    fprintf('frequency thresholded (AFT). Ne = 1/(3*(r^2 - r^2 correction factor)).\n');
    fprintf('The r^2 correction factor is median r^2 permute for r^2_{c} and harmonic mean S for r^2_{delta}.\n');
    fprintf('Compare n'' to number of loci, and number of locus pairs.\n');
    fprintf('A value of n'' > number of locus pairs suggests violation of assumptions\n');
    fprintf('and that confidence intervals based on n'' are too small.\n\n');
    fprintf('\t\testimate\tlower CI\t upper CI\n')
    fprintf('--------------------------------------------------------------------------\n');
end; %if print_output

%===================================
%------------- CIs and Ne's for r^2_{comp} AFW

%chi-squared CI and Ne using n' as sample size
chi_lower = chi2inv(alpha/2,effective_n_c(1,1));
chi_upper = chi2inv(1-(alpha/2),effective_n_c(1,1));

lower_bound = (effective_n_c(1,1)*r_squared_c_W)/chi_upper;
upper_bound = (effective_n_c(1,1)*r_squared_c_W)/chi_lower;    

Ne_r_squared_c_W = 1/(3*(r_squared_c_W - median_r_sqr_c_AFW_perm));
upper_bound_Ne = 1/(3*(lower_bound - median_r_sqr_c_AFW_perm));
lower_bound_Ne = 1/(3*(upper_bound - median_r_sqr_c_AFW_perm));

chi_CI_values(1,1) = lower_bound_Ne; %Ne from r_squared_c AFW chi-squared n' lower bound CI
chi_CI_values(1,2) = upper_bound_Ne; %Ne from r_squared_c AFW chi-squared n' upper bound CI

if print_output
    fprintf('r^2_{comp} AFW:\n');

    fprintf('\tn'' = %8.6g\n',effective_n_c(1,1));
    fprintf('\tr^2 = ');
    fprintf('\t%8.6f\t%8.6f\t%8.6g\n',r_squared_c_W,lower_bound,upper_bound);

    if upper_bound_Ne >= 0 
        fprintf('\tNe = \t%8.6g\t%8.6g\t%8.6g\n', Ne_r_squared_c_W, lower_bound_Ne, upper_bound_Ne);
    else
        fprintf('\tNe = \t%8.6g\t%8.6g\t%s\n', Ne_r_squared_c_W, lower_bound_Ne, 'infinity');
    end
    fprintf('\n');
end % if print_output

%chi-squared CI and Ne using locus pairs as sample size
chi_lower = chi2inv(alpha/2,c_locus_pairs);
chi_upper = chi2inv((1-alpha/2),c_locus_pairs);

lower_bound = (c_locus_pairs*r_squared_c_W)/chi_upper;
upper_bound = (c_locus_pairs*r_squared_c_W)/chi_lower;

upper_bound_Ne = 1/(3*(lower_bound - median_r_sqr_c_AFW_perm));
lower_bound_Ne = 1/(3*(upper_bound - median_r_sqr_c_AFW_perm));

chi_CI_values(1,3) = lower_bound_Ne; %Ne from r_squared_c AFW chi-squared locus pairs lower bound CI
chi_CI_values(1,4) = upper_bound_Ne; %Ne from r_squared_c AFW chi-squared locus pairs upper bound CI

if print_output
    fprintf('\tlocus pairs = %i\n',c_locus_pairs);
    fprintf('\tr^2 = ');
    fprintf('\t%8.6f\t%8.6f\t%8.6g\n',r_squared_c_W,lower_bound,upper_bound);

    if upper_bound_Ne >= 0 
        fprintf('\tNe = \t%8.6g\t%8.6g\t%8.6g\n', Ne_r_squared_c_W, lower_bound_Ne, upper_bound_Ne);
    else
        fprintf('\tNe = \t%8.6g\t%8.6g\t%s\n', Ne_r_squared_c_W, lower_bound_Ne, 'infinity');
    end
    fprintf('\n');
end % if print_output

%chi-squared CI and Ne using number of loci as sample size
chi_lower = chi2inv(alpha/2,num_loci);
chi_upper = chi2inv((1-alpha/2),num_loci);

lower_bound = (num_loci*r_squared_c_W)/chi_upper;
upper_bound = (num_loci*r_squared_c_W)/chi_lower;

upper_bound_Ne = 1/(3*(lower_bound - median_r_sqr_c_AFW_perm));
lower_bound_Ne = 1/(3*(upper_bound - median_r_sqr_c_AFW_perm));

chi_CI_values(1,5) = lower_bound_Ne; %Ne from r_squared_c AFW chi-squared number of loci lower bound CI
chi_CI_values(1,6) = upper_bound_Ne; %Ne from r_squared_c AFW chi-squared number of loci upper bound CI

if print_output
    fprintf('\tnumber of loci = %i\n',num_loci);
    fprintf('\tr^2 = ');
    fprintf('\t%8.6f\t%8.6f\t%8.6g\n',r_squared_c_W,lower_bound,upper_bound);

    if upper_bound_Ne >= 0 
        fprintf('\tNe = \t%8.6g\t%8.6g\t%8.6g\n', Ne_r_squared_c_W, lower_bound_Ne, upper_bound_Ne);
    else
        fprintf('\tNe = \t%8.6g\t%8.6g\t%s\n', Ne_r_squared_c_W, lower_bound_Ne, 'infinity');
    end
    fprintf('\n');
end % if print_output


%===================================
%------------- CIs and Ne's for r^2_{comp} AFT

%chi-squared CI and Ne using n' as sample size
chi_lower = chi2inv(alpha/2,effective_n_c(1,2));
chi_upper = chi2inv(1-(alpha/2),effective_n_c(1,2));

lower_bound = (effective_n_c(1,2)*r_squared_c_TH)/chi_upper;
upper_bound = (effective_n_c(1,2)*r_squared_c_TH)/chi_lower;

Ne_r_squared_c_TH = 1/(3*(r_squared_c_TH - median_r_sqr_c_AFT_perm));
upper_bound_Ne = 1/(3*(lower_bound - median_r_sqr_c_AFT_perm));
lower_bound_Ne = 1/(3*(upper_bound - median_r_sqr_c_AFT_perm));

chi_CI_values(1,7) = lower_bound_Ne; %Ne from r_squared_c AFT chi-squared n' lower bound CI
chi_CI_values(1,8) = upper_bound_Ne; %Ne from r_squared_c AFT chi-squared n' upper bound CI

if print_output
    fprintf('r^2_{comp} AFT:\n');

    fprintf('\tn'' = %8.6g\n',effective_n_c(1,2));
    fprintf('\tr^2 = ');
    fprintf('\t%8.6f\t%8.6f\t%8.6g\n',r_squared_c_TH,lower_bound,upper_bound);

    if upper_bound_Ne >= 0 
        fprintf('\tNe = \t%8.6g\t%8.6g\t%8.6g\n', Ne_r_squared_c_TH, lower_bound_Ne, upper_bound_Ne);
    else
        fprintf('\tNe = \t%8.6g\t%8.6g\t%s\n', Ne_r_squared_c_TH, lower_bound_Ne, 'infinity');
    end
    fprintf('\n');
end; % if print_output

%chi-squared CI and Ne using locus pairs as sample size
chi_lower = chi2inv(alpha/2,c_locus_pairs);
chi_upper = chi2inv((1-alpha/2),c_locus_pairs);

lower_bound = (c_locus_pairs*r_squared_c_TH)/chi_upper;
upper_bound = (c_locus_pairs*r_squared_c_TH)/chi_lower;

upper_bound_Ne = 1/(3*(lower_bound - median_r_sqr_c_AFT_perm));
lower_bound_Ne = 1/(3*(upper_bound - median_r_sqr_c_AFT_perm));

chi_CI_values(1,9) = lower_bound_Ne; %Ne from r_squared_c AFT chi-squared locus pairs lower bound CI
chi_CI_values(1,10) = upper_bound_Ne; %Ne from r_squared_c AFT chi-squared locus pairs upper bound CI

if print_output
    fprintf('\tlocus pairs = %i\n',c_locus_pairs);
    fprintf('\tr^2 = ');
    fprintf('\t%8.6f\t%8.6f\t%8.6g\n',r_squared_c_TH,lower_bound,upper_bound);

    if upper_bound_Ne >= 0 
        fprintf('\tNe = \t%8.6g\t%8.6g\t%8.6g\n', Ne_r_squared_c_TH, lower_bound_Ne, upper_bound_Ne);
    else
        fprintf('\tNe = \t%8.6g\t%8.6g\t%s\n', Ne_r_squared_c_TH, lower_bound_Ne, 'infinity');
    end
    fprintf('\n');
end; % if print_output

%chi-squared CI and Ne using number of loci as sample size
chi_lower = chi2inv(alpha/2,num_loci);
chi_upper = chi2inv((1-alpha/2),num_loci);

lower_bound = (num_loci*r_squared_c_TH)/chi_upper;
upper_bound = (num_loci*r_squared_c_TH)/chi_lower;

upper_bound_Ne = 1/(3*(lower_bound - median_r_sqr_c_AFT_perm));
lower_bound_Ne = 1/(3*(upper_bound - median_r_sqr_c_AFT_perm));

chi_CI_values(1,11) = lower_bound_Ne; %Ne from r_squared_c AFW chi-squared number of loci lower bound CI
chi_CI_values(1,12) = upper_bound_Ne; %Ne from r_squared_c AFW chi-squared number of loci upper bound CI

if print_output
    fprintf('\tnumber of loci = %i\n',num_loci);
    fprintf('\tr^2 = ');
    fprintf('\t%8.6f\t%8.6f\t%8.6g\n',r_squared_c_TH,lower_bound,upper_bound);

    if upper_bound_Ne >= 0 
        fprintf('\tNe = \t%8.6g\t%8.6g\t%8.6g\n', Ne_r_squared_c_TH, lower_bound_Ne, upper_bound_Ne);
    else
        fprintf('\tNe = \t%8.6g\t%8.6g\t%s\n', Ne_r_squared_c_TH, lower_bound_Ne, 'infinity');
    end
    fprintf('\n');
end; % if print_output

%===================================
%------------- CIs and Ne's for r^2_{delta} AFW

%chi-squared CI and Ne using n' as sample size
chi_lower = chi2inv(alpha/2,effective_n_delta(1,1));
chi_upper = chi2inv(1-(alpha/2),effective_n_delta(1,1));

lower_bound = (effective_n_delta(1,1)*r_squared_delta_W)/chi_upper;
upper_bound = (effective_n_delta(1,1)*r_squared_delta_W)/chi_lower;

Ne_r_squared_delta_W = 1/(3*(r_squared_delta_W - correction_factor_har_mean_S));
upper_bound_Ne = 1/(3*(lower_bound - correction_factor_har_mean_S));
lower_bound_Ne = 1/(3*(upper_bound - correction_factor_har_mean_S));

chi_CI_values(1,13) = lower_bound_Ne; %Ne from r_squared_delta AFW chi-squared n' lower bound CI
chi_CI_values(1,14) = upper_bound_Ne; %Ne from r_squared_delta AFW chi-squared n' upper bound CI

if print_output
    fprintf('r^2_{delta} AFW:\n');

    fprintf('\tn'' = %8.6g\n',effective_n_delta(1,1));
    fprintf('\tr^2 = ');
    fprintf('\t%8.6f\t%8.6f\t%8.6g\n',r_squared_c_W,lower_bound,upper_bound);

    if upper_bound_Ne >= 0 
        fprintf('\tNe = \t%8.6g\t%8.6g\t%8.6g\n', Ne_r_squared_delta_W, lower_bound_Ne, upper_bound_Ne);
    else
        fprintf('\tNe = \t%8.6g\t%8.6g\t%s\n', Ne_r_squared_delta_W, lower_bound_Ne, 'infinity');
    end
    fprintf('\n');
end; % if print_output

%chi-squared CI and Ne using locus pairs as sample size
chi_lower = chi2inv(alpha/2,delta_locus_pairs);
chi_upper = chi2inv((1-alpha/2),delta_locus_pairs);

lower_bound = (delta_locus_pairs*r_squared_delta_W)/chi_upper;
upper_bound = (delta_locus_pairs*r_squared_delta_W)/chi_lower;

upper_bound_Ne = 1/(3*(lower_bound - correction_factor_har_mean_S));
lower_bound_Ne = 1/(3*(upper_bound - correction_factor_har_mean_S));

chi_CI_values(1,15) = lower_bound_Ne; %Ne from r_squared_delta AFW chi-squared locus pairs lower bound CI
chi_CI_values(1,16) = upper_bound_Ne; %Ne from r_squared_delta AFW chi-squared locus pairs upper bound CI

if print_output
    fprintf('\tlocus pairs = %i\n',delta_locus_pairs);
    fprintf('\tr^2 = ');
    fprintf('\t%8.6f\t%8.6f\t%8.6g\n',r_squared_delta_W,lower_bound,upper_bound);

    if upper_bound_Ne >= 0 
        fprintf('\tNe = \t%8.6g\t%8.6g\t%8.6g\n', Ne_r_squared_delta_W, lower_bound_Ne, upper_bound_Ne);
    else
        fprintf('\tNe = \t%8.6g\t%8.6g\t%s\n', Ne_r_squared_delta_W, lower_bound_Ne, 'infinity');
    end
    fprintf('\n');
end; % if print_output

%chi-squared CI and Ne using number of loci as sample size
chi_lower = chi2inv(alpha/2,num_loci);
chi_upper = chi2inv((1-alpha/2),num_loci);

lower_bound = (num_loci*r_squared_delta_W)/chi_upper;
upper_bound = (num_loci*r_squared_delta_W)/chi_lower;

upper_bound_Ne = 1/(3*(lower_bound - correction_factor_har_mean_S));
lower_bound_Ne = 1/(3*(upper_bound - correction_factor_har_mean_S));

chi_CI_values(1,17) = lower_bound_Ne; %Ne from r_squared_delta AFW chi-squared number of loci lower bound CI
chi_CI_values(1,18) = upper_bound_Ne; %Ne from r_squared_delta AFW chi-squared number of loci upper bound CI

if print_output
    fprintf('\tnumber of loci = %i\n',num_loci);
    fprintf('\tr^2 = ');
    fprintf('\t%8.6f\t%8.6f\t%8.6g\n',r_squared_delta_W,lower_bound,upper_bound);

    if upper_bound_Ne >= 0 
        fprintf('\tNe = \t%8.6g\t%8.6g\t%8.6g\n', Ne_r_squared_delta_W, lower_bound_Ne, upper_bound_Ne);
    else
        fprintf('\tNe = \t%8.6g\t%8.6g\t%s\n', Ne_r_squared_delta_W, lower_bound_Ne, 'infinity');
    end
    fprintf('\n');
end; % if print_output

%===================================
%CIs and Ne's for r^2_{delta} AFT

%chi-squared CI and Ne using n' as sample size
chi_lower = chi2inv(alpha/2,effective_n_delta(1,2));
chi_upper = chi2inv(1-(alpha/2),effective_n_delta(1,2));

lower_bound = (effective_n_delta(1,2)*r_squared_delta_TH)/chi_upper;
upper_bound = (effective_n_delta(1,2)*r_squared_delta_TH)/chi_lower;

Ne_r_squared_delta_TH = 1/(3*(r_squared_delta_TH - correction_factor_har_mean_S));
upper_bound_Ne = 1/(3*(lower_bound - correction_factor_har_mean_S));
lower_bound_Ne = 1/(3*(upper_bound - correction_factor_har_mean_S));

chi_CI_values(1,19) = lower_bound_Ne; %Ne from r_squared_delta AFT chi-squared n' lower bound CI
chi_CI_values(1,20) = upper_bound_Ne; %Ne from r_squared_delta AFT chi-squared n' upper bound CI

if print_output
    fprintf('r^2_{delta} AFT:\n');

    fprintf('\tn'' = %8.6g\n',effective_n_delta(1,2));
    fprintf('\tr^2 = ');
    fprintf('\t%8.6f\t%8.6f\t%8.6g\n',r_squared_delta_TH,lower_bound,upper_bound);

    if upper_bound_Ne >= 0 
        fprintf('\tNe = \t%8.6g\t%8.6g\t%8.6g\n', Ne_r_squared_delta_TH, lower_bound_Ne, upper_bound_Ne);
    else
        fprintf('\tNe = \t%8.6g\t%8.6g\t%s\n', Ne_r_squared_delta_TH, lower_bound_Ne, 'infinity');
    end
    fprintf('\n');
end; % if print_output

%chi-squared CI and Ne using locus pairs as sample size
chi_lower = chi2inv(alpha/2,delta_locus_pairs);
chi_upper = chi2inv((1-alpha/2),delta_locus_pairs);

lower_bound = (c_locus_pairs*r_squared_delta_TH)/chi_upper;
upper_bound = (c_locus_pairs*r_squared_delta_TH)/chi_lower;

upper_bound_Ne = 1/(3*(lower_bound - correction_factor_har_mean_S));
lower_bound_Ne = 1/(3*(upper_bound - correction_factor_har_mean_S));

chi_CI_values(1,21) = lower_bound_Ne; %Ne from r_squared_delta AFT chi-squared locus pairs lower bound CI
chi_CI_values(1,22) = upper_bound_Ne; %Ne from r_squared_delta AFT chi-squared locus pairs upper bound CI

if print_output
    fprintf('\tlocus pairs = %i\n',delta_locus_pairs);
    fprintf('\tr^2 = ');
    fprintf('\t%8.6f\t%8.6f\t%8.6g\n',r_squared_delta_TH,lower_bound,upper_bound);

    if upper_bound_Ne >= 0 
        fprintf('\tNe = \t%8.6g\t%8.6g\t%8.6g\n', Ne_r_squared_delta_TH, lower_bound_Ne, upper_bound_Ne);
    else
        fprintf('\tNe = \t%8.6g\t%8.6g\t%s\n', Ne_r_squared_delta_TH, lower_bound_Ne, 'infinity');
    end
    fprintf('\n');
end; % if print_output

%chi-squared CI and Ne using number of loci as sample size
chi_lower = chi2inv(alpha/2,num_loci);
chi_upper = chi2inv((1-alpha/2),num_loci);

lower_bound = (num_loci*r_squared_delta_TH)/chi_upper;
upper_bound = (num_loci*r_squared_delta_TH)/chi_lower;

upper_bound_Ne = 1/(3*(lower_bound - correction_factor_har_mean_S));
lower_bound_Ne = 1/(3*(upper_bound - correction_factor_har_mean_S));

chi_CI_values(1,23) = lower_bound_Ne; %Ne from r_squared_delta AFT chi-squared number of loci lower bound CI
chi_CI_values(1,24) = upper_bound_Ne; %Ne from r_squared_delta AFT chi-squared number of loci upper bound CI

if print_output
    fprintf('\tnumber of loci = %i\n',num_loci);
    fprintf('\tr^2 = ');
    fprintf('\t%8.6f\t%8.6f\t%8.6g\n',r_squared_delta_TH,lower_bound,upper_bound);

    if upper_bound_Ne >= 0 
        fprintf('\tNe = \t%8.6g\t%8.6g\t%8.6g\n', Ne_r_squared_delta_TH, lower_bound_Ne, upper_bound_Ne);
    else
        fprintf('\tNe = \t%8.6g\t%8.6g\t%s\n', Ne_r_squared_delta_TH, lower_bound_Ne, 'infinity');
    end
    fprintf('\n');

    fprintf('==========================================================================\n');
end; % if print_output


