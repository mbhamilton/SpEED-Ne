function post_process_ci

file_name = '/Users/matthewhamilton/Dropbox/R project/CI sims/Ne_250_L_12_1000_iter_TF.mat';
load(file_name);

fprintf('output for:\n');
fprintf('%s\n', file_name);

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


% Ne(iter,1) - r_squared_c AFW
% Ne(iter,2) - r_squared_c AFT
% Ne(iter,3) - r_squared_c AFW UB
% Ne(iter,4) - r_squared_c AFT UB
% Ne(iter,5) - r_squared_delta AFW
% Ne(iter,6) - r_squared_delta AFT
% Ne(iter,7) - r_squared_delta AFW UB
% Ne(iter,8) - r_squared_delta AFT UB
% Ne(iter,9) - Waples r_squared_delta AFW
% Ne(iter,10) - Waples r_squared_delta AFT
% Ne(iter,11) - Waples r_squared_delta AFW UB
% Ne(iter,12) - Waples r_squared_delta AFT UB
% Ne(iter,13) - r_squared_c AFW permute
% Ne(iter,14) - r_squared_c AFT permute
% Ne(iter,15) - r_squared_c AFW UB permute
% Ne(iter,16) - r_squared_c AFT UB permute


    % Outputs:
%   jack_allele_pct_CI_values - an 1 by y vector of confidence interval values
%           values are percentile confidence intervals of Ne estimates based on jackknifing over r_squared estimates from allele pairs
%           1,5 - Ne from r_squared_delta AFW lower bound CI
%           1,6 - Ne from r_squared_delta AFW upper bound CI
%           1,7 - Ne from r_squared_delta AFT lower bound CI
%           1,8 - Ne from r_squared_delta AFT upper bound CI
%
%   jack_allele_norm_CI_values - an 1 by y vector of confidence interval values
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

    
    
    
    
    
    
    
    % compute Waples Ne with jackknife over individuals r_sqr_delta upper and lower CIs
    % first column contains upper bound and second column contains lower bound
    [Waples_Ne_AFW_ji_pct_CI(:,1)] = Waples_Ne_estimate(har_mean_S_delta,r_sqr_d_AFW_jack_individ_CI_pct(:,1)); % col vector of lower CI Waples Ne from r_sqr_delta_AFW_jack_individ percentile
    [Waples_Ne_AFW_ji_pct_CI(:,2)] = Waples_Ne_estimate(har_mean_S_delta,r_sqr_d_AFW_jack_individ_CI_pct(:,2)); % col vector of upper CI Waples Ne from r_sqr_delta_AFW_jack_individ percentile

    [Waples_Ne_AFT_ji_pct_CI(:,1)] = Waples_Ne_estimate(har_mean_S_delta,r_sqr_d_AFT_jack_individ_CI_pct(:,1)); % col vector of lower CI Waples Ne from r_sqr_delta_AFW_jack_individ percentile
    [Waples_Ne_AFT_ji_pct_CI(:,2)] = Waples_Ne_estimate(har_mean_S_delta,r_sqr_d_AFT_jack_individ_CI_pct(:,2)); % col vector of upper CI Waples Ne from r_sqr_delta_AFW_jack_individ percentile

    [Waples_Ne_AFW_ji_norm_CI(:,1)] = Waples_Ne_estimate(har_mean_S_delta,r_sqr_d_AFW_jack_individ_CI_norm(:,1)); % col vector of lower CI Waples Ne from r_sqr_delta_AFW_jack_individ normal
    [Waples_Ne_AFW_ji_norm_CI(:,2)] = Waples_Ne_estimate(har_mean_S_delta,r_sqr_d_AFW_jack_individ_CI_norm(:,2)); % col vector of upper CI Waples Ne from r_sqr_delta_AFW_jack_individ normal
    
    [Waples_Ne_AFT_ji_norm_CI(:,1)] = Waples_Ne_estimate(har_mean_S_delta,r_sqr_d_AFT_jack_individ_CI_norm(:,1)); % col vector of lower CI Waples Ne from r_sqr_delta_AFW_jack_individ normal
    [Waples_Ne_AFT_ji_norm_CI(:,2)] = Waples_Ne_estimate(har_mean_S_delta,r_sqr_d_AFT_jack_individ_CI_norm(:,2)); % col vector of upper CI Waples Ne from r_sqr_delta_AFW_jack_individ normal

    % UB_factor(iter,2) is unbiased factor for r_squared_delta
    [Waples_Ne_AFW_UB_ji_pct_CI(:,1)] = Waples_Ne_estimate(har_mean_S_delta,r_sqr_d_AFW_jack_individ_CI_pct(:,1) .* UB_factor(:,2)); % col vector of lower CI Waples Ne from r_sqr_delta_AFW_jack_individ times UB_factor percentile
    [Waples_Ne_AFW_UB_ji_pct_CI(:,2)] = Waples_Ne_estimate(har_mean_S_delta,r_sqr_d_AFW_jack_individ_CI_pct(:,2) .* UB_factor(:,2)); % col vector of upper CI Waples Ne from r_sqr_delta_AFW_jack_individ times UB_factor percentile

    [Waples_Ne_AFT_UB_ji_pct_CI(:,1)] = Waples_Ne_estimate(har_mean_S_delta,r_sqr_d_AFT_jack_individ_CI_pct(:,1) .* UB_factor(:,2)); % col vector of lower CI Waples Ne from r_sqr_delta_AFW_jack_individ times UB_factor percentile
    [Waples_Ne_AFT_UB_ji_pct_CI(:,2)] = Waples_Ne_estimate(har_mean_S_delta,r_sqr_d_AFT_jack_individ_CI_pct(:,2) .* UB_factor(:,2)); % col vector of upper CI Waples Ne from r_sqr_delta_AFW_jack_individ times UB_factor percentile

    [Waples_Ne_AFW_UB_ji_norm_CI(:,1)] = Waples_Ne_estimate(har_mean_S_delta,r_sqr_d_AFW_jack_individ_CI_norm(:,1) .* UB_factor(:,2)); % col vector of lower CI Waples Ne from r_sqr_delta_AFW_jack_individ times UB_factor normal
    [Waples_Ne_AFW_UB_ji_norm_CI(:,2)] = Waples_Ne_estimate(har_mean_S_delta,r_sqr_d_AFW_jack_individ_CI_norm(:,2) .* UB_factor(:,2)); % col vector of upper CI Waples Ne from r_sqr_delta_AFW_jack_individ times UB_factor normal
    
    [Waples_Ne_AFT_UB_ji_norm_CI(:,1)] = Waples_Ne_estimate(har_mean_S_delta,r_sqr_d_AFT_jack_individ_CI_norm(:,1) .* UB_factor(:,2)); % col vector of lower CI Waples Ne from r_sqr_delta_AFW_jack_individ times UB_factor normal
    [Waples_Ne_AFT_UB_ji_norm_CI(:,2)] = Waples_Ne_estimate(har_mean_S_delta,r_sqr_d_AFT_jack_individ_CI_norm(:,2) .* UB_factor(:,2)); % col vector of upper CI Waples Ne from r_sqr_delta_AFW_jack_individ times UB_factor normal
    
    
    % compute Waples Ne with jackknife over loci r_sqr_delta upper and lower CIs
    % first column contains upper bound and second column contains lower bound
    [Waples_Ne_AFW_jl_pct_CI(:,1)] = Waples_Ne_estimate(har_mean_S_delta,r_sqr_d_AFW_jack_loci_CI_pct(:,1)); % col vector of lower CI Waples Ne from r_sqr_delta_AFW_jack_loci percentile
    [Waples_Ne_AFW_jl_pct_CI(:,2)] = Waples_Ne_estimate(har_mean_S_delta,r_sqr_d_AFW_jack_loci_CI_pct(:,2)); % col vector of upper CI Waples Ne from r_sqr_delta_AFW_jack_loci percentile

    [Waples_Ne_AFT_jl_pct_CI(:,1)] = Waples_Ne_estimate(har_mean_S_delta,r_sqr_d_AFT_jack_loci_CI_pct(:,1)); % col vector of lower CI Waples Ne from r_sqr_delta_AFW_jack_loci percentile
    [Waples_Ne_AFT_jl_pct_CI(:,2)] = Waples_Ne_estimate(har_mean_S_delta,r_sqr_d_AFT_jack_loci_CI_pct(:,2)); % col vector of upper CI Waples Ne from r_sqr_delta_AFW_jack_loci percentile

    [Waples_Ne_AFW_jl_norm_CI(:,1)] = Waples_Ne_estimate(har_mean_S_delta,r_sqr_d_AFW_jack_loci_CI_norm(:,1)); % col vector of lower CI Waples Ne from r_sqr_delta_AFW_jack_loci normal
    [Waples_Ne_AFW_jl_norm_CI(:,2)] = Waples_Ne_estimate(har_mean_S_delta,r_sqr_d_AFW_jack_loci_CI_norm(:,2)); % col vector of upper CI Waples Ne from r_sqr_delta_AFW_jack_loci normal
    
    [Waples_Ne_AFT_jl_norm_CI(:,1)] = Waples_Ne_estimate(har_mean_S_delta,r_sqr_d_AFT_jack_loci_CI_norm(:,1)); % col vector of lower CI Waples Ne from r_sqr_delta_AFW_jack_loci normal
    [Waples_Ne_AFT_jl_norm_CI(:,2)] = Waples_Ne_estimate(har_mean_S_delta,r_sqr_d_AFT_jack_loci_CI_norm(:,2)); % col vector of upper CI Waples Ne from r_sqr_delta_AFW_jack_loci normal
    
    % UB_factor(iter,2) is unbiased factor for r_squared_delta
    [Waples_Ne_AFW_UB_jl_pct_CI(:,1)] = Waples_Ne_estimate(har_mean_S_delta,r_sqr_d_AFW_jack_loci_CI_pct(:,1) .* UB_factor(:,2)); % col vector of lower CI Waples Ne from r_sqr_delta_AFW_jack_individ times UB_factor percentile
    [Waples_Ne_AFW_UB_jl_pct_CI(:,2)] = Waples_Ne_estimate(har_mean_S_delta,r_sqr_d_AFW_jack_loci_CI_pct(:,2) .* UB_factor(:,2)); % col vector of upper CI Waples Ne from r_sqr_delta_AFW_jack_individ times UB_factor percentile

    [Waples_Ne_AFT_UB_jl_pct_CI(:,1)] = Waples_Ne_estimate(har_mean_S_delta,r_sqr_d_AFT_jack_loci_CI_pct(:,1) .* UB_factor(:,2)); % col vector of lower CI Waples Ne from r_sqr_delta_AFW_jack_individ times UB_factor percentile
    [Waples_Ne_AFT_UB_jl_pct_CI(:,2)] = Waples_Ne_estimate(har_mean_S_delta,r_sqr_d_AFT_jack_loci_CI_pct(:,2) .* UB_factor(:,2)); % col vector of upper CI Waples Ne from r_sqr_delta_AFW_jack_individ times UB_factor percentile

    [Waples_Ne_AFW_UB_jl_norm_CI(:,1)] = Waples_Ne_estimate(har_mean_S_delta,r_sqr_d_AFW_jack_loci_CI_norm(:,1) .* UB_factor(:,2)); % col vector of lower CI Waples Ne from r_sqr_delta_AFW_jack_individ times UB_factor normal
    [Waples_Ne_AFW_UB_jl_norm_CI(:,2)] = Waples_Ne_estimate(har_mean_S_delta,r_sqr_d_AFW_jack_loci_CI_norm(:,2) .* UB_factor(:,2)); % col vector of upper CI Waples Ne from r_sqr_delta_AFW_jack_individ times UB_factor normal
    
    [Waples_Ne_AFT_UB_jl_norm_CI(:,1)] = Waples_Ne_estimate(har_mean_S_delta,r_sqr_d_AFT_jack_loci_CI_norm(:,1) .* UB_factor(:,2)); % col vector of lower CI Waples Ne from r_sqr_delta_AFW_jack_individ times UB_factor normal
    [Waples_Ne_AFT_UB_jl_norm_CI(:,2)] = Waples_Ne_estimate(har_mean_S_delta,r_sqr_d_AFT_jack_loci_CI_norm(:,2) .* UB_factor(:,2)); % col vector of upper CI Waples Ne from r_sqr_delta_AFW_jack_individ times UB_factor normal


    % determine how many times Ne estimates fall inside the various jackknife CIs
    vec9_ji = Ne(:,9) >= Waples_Ne_AFW_ji_pct_CI(:,2) & Ne(:,9) <= Waples_Ne_AFW_ji_pct_CI(:,1);
    proportion9_ji_pct = sum(vec9_ji)/num_sims;

    vec10_ji = Ne(:,10) >= Waples_Ne_AFT_ji_pct_CI(:,2) & Ne(:,10) <= Waples_Ne_AFT_ji_pct_CI(:,1);
    proportion10_ji_pct = sum(vec10_ji)/num_sims;
    
    vec11_ji = Ne(:,9) >= Waples_Ne_AFW_ji_norm_CI(:,2) & Ne(:,9) <= Waples_Ne_AFW_ji_norm_CI(:,1);
    proportion11_ji_norm = sum(vec11_ji)/num_sims;

    vec12_ji = Ne(:,10) >= Waples_Ne_AFT_ji_norm_CI(:,2) & Ne(:,10) <= Waples_Ne_AFT_ji_norm_CI(:,1);
    proportion12_ji_norm = sum(vec12_ji)/num_sims;
    
    vec9_ji_UB = Ne(:,11) >= Waples_Ne_AFW_UB_ji_pct_CI(:,2) & Ne(:,11) <= Waples_Ne_AFW_UB_ji_pct_CI(:,1);
    proportion9_ji_pct_UB = sum(vec9_ji_UB)/num_sims;

    vec10_ji_UB = Ne(:,12) >= Waples_Ne_AFT_UB_ji_pct_CI(:,2) & Ne(:,12) <= Waples_Ne_AFT_UB_ji_pct_CI(:,1);
    proportion10_ji_pct_UB = sum(vec10_ji_UB)/num_sims;
    
    vec11_ji_UB = Ne(:,11) >= Waples_Ne_AFW_UB_ji_norm_CI(:,2) & Ne(:,11) <= Waples_Ne_AFW_UB_ji_norm_CI(:,1);
    proportion11_ji_norm_UB = sum(vec11_ji_UB)/num_sims;

    vec12_ji_UB = Ne(:,12) >= Waples_Ne_AFT_UB_ji_norm_CI(:,2) & Ne(:,12) <= Waples_Ne_AFT_UB_ji_norm_CI(:,1);
    proportion12_ji_norm_UB = sum(vec12_ji_UB)/num_sims;
    

    vec9_jl = Ne(:,9) >= Waples_Ne_AFW_jl_pct_CI(:,2) & Ne(:,9) <= Waples_Ne_AFW_jl_pct_CI(:,1);
    proportion9_jl_pct = sum(vec9_jl)/num_sims;

    vec10_jl = Ne(:,10) >= Waples_Ne_AFT_jl_pct_CI(:,2) & Ne(:,10) <= Waples_Ne_AFT_jl_pct_CI(:,1);
    proportion10_jl_pct = sum(vec10_jl)/num_sims;
    
    vec11_jl = Ne(:,9) >= Waples_Ne_AFW_jl_norm_CI(:,2) & Ne(:,9) <= Waples_Ne_AFW_jl_norm_CI(:,1);
    proportion11_jl_norm = sum(vec11_jl)/num_sims;

    vec12_jl = Ne(:,10) >= Waples_Ne_AFT_jl_norm_CI(:,2) & Ne(:,10) <= Waples_Ne_AFT_jl_norm_CI(:,1);
    proportion12_jl_norm = sum(vec12_jl)/num_sims;
    
    vec9_jl_UB = Ne(:,11) >= Waples_Ne_AFW_UB_jl_pct_CI(:,2) & Ne(:,11) <= Waples_Ne_AFW_UB_jl_pct_CI(:,1);
    proportion9_jl_pct_UB = sum(vec9_jl_UB)/num_sims;

    vec10_jl_UB = Ne(:,12) >= Waples_Ne_AFT_UB_jl_pct_CI(:,2) & Ne(:,12) <= Waples_Ne_AFT_UB_jl_pct_CI(:,1);
    proportion10_jl_pct_UB = sum(vec10_jl_UB)/num_sims;
    
    vec11_jl_UB = Ne(:,11) >= Waples_Ne_AFW_UB_jl_norm_CI(:,2) & Ne(:,11) <= Waples_Ne_AFW_UB_jl_norm_CI(:,1);
    proportion11_jl_norm_UB = sum(vec11_jl_UB)/num_sims;

    vec12_jl_UB = Ne(:,12) >= Waples_Ne_AFT_UB_jl_norm_CI(:,2) & Ne(:,12) <= Waples_Ne_AFT_UB_jl_norm_CI(:,1);
    proportion12_jl_norm_UB = sum(vec12_jl_UB)/num_sims;
    
    % table of output
    fprintf('\n');
    fprintf('____________________________________________________________________\n');
    fprintf('Proportion of times Ne estimate falls inside of estimated CI with .\n');
    fprintf('jackknife over loci and jackknife over individuals.\n');
    fprintf('Ne estimator is Waples regression equation using r_sqr_delta with AFW or AFT\n');
    fprintf('and with or without multiplations by UB factor.\n');
    fprintf('True Ne = %d\n', pop_size);
    fprintf('Number of independent simulated data sets = %d\n', num_sims);
    fprintf('____________________________________________________________________\n');
    fprintf('Estimator\n');
    fprintf('\n');
    fprintf('jackknife over individuals with percentile CI:\n');
    fprintf('r_sqr_delta_AFW\t\t%f\n', proportion9_ji_pct);
    fprintf('r_sqr_delta_AFT\t\t%f\n', proportion10_ji_pct);
    fprintf('r_sqr_delta_AFW UB\t%f\n', proportion9_ji_pct_UB);
    fprintf('r_sqr_delta_AFT UB\t%f\n', proportion9_ji_pct_UB);
    
    fprintf('\n');
    
    fprintf('jackknife over individuals with normal CI:\n');
    fprintf('r_sqr_delta_AFW\t\t%f\n', proportion11_ji_norm);
    fprintf('r_sqr_delta_AFT\t\t%f\n', proportion12_ji_norm);
    fprintf('r_sqr_delta_AFW UB\t%f\n', proportion11_ji_norm_UB);
    fprintf('r_sqr_delta_AFT UB\t%f\n', proportion12_ji_norm_UB);
    
    fprintf('\n');
    
    fprintf('jackknife over loci with percentile CI:\n');
    fprintf('r_sqr_delta_AFW\t\t%f\n', proportion9_jl_pct);
    fprintf('r_sqr_delta_AFT\t\t%f\n', proportion10_jl_pct);    
    fprintf('r_sqr_delta_AFW UB\t%f\n', proportion9_jl_pct_UB);
    fprintf('r_sqr_delta_AFT UB\t%f\n', proportion10_jl_pct_UB);
    
    fprintf('\n');
    
    fprintf('jackknife over loci with normal CI:\n');
    fprintf('r_sqr_delta_AFW\t\t%f\n', proportion11_jl_norm);
    fprintf('r_sqr_delta_AFT\t\t%f\n', proportion12_jl_norm);    
    fprintf('r_sqr_delta_AFW UB\t%f\n', proportion11_jl_norm_UB);
    fprintf('r_sqr_delta_AFT UB\t%f\n', proportion12_jl_norm_UB);
    fprintf('\n');

    fprintf('____________________________________________________________________\n');
    fprintf('\n');



%**********
    
    % determine how many times Ne estimates fall inside the various jackknife CIs
    % jackknife over loci 
    
    % percentile CI    
    % Ne estimated with r_squared_d AFW and r_squared_d AFT
    vec3 = Ne(:,5) >= r_sqr_d_AFW_Ne_jack_loci_CI_pct(:,1) & Ne(:,5) <= r_sqr_d_AFW_Ne_jack_loci_CI_pct(:,2);
    proportion3 = sum(vec3)/num_sims;
    
    vec4 = Ne(:,6) >= r_sqr_d_AFT_Ne_jack_loci_CI_pct(:,1) & Ne(:,6) <= r_sqr_d_AFT_Ne_jack_loci_CI_pct(:,2);
    proportion4 = sum(vec4)/num_sims;
    
    vec3_UB = Ne(:,7) >= r_sqr_d_AFW_Ne_jack_loci_CI_pct(:,1) & Ne(:,7) <= r_sqr_d_AFW_Ne_jack_loci_CI_pct(:,2);
    proportion3_UB = sum(vec3_UB)/num_sims;
    
    vec4_UB = Ne(:,8) >= r_sqr_d_AFT_Ne_jack_loci_CI_pct(:,1) & Ne(:,8) <= r_sqr_d_AFT_Ne_jack_loci_CI_pct(:,2);
    proportion4_UB = sum(vec4_UB)/num_sims;   

    % normal CI    
    % Ne estimated with r_squared_d AFW and r_squared_d AFT
    vec7 = Ne(:,5) >= r_sqr_d_AFW_Ne_jack_loci_CI_norm(:,1) & Ne(:,5) <= r_sqr_d_AFW_Ne_jack_loci_CI_norm(:,2);
    proportion7 = sum(vec7)/num_sims;
    
    vec8 = Ne(:,6) >= r_sqr_d_AFT_Ne_jack_loci_CI_norm(:,1) & Ne(:,6) <= r_sqr_d_AFT_Ne_jack_loci_CI_norm(:,2);
    proportion8 = sum(vec8)/num_sims;
    
    vec7_UB = Ne(:,7) >= r_sqr_d_AFW_Ne_jack_loci_CI_norm(:,1) & Ne(:,7) <= r_sqr_d_AFW_Ne_jack_loci_CI_norm(:,2);
    proportion7_UB = sum(vec7_UB)/num_sims;
    
    vec8_UB = Ne(:,8) >= r_sqr_d_AFT_Ne_jack_loci_CI_norm(:,1) & Ne(:,8) <= r_sqr_d_AFT_Ne_jack_loci_CI_norm(:,2);
    proportion8_UB = sum(vec8_UB)/num_sims;   

    % jackknife over loci with percentile CI
    % Ne estimated with r_squared_c AFW permute and r_squared_c AFT permute
    vec13_pct = Ne(:,13) >= r_sqr_c_AFW_Ne_jack_loci_CI_pct(:,1) & Ne(:,13) <= r_sqr_c_AFW_Ne_jack_loci_CI_pct(:,2);
    proportion13_pct = sum(vec13_pct)/num_sims;
    
    vec14_pct = Ne(:,14) >= r_sqr_c_AFT_Ne_jack_loci_CI_pct(:,1) & Ne(:,14) <= r_sqr_c_AFT_Ne_jack_loci_CI_pct(:,2);
    proportion14_pct = sum(vec14_pct)/num_sims;
    
    vec13_UB_pct = Ne(:,15) >= r_sqr_c_AFW_Ne_jack_loci_CI_pct(:,1) & Ne(:,15) <= r_sqr_c_AFW_Ne_jack_loci_CI_pct(:,2);
    proportion13_UB_pct = sum(vec13_UB_pct)/num_sims;
    
    vec14_UB_pct = Ne(:,16) >= r_sqr_c_AFT_Ne_jack_loci_CI_pct(:,1) & Ne(:,16) <= r_sqr_c_AFT_Ne_jack_loci_CI_pct(:,2);
    proportion14_UB_pct = sum(vec14_UB_pct)/num_sims;   
    
    % jackknife over loci with norm CI
    % Ne estimated with r_squared_c AFW permute and r_squared_c AFT permute
    vec13_norm = Ne(:,13) >= r_sqr_c_AFW_Ne_jack_loci_CI_norm(:,1) & Ne(:,13) <= r_sqr_c_AFW_Ne_jack_loci_CI_norm(:,2);
    proportion13_norm = sum(vec13_norm)/num_sims;
    
    vec14_norm = Ne(:,14) >= r_sqr_c_AFT_Ne_jack_loci_CI_norm(:,1) & Ne(:,14) <= r_sqr_c_AFT_Ne_jack_loci_CI_norm(:,2);
    proportion14_norm = sum(vec14_norm)/num_sims;
    
    vec13_UB_norm = Ne(:,15) >= r_sqr_c_AFW_Ne_jack_loci_CI_norm(:,1) & Ne(:,15) <= r_sqr_c_AFW_Ne_jack_loci_CI_norm(:,2);
    proportion13_UB_norm = sum(vec13_UB_norm)/num_sims;
    
    vec14_UB_norm = Ne(:,16) >= r_sqr_c_AFT_Ne_jack_loci_CI_norm(:,1) & Ne(:,16) <= r_sqr_c_AFT_Ne_jack_loci_CI_norm(:,2);
    proportion14_UB_norm = sum(vec14_UB_norm)/num_sims;   


    % table of output
    fprintf('\n');
    fprintf('____________________________________________________________________\n');
    fprintf('Proportion of times Ne estimate falls inside of estimated CI\n');
    fprintf('for jackknife over loci. \n');
    fprintf('True Ne = %d\n', pop_size);
    fprintf('Number of independent simulated data sets = %d\n', num_sims);
    fprintf('____________________________________________________________________\n');
    fprintf('Estimator\n');
    fprintf('\n');
    fprintf('jackknife over loci with percentile CI:\n');
    fprintf('r_sqr_delta_AFW\t\t%f\n', proportion3);
    fprintf('r_sqr_delta_AFT\t\t%f\n', proportion4);
    fprintf('r_sqr_delta_AFW UB\t%f\n', proportion3_UB);
    fprintf('r_sqr_delta_AFT UB\t%f\n', proportion4_UB);
    fprintf('\n');

    fprintf('jackknife over loci with normal CI:\n');
    fprintf('r_sqr_delta_AFW\t\t%f\n', proportion7);
    fprintf('r_sqr_delta_AFT\t\t%f\n', proportion8);
    fprintf('r_sqr_delta_AFW UB\t%f\n', proportion7_UB);
    fprintf('r_sqr_delta_AFT UB\t%f\n', proportion8_UB);
    fprintf('\n');     
   
    fprintf('jackknife over loci with percentile CI:\n');
    fprintf('r_sqr_c_AFW permute\t%f\n', proportion13_pct);
    fprintf('r_sqr_c_AFT permute\t%f\n', proportion14_pct);    
    fprintf('r_sqr_c_AFW UB permute\t%f\n', proportion13_UB_pct);
    fprintf('r_sqr_c_AFT UB permute\t%f\n', proportion14_UB_pct);
    fprintf('\n');

    fprintf('jackknife over loci with normal CI:\n');
    fprintf('r_sqr_c_AFW permute\t%f\n', proportion13_norm);
    fprintf('r_sqr_c_AFT permute\t%f\n', proportion14_norm);    
    fprintf('r_sqr_c_AFW UB permute\t%f\n', proportion13_UB_norm);
    fprintf('r_sqr_c_AFT UB permute\t%f\n', proportion14_UB_norm);
    fprintf('\n');    
    fprintf('____________________________________________________________________\n');
    fprintf('\n');


    % determine how many times Ne estimates fall inside the various jackknife CIs
    % jackknife over individuals 
    
    % percentile CI
    % Ne estimated with r_squared_d AFW and r_squared_d AFT
    vec3 = Ne(:,5) >= r_sqr_d_AFW_Ne_jack_loci_CI_pct(:,1) & Ne(:,5) <= r_sqr_d_AFW_Ne_jack_individ_CI_pct(:,2);
    proportion3 = sum(vec3)/num_sims;
    
    vec4 = Ne(:,6) >= r_sqr_d_AFT_Ne_jack_individ_CI_pct(:,1) & Ne(:,6) <= r_sqr_d_AFT_Ne_jack_individ_CI_pct(:,2);
    proportion4 = sum(vec4)/num_sims;
    
    vec3_UB = Ne(:,7) >= r_sqr_d_AFW_Ne_jack_individ_CI_pct(:,1) & Ne(:,7) <= r_sqr_d_AFW_Ne_jack_individ_CI_pct(:,2);
    proportion3_UB = sum(vec3_UB)/num_sims;
    
    vec4_UB = Ne(:,8) >= r_sqr_d_AFT_Ne_jack_individ_CI_pct(:,1) & Ne(:,8) <= r_sqr_d_AFT_Ne_jack_individ_CI_pct(:,2);
    proportion4_UB = sum(vec4_UB)/num_sims;   

    % normal CI
    % Ne estimated with r_squared_d AFW and r_squared_d AFT
    vec7 = Ne(:,5) >= r_sqr_d_AFW_Ne_jack_individ_CI_norm(:,1) & Ne(:,5) <= r_sqr_d_AFW_Ne_jack_individ_CI_norm(:,2);
    proportion7 = sum(vec7)/num_sims;
    
    vec8 = Ne(:,6) >= r_sqr_d_AFT_Ne_jack_individ_CI_norm(:,1) & Ne(:,6) <= r_sqr_d_AFT_Ne_jack_individ_CI_norm(:,2);
    proportion8 = sum(vec8)/num_sims;
    
    vec7_UB = Ne(:,7) >= r_sqr_d_AFW_Ne_jack_individ_CI_norm(:,1) & Ne(:,7) <= r_sqr_d_AFW_Ne_jack_individ_CI_norm(:,2);
    proportion7_UB = sum(vec7_UB)/num_sims;
    
    vec8_UB = Ne(:,8) >= r_sqr_d_AFT_Ne_jack_individ_CI_norm(:,1) & Ne(:,8) <= r_sqr_d_AFT_Ne_jack_individ_CI_norm(:,2);
    proportion8_UB = sum(vec8_UB)/num_sims;   

    


%     % percentile CI
%     % Ne estimated with Waples regression equations fron r_squared_delta AFW and r_squared_delta AFT 
%     vec9_pct = Ne(:,9) >= Waples_Ne_AFW_ji_pct_CI(:,2) & Ne(:,9) <= Waples_Ne_AFW_ji_pct_CI(:,1);
%     proportion9_pct = sum(vec9_pct)/num_sims;
%     
%     vec10_pct = Ne(:,10) >= Waples_Ne_AFT_ji_pct_CI(:,2) & Ne(:,10) <= Waples_Ne_AFT_ji_pct_CI(:,1);
%     proportion10_pct = sum(vec10_pct)/num_sims;
%     
%     vec13_UB_pct = Ne(:,15) >= Waples_Ne_AFW_jl_pct_CI(:,1) & Ne(:,15) <= Waples_Ne_AFW_jl_pct_CI(:,1);
%     proportion13_UB_pct = sum(vec13_UB_pct)/num_sims;
%     
%     vec14_UB_pct = Ne(:,16) >= Waples_Ne_AFT_jl_pct_CI(:,2) & Ne(:,16) <= Waples_Ne_AFT_jl_pct_CI(:,1);
%     proportion14_UB_pct = sum(vec14_UB_pct)/num_sims;   
%     
%     % normal CI
%     % Ne estimated with r_squared_c AFW permute and r_squared_c AFT permute
%     vec13_norm = Ne(:,13) >= r_sqr_c_AFW_Ne_jack_individ_CI_norm(:,1) & Ne(:,13) <= r_sqr_c_AFW_Ne_jack_individ_CI_norm(:,2);
%     proportion13_norm = sum(vec13_norm)/num_sims;
%     
%     vec14_norm = Ne(:,14) >= r_sqr_c_AFT_Ne_jack_individ_CI_norm(:,1) & Ne(:,14) <= r_sqr_c_AFT_Ne_jack_individ_CI_norm(:,2);
%     proportion14_norm = sum(vec14_norm)/num_sims;
%     
%     vec13_UB_norm = Ne(:,15) >= r_sqr_c_AFW_Ne_jack_individ_CI_norm(:,1) & Ne(:,15) <= r_sqr_c_AFW_Ne_jack_individ_CI_norm(:,2);
%     proportion13_UB_norm = sum(vec13_UB_norm)/num_sims;
%     
%     vec14_UB_norm = Ne(:,16) >= r_sqr_c_AFT_Ne_jack_individ_CI_norm(:,1) & Ne(:,16) <= r_sqr_c_AFT_Ne_jack_individ_CI_norm(:,2);
%     proportion14_UB_norm = sum(vec14_UB_norm)/num_sims;   


    % percentile CI
    % Ne estimated with r_squared_c AFW permute and r_squared_c AFT permute
    vec13_pct = Ne(:,13) >= r_sqr_c_AFW_Ne_jack_individ_CI_pct(:,1) & Ne(:,13) <= r_sqr_c_AFW_Ne_jack_individ_CI_pct(:,2);
    proportion13_pct = sum(vec13_pct)/num_sims;
    
    vec14_pct = Ne(:,14) >= r_sqr_c_AFT_Ne_jack_individ_CI_pct(:,1) & Ne(:,14) <= r_sqr_c_AFT_Ne_jack_individ_CI_pct(:,2);
    proportion14_pct = sum(vec14_pct)/num_sims;
    
    vec13_UB_pct = Ne(:,15) >= r_sqr_c_AFW_Ne_jack_individ_CI_pct(:,1) & Ne(:,15) <= r_sqr_c_AFW_Ne_jack_individ_CI_pct(:,2);
    proportion13_UB_pct = sum(vec13_UB_pct)/num_sims;
    
    vec14_UB_pct = Ne(:,16) >= r_sqr_c_AFT_Ne_jack_individ_CI_pct(:,1) & Ne(:,16) <= r_sqr_c_AFT_Ne_jack_individ_CI_pct(:,2);
    proportion14_UB_pct = sum(vec14_UB_pct)/num_sims;   
    
    % normal CI
    % Ne estimated with r_squared_c AFW permute and r_squared_c AFT permute
    vec13_norm = Ne(:,13) >= r_sqr_c_AFW_Ne_jack_individ_CI_norm(:,1) & Ne(:,13) <= r_sqr_c_AFW_Ne_jack_individ_CI_norm(:,2);
    proportion13_norm = sum(vec13_norm)/num_sims;
    
    vec14_norm = Ne(:,14) >= r_sqr_c_AFT_Ne_jack_individ_CI_norm(:,1) & Ne(:,14) <= r_sqr_c_AFT_Ne_jack_individ_CI_norm(:,2);
    proportion14_norm = sum(vec14_norm)/num_sims;
    
    vec13_UB_norm = Ne(:,15) >= r_sqr_c_AFW_Ne_jack_individ_CI_norm(:,1) & Ne(:,15) <= r_sqr_c_AFW_Ne_jack_individ_CI_norm(:,2);
    proportion13_UB_norm = sum(vec13_UB_norm)/num_sims;
    
    vec14_UB_norm = Ne(:,16) >= r_sqr_c_AFT_Ne_jack_individ_CI_norm(:,1) & Ne(:,16) <= r_sqr_c_AFT_Ne_jack_individ_CI_norm(:,2);
    proportion14_UB_norm = sum(vec14_UB_norm)/num_sims;   

    % table of output
    fprintf('\n');
    fprintf('____________________________________________________________________\n');
    fprintf('Proportion of times Ne estimate falls inside of estimated CI\n');
    fprintf('for jackknife over individuals. \n');
    fprintf('True Ne = %d\n', pop_size);
    fprintf('Number of independent simulated data sets = %d\n', num_sims);
    fprintf('____________________________________________________________________\n');
    fprintf('Estimator\n');
    fprintf('\n');
    fprintf('jackknife over individuals with percentile CI:\n');
    fprintf('r_sqr_delta_AFW\t\t%f\n', proportion3);
    fprintf('r_sqr_delta_AFT\t\t%f\n', proportion4);
    fprintf('r_sqr_delta_AFW UB\t%f\n', proportion3_UB);
    fprintf('r_sqr_delta_AFT UB\t%f\n', proportion4_UB);
    
    fprintf('\n');
    
    fprintf('jackknife over individuals with normal CI:\n');
    fprintf('r_sqr_delta_AFW\t%f\n', proportion7);
    fprintf('r_sqr_delta_AFT\t%f\n', proportion8);
    fprintf('r_sqr_delta_AFW UB\t%f\n', proportion7_UB);
    fprintf('r_sqr_delta_AFT UB\t%f\n', proportion8_UB);
    
    fprintf('\n\n');
    
    fprintf('jackknife over individuals with percentile CI:\n');
    fprintf('r_sqr_c_AFW permute\t%f\n', proportion13_pct);
    fprintf('r_sqr_c_AFT permute\t%f\n', proportion14_pct);    
    fprintf('r_sqr_c_AFW UB permute\t%f\n', proportion13_UB_pct);
    fprintf('r_sqr_c_AFT UB permute\t%f\n', proportion14_UB_pct);
    
    fprintf('\n');
    
    fprintf('jackknife over individuals with normal CI:\n');
    fprintf('r_sqr_c_AFW permute\t%f\n', proportion13_norm);
    fprintf('r_sqr_c_AFT permute\t%f\n', proportion14_norm);    
    fprintf('r_sqr_c_AFW UB permute\t%f\n', proportion13_UB_norm);
    fprintf('r_sqr_c_AFT UB permute\t%f\n', proportion14_UB_norm);
    fprintf('\n');

    fprintf('____________________________________________________________________\n');
    fprintf('\n');


%**********
    
%     % compare variance in Ne estimates from all simulated data sets with
%     % mean jackknife variances 
%     
%     % compute variances from Ne estimates from all simulated data sets    
%     var_Ne_r_squared_delta_AFW = var(Ne(:,5));
%     var_Ne_r_squared_delta_AFT = var(Ne(:,6));
%     var_Ne_r_squared_delta_AFW_UB = var(Ne(:,7));
%     var_Ne_r_squared_delta_AFT_UB = var(Ne(:,8));    
%     
%     var_Ne_Waples_r_squared_delta_AFW = var(Ne(:,9));
%     var_Ne_Waples_r_squared_delta_AFT = var(Ne(:,10));
%     var_Ne_Waples_r_squared_delta_AFW_UB = var(Ne(:,11));
%     var_Ne_Waples_r_squared_delta_AFT_UB = var(Ne(:,12));
% 
%     var_Ne_r_squared_c_AFW_perm = var(Ne(:,13));
%     var_Ne_r_squared_c_AFT_perm = var(Ne(:,14));
%     var_Ne_r_squared_c_AFW_UB_perm = var(Ne(:,15));
%     var_Ne_r_squared_c_AFT_UB_perm = var(Ne(:,16));    
    
    


%**********
            
    % compute variances from r^2 estimates from all simulated data sets
    var_r_squared_c_AFW = var(r_squared_c(:,1));
    var_r_squared_c_AFT = var(r_squared_c(:,2));
    var_r_squared_d_AFW = var(r_squared_delta(:,1));
    var_r_squared_d_AFT = var(r_squared_delta(:,2));

    % *** jackknife over individuals
    % compute average jackknife estimated variances from all simulated data sets
    avg_j_individ_var_r_sqr_c_AFW = mean(jack_individ_var_r_sqr_c_AFW);
    avg_j_individ_var_r_sqr_c_AFT = mean(jack_individ_var_r_sqr_c_AFT);
    avg_j_individ_var_r_sqr_delta_AFW = mean(jack_individ_var_r_sqr_delta_AFW);
    avg_j_individ_var_r_sqr_delta_AFT = mean(jack_individ_var_r_sqr_delta_AFT);
    
    % *** jackknife over loci    
    % compute average jackknife estimated variances from all simulated data sets
    avg_j_loci_var_r_sqr_c_AFW = mean(jack_loci_var_r_sqr_c_AFW);
    avg_j_loci_var_r_sqr_c_AFT = mean(jack_loci_var_r_sqr_c_AFT);
    avg_j_loci_var_r_sqr_delta_AFW = mean(jack_loci_var_r_sqr_delta_AFW);
    avg_j_loci_var_r_sqr_delta_AFT = mean(jack_loci_var_r_sqr_delta_AFT);
    
    fprintf('\n');
    fprintf('____________________________________________________________________\n');
    fprintf('Variances in r^2 from independent simulations and from the average of variances\n');
    fprintf('estimated with jackknife over individuals or loci for each data set.\n');
    fprintf('True Ne = %d\n', pop_size);
    fprintf('Number of independent simulated data sets = %d\n', num_sims);
    fprintf('_________________________________________________________________________________________\n');
    fprintf('estimator\tvar over independent\tmean variance\t\tmean variance\n');
    fprintf('\t\tdata replicates\t\tjackknife individuals\tjackknife loci\n');
    fprintf('r_sqr_c_AFW\t%f\t\t%f\t\t%f\n', var_r_squared_c_AFW, avg_j_individ_var_r_sqr_c_AFW, avg_j_loci_var_r_sqr_c_AFW);
    fprintf('r_sqr_c_AFT\t%f\t\t%f\t\t%f\n', var_r_squared_c_AFT, avg_j_individ_var_r_sqr_c_AFT, avg_j_loci_var_r_sqr_c_AFT);
    fprintf('\n');
    fprintf('r_sqr_delta_AFW\t%f\t\t%f\t\t%f\n', var_r_squared_d_AFW, avg_j_individ_var_r_sqr_delta_AFW, avg_j_loci_var_r_sqr_delta_AFW);
    fprintf('r_sqr_delta_AFT\t%f\t\t%f\t\t%f\n', var_r_squared_d_AFT, avg_j_individ_var_r_sqr_delta_AFT, avg_j_loci_var_r_sqr_delta_AFT);
    fprintf('_________________________________________________________________________________________\n');
    fprintf('\n');

    
