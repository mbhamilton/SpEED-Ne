function simulate_r
% Program designed as a simulator that calls functions to estimate
% r_squared_comp and r_squared_delta using simulated data. 
% Utilizes the same underlying estimation functions as the GUI version of SpEED-Ne.
%
% This version simulates genotype data for each iteration using generate_pop.m. 
%
%   version 1.6  03 January 2017
%
% Makes 16 different estimates of Ne:
%
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
%
% where:
%   r_squared_c = r^2 composite
%	r_squared_delta = r^2 delta method to account for within locus disequilibrium
%	AFT = allele freq thresholded
%	AFW = allele freq weighted
%	UB = unbiased becuase of multiplation by S/(S-1)
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

    rng('shuffle'); %initialize random number generator with a seed based on the current time
    
    %set key parameter values 
    output_filename = '/Users/matthewhamilton/Dropbox/R project/manuscript/time forward Ne simulation results/03_Jan_no_selfing_Ne_'; % file name for saved workspace
    file_name_ext = '.mat'; % file name extension
    
    pop_size = 50; % true Ne value
    
    num_sims = 250;   % number of simulation iterations for a given true Ne
    
    threshold_allele_freq = 0.05;
    n_permutations = 2500; % number of permutations of genotype data to estimate null value of r_squared
    alpha = 0.05;   % total width of confidence interval tails in jackknife and bootstrap
    
    % switches for output, supresses output during iterative simulations

    % permute genotype data for correction factor
    permute = 1; % 1=yes, 0=no
    
    % output table of r^2 estimates for all pairs of loci
    print_table = 0; % 1=yes, 0=no
    
    % print table of allele states and frequencies
    print_allele_table = 0; % 1=yes, 0=no
    
    % print table of r^2 estimates for each pair of alleles within all loci
    print_locus_table = 0; % 1=yes, 0=no
    
    % output of graphs
    make_graphs = 0; % 1=yes, 0=no
    
    print_notes = false;  % print error output strings to console
    
    make_files = false;  % output files of estimates from all locus pairs and all allele pairs


    % parameters for simulated populations and genotype data
    mating_model = 1; % random mating w/ possible selfing =1, random without selfing = 2
    num_loci = 12;
    num_allele_locus = 4;
    gens_2_sim = 7;
    
    
    % allocate space
    num_locus_pairs = (num_loci*(num_loci-1))/2;
    
    r_squared_c = zeros(num_sims,2); % 1st column is allele frequency weighted (AFW), 2nd column is allele frequency thresholded (AFT)
    r_squared_delta = zeros(num_sims,2); % 1st column is allele frequency weighted (AFW), 2nd column is allele frequency thresholded (AFT)
    
    r_squared_c_UB = zeros(num_sims,2); % unbiased estimates weighted by S/(S-1), 1st column is allele frequency weighted (AFW), 2nd column is allele frequency thresholded (AFT)
    r_squared_delta_UB = zeros(num_sims,2); % unbiased estimates weighted by S/(S-1), 1st column is allele frequency weighted (AFW), 2nd column is allele frequency thresholded (AFT)
    
    Ne = zeros(num_sims,7);
    UB_factor = ones(num_sims,2);  % S/(S-1) weight for unbiased estimate based on Weir(1979). 1st column is r_squared_c, 2nd column is r_squared_delta.
    

    for iter = 1:num_sims
        
        % get simulated genotype data
        [genotypes] = generate_pop(mating_model,num_loci,num_allele_locus,pop_size,gens_2_sim,false,'null'); 

        [num_individuals, ~] = size(genotypes);   %get dimensions of genotype table

        [orig_num_alleles,orig_allele_freqs,orig_allele_states] = alleles(genotypes,num_individuals,num_loci); % get allele frequencies and states

        if print_allele_table == 1
            fprintf('Allele frequencies in genotype data:\n');
            allele_freq_table(num_loci, orig_num_alleles, orig_allele_freqs, orig_allele_states); % print table of allele states and frequencies for all loci
        end


        [r_squared_c(iter,1),r_squared_c(iter,2),r_squared_delta(iter,1),r_squared_delta(iter,2),r_squared_c_values,r_squared_delta_values,...
            S_c_values,S_delta_values,c_locus_pairs,delta_locus_pairs] = ...
            make_estimate(print_notes,genotypes,num_individuals,num_loci,print_table,threshold_allele_freq,print_locus_table,make_files);


        % compute an average for S, the number of individuals for each estimate of r. There is one version for r^2_c and another for
        % r^2_delta because of div zero instances that may impact one but not the other.
        % Waples and Do 2008 used "weighted harmonic mean with weights proportional to the nij"

        sum_S = sum(S_c_values);
        mean_S_c = sum_S/c_locus_pairs;

        inv_S_c = 1./S_c_values;
        sum_inv_S_c = sum(inv_S_c);
        har_mean_S_c = 1/((1./c_locus_pairs)*sum_inv_S_c);

        UB_factor(iter,1) = har_mean_S_c/(har_mean_S_c - 1);

        r_squared_c_UB(iter,:) = UB_factor(iter,1).*r_squared_c(iter,:);    % rescale to get unbiased r_squared_c

        sum_S_delta = sum(S_delta_values);
        mean_S_delta = sum_S_delta/delta_locus_pairs;

        inv_S_delta = 1./S_delta_values;
        sum_inv_S_delta = sum(inv_S_delta);
        har_mean_S_delta = 1/((1./delta_locus_pairs)*sum_inv_S_delta);

        UB_factor(iter,2) = har_mean_S_delta/(har_mean_S_delta - 1);

        r_squared_delta_UB(iter,:) = UB_factor(iter,2).*r_squared_delta(iter,:);    % rescale to get unbiased r_squared_delta

        % Estimate genetic association caused by finite sampling, see equations 13 and 14 in Sved et al. 2013. 
        correction_factor = (1/num_individuals)*(1-(1/((2*num_individuals - 1)^2))); %due to finite sampling of individuals
        correction_factor_avg_S = (1/mean_S_c)*(1-(1/((2*mean_S_c - 1)^2))); % arithmetic mean of sample sizes for locus pairs
        correction_factor_har_mean_S = (1/har_mean_S_c)*(1-(1/((2*har_mean_S_c - 1)^2))); % harmonic mean of sample sizes for locus pairs

        r_delta_correction_factor = (1/num_individuals)*(1-(1/((2*num_individuals - 1)^2))); %due to finite sampling of individuals
        r_delta_correction_factor_avg_S = (1/mean_S_delta)*(1-(1/((2*mean_S_delta - 1)^2))); % arithmetic mean of sample sizes for locus pairs
        r_delta_correction_factor_har_mean_S = (1/har_mean_S_delta)*(1-(1/((2*har_mean_S_delta - 1)^2))); % harmonic mean of sample sizes for locus pairs


        % compute Ne with simple expectation and finite sampling correction based on harmonic mean sample size
        Ne(iter,1) = 1/(3*(r_squared_c(iter,1) - correction_factor_har_mean_S)); % r_squared_c AFW
        Ne(iter,2) = 1/(3*(r_squared_c(iter,2) - correction_factor_har_mean_S)); % r_squared_c AFT
        Ne(iter,3) = 1/(3*(r_squared_c_UB(iter,1) - correction_factor_har_mean_S)); % r_squared_c AFW UB
        Ne(iter,4) = 1/(3*(r_squared_c_UB(iter,2) - correction_factor_har_mean_S)); % r_squared_c AFT UB

        Ne(iter,5) = 1/(3*(r_squared_delta(iter,1) - r_delta_correction_factor_har_mean_S)); % r_squared_delta AFW
        Ne(iter,6) = 1/(3*(r_squared_delta(iter,2) - r_delta_correction_factor_har_mean_S)); % r_squared_delta AFT
        Ne(iter,7) = 1/(3*(r_squared_delta_UB(iter,1) - r_delta_correction_factor_har_mean_S)); % r_squared_delta AFW UB
        Ne(iter,8) = 1/(3*(r_squared_delta_UB(iter,2) - r_delta_correction_factor_har_mean_S)); % r_squared_delta AFT UB


        % compute Waples (2006) statistical-fit corrected Ne estimate with each version of the sample size 
        % equations based on Table 1 in Waples and Do (2007). Includes
        % multiplication by UB_factor(iter,2) to give unbiased estimate of r_squared delta. 
        if har_mean_S_delta < 30
            E_r2_delta_sample = 0.0018 + 0.907/har_mean_S_delta + 4.44/(har_mean_S_delta^2);

            r_2_prime = r_squared_delta(iter,1) - E_r2_delta_sample;
            Ne(iter,9) = (0.308 + sqrt(0.308^2 + 2.08*r_2_prime^2) )/ (2*r_2_prime); % Waples r_squared_delta AFW

            r_2_prime = r_squared_delta(iter,2) - E_r2_delta_sample;
            Ne(iter,10) = (0.308 + sqrt(0.308^2 + 2.08*r_2_prime^2) )/ (2*r_2_prime); % Waples r_squared_delta AFT               

            r_2_prime = r_squared_delta_UB(iter,1) - E_r2_delta_sample;
            Ne(iter,11) = (0.308 + sqrt(0.308^2 + 2.08*r_2_prime^2) )/ (2*r_2_prime); % Waples r_squared_delta AFW UB

            r_2_prime = r_squared_delta_UB(iter,2) - E_r2_delta_sample;
            Ne(iter,12) = (0.308 + sqrt(0.308^2 + 2.08*r_2_prime^2) )/ (2*r_2_prime); % Waples r_squared_delta AFT UB

        else % har_mean_S_delta >= 30
            E_r2_delta_sample = 1/har_mean_S_delta + 3.19/(har_mean_S_delta^2);

            r_2_prime = r_squared_delta(iter,1) - E_r2_delta_sample;
            Ne(iter,9) = (1/3 + sqrt(1/9 + 2.76*r_2_prime^2) )/ (2*r_2_prime); % Waples r_squared_delta AFW

            r_2_prime = r_squared_delta(iter,2) - E_r2_delta_sample;
            Ne(iter,10) = (1/3 + sqrt(1/9 + 2.76*r_2_prime^2) )/ (2*r_2_prime); % Waples r_squared_delta AFT              

            r_2_prime = r_squared_delta_UB(iter,1) - E_r2_delta_sample;
            Ne(iter,11) = (1/3 + sqrt(1/9 + 2.76*r_2_prime^2) )/ (2*r_2_prime); % Waples r_squared_delta AFW UB

            r_2_prime = r_squared_delta_UB(iter,2) - E_r2_delta_sample;
            Ne(iter,12) = (1/3 + sqrt(1/9 + 2.76*r_2_prime^2) )/ (2*r_2_prime); % Waples r_squared_delta AFT UB

        end;

        % iteratively permute genotype data to estimate correction factor and associated Ne
        print_output = false;   % turn off text outputs
        make_table = false;   % turn off table outputs
        make_graphs = false;   % turn off graph outputs

        [~,~,~,~,median_r_sqr_c_AFW_perm,median_r_sqr_c_AFT_perm,~,~]= ...
        permute_data(genotypes,n_permutations,num_individuals,num_loci,threshold_allele_freq,r_squared_c(iter,1),r_squared_c(iter,2),...
        r_squared_delta(iter,1),r_squared_delta(iter,2),alpha,make_table,make_graphs,make_files);

        % compute Ne with simple expectation and finite sampling correction based on median permuted r^2
        Ne(iter,13) = 1/(3*(r_squared_c(iter,1) - median_r_sqr_c_AFW_perm)); % r_squared_c AFW permute
        Ne(iter,14) = 1/(3*(r_squared_c(iter,2) - median_r_sqr_c_AFT_perm)); % r_squared_c AFT permute

        Ne(iter,15) = 1/(3*(r_squared_c_UB(iter,1) - median_r_sqr_c_AFW_perm)); % r_squared_c AFW UB permute
        Ne(iter,16) = 1/(3*(r_squared_c_UB(iter,2) - median_r_sqr_c_AFT_perm)); % r_squared_c AFT UB permute


        % Compute Waples & Do chi-square CIs, "jackknife" CIs, and estimate of effective sample size n' and associated CIs.  
        % requires >= 2 loci
        [jack_pct_CI_values(iter,:),jack_norm_CI_values(iter,:),chi_CI_values(iter,:),effective_n_c(iter,:),effective_n_delta(iter,:)] = ...
        wd_confidence_intervals_v2(num_loci,r_squared_c_values,r_squared_delta_values,c_locus_pairs,delta_locus_pairs,...
        r_squared_c(iter,1),r_squared_c(iter,2),r_squared_delta(iter,1),r_squared_delta(iter,2),correction_factor_har_mean_S,...
        median_r_sqr_c_AFW_perm,median_r_sqr_c_AFT_perm,alpha,print_output,make_graphs);

        % jackknife over individuals
        [r_sqr_c_AFW_jack_individ_CI_pct(iter,:),r_sqr_c_AFT_jack_individ_CI_pct(iter,:),r_sqr_d_AFW_jack_individ_CI_pct(iter,:),...
        r_sqr_d_AFT_jack_individ_CI_pct(iter,:),r_sqr_c_AFW_jack_individ_CI_norm(iter,:),r_sqr_c_AFT_jack_individ_CI_norm(iter,:),...
        r_sqr_d_AFW_jack_individ_CI_norm(iter,:),r_sqr_d_AFT_jack_individ_CI_norm(iter,:),r_sqr_c_AFW_Ne_jack_individ_CI_pct(iter,:),...
        r_sqr_c_AFW_Ne_jack_individ_CI_norm(iter,:),r_sqr_c_AFT_Ne_jack_individ_CI_pct(iter,:),r_sqr_c_AFT_Ne_jack_individ_CI_norm(iter,:),...
        r_sqr_d_AFW_Ne_jack_individ_CI_pct(iter,:),r_sqr_d_AFT_Ne_jack_individ_CI_pct(iter,:),r_sqr_d_AFW_Ne_jack_individ_CI_norm(iter,:),...
        r_sqr_d_AFT_Ne_jack_individ_CI_norm(iter,:)] = ...
        jackknife_individuals(r_squared_c(iter,1),r_squared_c(iter,2),r_squared_delta(iter,1),r_squared_delta(iter,2),median_r_sqr_c_AFW_perm,...
        median_r_sqr_c_AFT_perm,correction_factor_har_mean_S,genotypes,num_individuals,num_loci,threshold_allele_freq,alpha,print_output,...
        make_graphs,make_files);

        % jackknife over loci
        [r_sqr_c_AFW_jack_loci_CI_pct(iter,:),r_sqr_c_AFT_jack_loci_CI_pct(iter,:),r_sqr_d_AFW_jack_loci_CI_pct(iter,:),...
        r_sqr_d_AFT_jack_loci_CI_pct(iter,:),r_sqr_c_AFW_jack_loci_CI_norm(iter,:),r_sqr_c_AFT_jack_loci_CI_norm(iter,:),...
        r_sqr_d_AFW_jack_loci_CI_norm(iter,:),r_sqr_d_AFT_jack_loci_CI_norm(iter,:),r_sqr_c_AFW_Ne_jack_loci_CI_pct(iter,:),...
        r_sqr_c_AFW_Ne_jack_loci_CI_norm(iter,:),r_sqr_c_AFT_Ne_jack_loci_CI_pct(iter,:),r_sqr_c_AFT_Ne_jack_loci_CI_norm(iter,:),...
        r_sqr_d_AFW_Ne_jack_loci_CI_pct(iter,:),r_sqr_d_AFT_Ne_jack_loci_CI_pct(iter,:),r_sqr_d_AFW_Ne_jack_loci_CI_norm(iter,:),...
        r_sqr_d_AFT_Ne_jack_loci_CI_norm(iter,:)] = ...
        jackknife_loci(r_squared_c(iter,1),r_squared_c(iter,2),r_squared_delta(iter,1),r_squared_delta(iter,2),median_r_sqr_c_AFW_perm,...
        median_r_sqr_c_AFT_perm,correction_factor_har_mean_S,genotypes,num_individuals,num_loci,threshold_allele_freq,alpha,print_output,...
        make_graphs,make_files);

        fprintf('simuation iteration %d of %d completed.\n',iter,num_sims);

    end; % for sim=1:num_sims

    % save simulation results as a matlab file for later analysis and graphing
    tNe = num2str(pop_size);
    full_file_name = strcat(output_filename,tNe,file_name_ext); 
    save(full_file_name); % save all variables from the current workspace to a .mat file. If filename exists, save overwrites the file.

    % some output to console
    fprintf('True Ne for this simulation = %d\n', pop_size);

    fprintf('Mean Ne for all 16 estimators: \n');
    mean(Ne)

    fprintf('Number of estimates less than zero for all 16 estimators: \n');
    counts = Ne<0;
    sum(counts) % show number of <=zero for all 16 estimators        
    
    
    fprintf('Simulation completed.\n');
    
    
