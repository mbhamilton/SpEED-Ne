function [r_sqr_c_AFW_pct_CI,r_sqr_c_AFT_pct_CI,r_sqr_c_AFW_norm_CI,r_sqr_c_AFT_norm_CI,...
    r_sqr_delta_AFW_pct_CI,r_sqr_delta_AFT_pct_CI,r_sqr_delta_AFW_norm_CI,r_sqr_delta_AFT_norm_CI,...
    r_sqr_c_AFW_Ne_pct_CI,r_sqr_c_AFW_Ne_norm_CI,r_sqr_c_AFT_Ne_pct_CI,r_sqr_c_AFT_Ne_norm_CI,...
    r_sqr_delta_AFW_Ne_pct_CI,r_sqr_delta_AFT_Ne_pct_CI,r_sqr_delta_AFW_Ne_norm_CI,r_sqr_delta_AFT_Ne_norm_CI] = ...
    jackknife_individuals(r_squared_c_W,r_squared_delta_W,r_squared_c_TH,r_squared_delta_TH,median_r_sqr_c_AFW_perm,median_r_sqr_c_AFT_perm,...
    correction_factor_har_mean_S,genotypes,num_individuals,num_loci,UB_weight_comp,UB_weight_delta,threshold_allele_freq,alpha,print_output,make_graphs,make_files)

%	Carry out a delete-one jackknife over individuals (rows) in genotype data.
%
%	version 1.3 17 Aug 2022
%
%   Inputs:
%   r_squared_c_W - estimate of r^2_{c} with weighting by allele frequency
%   r_squared_c_TH - estimate of r^2_{c} with allele frequency thresholding
%   r_squared_delta_W - estimate of r^2_{delta} with weighting by allele frequency
%   r_squared_delta_TH - estimate of r^2_{delta} with allele frequency thresholding
%   median_r_sqr_c_AFW_perm - median value from r_sqr_c_AFW from genotype permutation
%   median_r_sqr_c_AFT_perm - median value from r_sqr_c_AFT from genotype permutation
%   correction_factor_har_mean_S - r^2 correction factor from harmonic mean sample size used for r_squared_delta
%   genotypes - num_individuals by 2*num_loci matirx of genetic data with pairs of columns containing the two
%   	integer alleles for a diploid locus. All missing alleles should be
%   	coded as NaN.
%   num_individuals - integer number of rows in the genotype data matrix
%   UB_weight_comp - weighting for unbiased r_squared_comp har_mean_S_c/(har_mean_S_c - 1);
%   UB_weight_delta - weighting for unbiased r_squared_delta har_mean_S_delta/(har_mean_S_delta - 1)
%   num_loci - integer half the number of columns in the genotype data matrix
%   threshold_allele_freq - minimum allele frequency for estimate of r^2 to be used in average
%   alpha - probability mass in tail of confidence intervals such that lower
%       CI = (alpha/2) and upper %CI = 1 - (alpha/2).
%   print_output - boolean switch for fprintf output of values
%   make_graphs - boolean switch to output histograms of jackknife estimates in a four panel plot
%   make_files - boolean switch for output of results to files - requires
%       prep of results tables
%
%   Outputs:
%   columns 1 and 2 contain unweighted estimates, columns 3 and 4 contain unbiased weighted estimates
%   r_sqr_c_AFW_pct_CI - 1 by 4 vector of percentile alpha/2, 1-(alpha/2) CI for r_sqr_c AFW and r_sqr_c AFW UB
%   r_sqr_c_AFT_pct_CI - 1 by 4 vector of percentile alpha/2,1-(alpha/2) CI for r_sqr_c AFT and r_sqr_c AFT UB
%   r_sqr_delta_AFW_pct_CI - 1 by 4 vector of percentile alpha/2, 1-(alpha/2) CI for r_sqr_delta AFW and r_sqr_delta AFW UB
%   r_sqr_delta_AFT_pct_CI - 1 by 4 vector of percentile alpha/2,1-(alpha/2) CI for r_sqr_delta AFT and r_sqr_delta AFT UB
%   r_sqr_c_AFW_norm_CI - 1 by 4 vector of normal distribution alpha/2, 1-(alpha/2) CI for r_sqr_c AFW and r_sqr_c AFW UB
%   r_sqr_c_AFT_norm_CI - 1 by 4 vector of normal distribution alpha/2, 1-(alpha/2) CI for r_sqr_c AFT and r_sqr_c AFT UB
%   r_sqr_delta_AFW_norm_CI - 1 by 4 vector of normal distribution alpha/2, 1-(alpha/2) CI for r_sqr_delta AFW and r_sqr_delta AFW UB
%   r_sqr_delta_AFT_norm_CI - 1 by 4 vector of normal distribution alpha/2, 1-(alpha/2) CI for r_sqr_delta AFT and r_sqr_delta AFT UB
%
%   columns 1 and 2 contain unweighted estimates, columns 3 and 4 contain unbiased weighted estimates
%   r_sqr_c_AFW_Ne_pct_CI - 1 by 4 vector of percentile alpha/2,1-(alpha/2) CI for Ne estimated from r_sqr_c AFW and r_sqr_c AFW UB
%   r_sqr_c_AFW_Ne_norm_CI - 1 by 4 vector of percentile alpha/2,1-(alpha/2) CI for Ne estimated from r_sqr_c AFT and r_sqr_c AFT UB
%   r_sqr_c_AFT_Ne_pct_CI - 1 by 4 vector of percentile alpha/2,1-(alpha/2) CI for Ne estimated from r_sqr_delta AFW and r_sqr_delta AFW UB
%   r_sqr_c_AFT_Ne_norm_CI - 1 by 4 vector of percentile alpha/2,1-(alpha/2) CI for Ne estimated from r_sqr_delta AFT and r_sqr_delta AFT UB
%   r_sqr_delta_AFW_Ne_pct_CI - 1 by 4 vector of normal distribution alpha/2, 1-(alpha/2) CI for Ne estimated from r_sqr_c AFW and r_sqr_c AFW UB
%   r_sqr_delta_AFT_Ne_pct_CI - 1 by 4 vector of normal distribution alpha/2, 1-(alpha/2) CI for Ne estimated from r_sqr_c AFT and r_sqr_c AFT UB
%   r_sqr_delta_AFW_Ne_norm_CI - 1 by 4 vector of normal distribution alpha/2, 1-(alpha/2) CI for Ne estimated from r_sqr_delta AFW and r_sqr_delta AFW UB
%   r_sqr_delta_AFT_Ne_norm_CI - 1 by 4 vector of normal distribution alpha/2, 1-(alpha/2) CI for Ne estimated from r_sqr_delta AFT and r_sqr_delta AFT UB
%
%**************
%   Copyright 2022 by Matthew B Hamilton. 
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


	print_notes = false; % print error output strings to console
    print_table = false; 
    print_locus_table = false;

    %del_one_individual_genotypes = zeros(num_individuals-1,2*num_loci); % initialize

    fprintf('\n\n');

    %fprintf('jackknifing by deleting one individual of %d ...\n',num_individuals);

    h = waitbar(0,'jackknifing over individuals ...');

    
    for jack_iter = 1:num_individuals

        %fprintf('jackknife delete individual %d of %d...\n',jack_iter,num_individuals);

        waitbar(jack_iter / num_individuals);
        
        del_one_individual_genotypes = genotypes(~ismember(1:num_individuals, jack_iter),:); % get subset of individuals with jack_iter row deleted

        [jack_r_squared_c(jack_iter,1),jack_r_squared_c(jack_iter,2),jack_r_squared_delta(jack_iter,1),jack_r_squared_delta(jack_iter,2),~,~,~,~,~,~] = ...
            make_estimate(print_notes,del_one_individual_genotypes,num_individuals-1,num_loci,print_table,threshold_allele_freq,print_locus_table,make_files);

    end %for jack_iter = 1:num_loci

    close(h);

    fprintf('\n\n');
    
    % compute lower and upper bounds for all types of estimates (comp and delta, AFW and AFT, unweighted and UB weigthed)
    [r_sqr_c_AFW_pct_CI(1,1),r_sqr_c_AFW_pct_CI(1,2)] = percentile_CI(jack_r_squared_c(:,1),num_loci,alpha);     
    r_sqr_c_AFW_pct_CI(1,3) = r_sqr_c_AFW_pct_CI(1,1)*UB_weight_comp;
    r_sqr_c_AFW_pct_CI(1,4) = r_sqr_c_AFW_pct_CI(1,2)*UB_weight_comp;
    
    [r_sqr_c_AFT_pct_CI(1,1),r_sqr_c_AFT_pct_CI(1,2)] = percentile_CI(jack_r_squared_c(:,2),num_loci,alpha);
    r_sqr_c_AFT_pct_CI(1,3) = r_sqr_c_AFT_pct_CI(1,1)*UB_weight_comp;
    r_sqr_c_AFT_pct_CI(1,4) = r_sqr_c_AFT_pct_CI(1,2)*UB_weight_comp;

    [r_sqr_delta_AFW_pct_CI(1,1),r_sqr_delta_AFW_pct_CI(1,2)] = percentile_CI(jack_r_squared_delta(:,1),num_loci,alpha);
    r_sqr_delta_AFW_pct_CI(1,3) = r_sqr_delta_AFW_pct_CI(1,1)*UB_weight_delta;
    r_sqr_delta_AFW_pct_CI(1,4) = r_sqr_delta_AFW_pct_CI(1,2)*UB_weight_delta;
    
    [r_sqr_delta_AFT_pct_CI(1,1),r_sqr_delta_AFT_pct_CI(1,2)] = percentile_CI(jack_r_squared_delta(:,2),num_loci,alpha);
    r_sqr_delta_AFT_pct_CI(1,3) = r_sqr_delta_AFT_pct_CI(1,1)*UB_weight_delta;
    r_sqr_delta_AFT_pct_CI(1,4) = r_sqr_delta_AFT_pct_CI(1,2)*UB_weight_delta;
    
    [r_sqr_c_AFW_norm_CI(1,1),r_sqr_c_AFW_norm_CI(1,2)] = normal_dist_CI(r_squared_c_W,jack_r_squared_c(:,1),num_loci,alpha);
    r_sqr_c_AFW_norm_CI(1,3) = r_sqr_c_AFW_norm_CI(1,1)*UB_weight_comp;
    r_sqr_c_AFW_norm_CI(1,4) = r_sqr_c_AFW_norm_CI(1,2)*UB_weight_comp;
    
    [r_sqr_c_AFT_norm_CI(1,1),r_sqr_c_AFT_norm_CI(1,2)] = normal_dist_CI(r_squared_c_TH,jack_r_squared_c(:,2),num_loci,alpha);
    r_sqr_c_AFT_norm_CI(1,3) = r_sqr_c_AFT_norm_CI(1,1)*UB_weight_comp;
    r_sqr_c_AFT_norm_CI(1,4) = r_sqr_c_AFT_norm_CI(1,2)*UB_weight_comp;

    [r_sqr_delta_AFW_norm_CI(1,1),r_sqr_delta_AFW_norm_CI(1,2)] = normal_dist_CI(r_squared_delta_W,jack_r_squared_delta(:,1),num_loci,alpha);
    r_sqr_delta_AFW_norm_CI(1,3) = r_sqr_delta_AFW_norm_CI(1,1)*UB_weight_delta;
    r_sqr_delta_AFW_norm_CI(1,4) = r_sqr_delta_AFW_norm_CI(1,2)*UB_weight_delta;
    
    [r_sqr_delta_AFT_norm_CI(1,1),r_sqr_delta_AFT_norm_CI(1,2)] = normal_dist_CI(r_squared_delta_TH,jack_r_squared_delta(:,2),num_loci,alpha);
    r_sqr_delta_AFT_norm_CI(1,3) = r_sqr_delta_AFT_norm_CI(1,1)*UB_weight_delta;
    r_sqr_delta_AFT_norm_CI(1,4) = r_sqr_delta_AFT_norm_CI(1,2)*UB_weight_delta;
    
    % estimate bounds for Ne based on various estimates of r^2 and CIs

    % r^2_{comp} AFW:
    r_sqr_c_AFW_Ne_pct_CI(1,1) = 1/(3*(r_sqr_c_AFW_pct_CI(1,2) - median_r_sqr_c_AFW_perm));
    r_sqr_c_AFW_Ne_pct_CI(1,2) = 1/(3*(r_sqr_c_AFW_pct_CI(1,1) - median_r_sqr_c_AFW_perm));
    
    r_sqr_c_AFW_Ne_pct_CI(1,3) = 1/(3*(r_sqr_c_AFW_pct_CI(1,2)*UB_weight_comp - median_r_sqr_c_AFW_perm*UB_weight_comp));
    r_sqr_c_AFW_Ne_pct_CI(1,4) = 1/(3*(r_sqr_c_AFW_pct_CI(1,1)*UB_weight_comp - median_r_sqr_c_AFW_perm*UB_weight_comp));
    
    r_sqr_c_AFW_Ne_norm_CI(1,1) = 1/(3*(r_sqr_c_AFW_norm_CI(1,2) - median_r_sqr_c_AFW_perm));
    r_sqr_c_AFW_Ne_norm_CI(1,2) = 1/(3*(r_sqr_c_AFW_norm_CI(1,1) - median_r_sqr_c_AFW_perm));
    
    r_sqr_c_AFW_Ne_norm_CI(1,3) = 1/(3*(r_sqr_c_AFW_norm_CI(1,2)*UB_weight_comp - median_r_sqr_c_AFW_perm*UB_weight_comp));
    r_sqr_c_AFW_Ne_norm_CI(1,4) = 1/(3*(r_sqr_c_AFW_norm_CI(1,1)*UB_weight_comp - median_r_sqr_c_AFW_perm*UB_weight_comp));
    
    % r^2_{comp} AFT:
    r_sqr_c_AFT_Ne_pct_CI(1,1) = 1/(3*(r_sqr_c_AFT_pct_CI(1,2) - median_r_sqr_c_AFT_perm));
    r_sqr_c_AFT_Ne_pct_CI(1,2) = 1/(3*(r_sqr_c_AFT_pct_CI(1,1) - median_r_sqr_c_AFT_perm));
    
    r_sqr_c_AFT_Ne_pct_CI(1,3) = 1/(3*(r_sqr_c_AFT_pct_CI(1,2)*UB_weight_comp - median_r_sqr_c_AFT_perm*UB_weight_comp));
    r_sqr_c_AFT_Ne_pct_CI(1,4) = 1/(3*(r_sqr_c_AFT_pct_CI(1,1)*UB_weight_comp - median_r_sqr_c_AFT_perm*UB_weight_comp));
    
    r_sqr_c_AFT_Ne_norm_CI(1,1) = 1/(3*(r_sqr_c_AFT_norm_CI(1,2) - median_r_sqr_c_AFT_perm));
    r_sqr_c_AFT_Ne_norm_CI(1,2) = 1/(3*(r_sqr_c_AFT_norm_CI(1,1) - median_r_sqr_c_AFT_perm));
    
    r_sqr_c_AFT_Ne_norm_CI(1,3) = 1/(3*(r_sqr_c_AFT_norm_CI(1,2)*UB_weight_comp - median_r_sqr_c_AFT_perm*UB_weight_comp));
    r_sqr_c_AFT_Ne_norm_CI(1,4) = 1/(3*(r_sqr_c_AFT_norm_CI(1,1)*UB_weight_comp - median_r_sqr_c_AFT_perm*UB_weight_comp));
    
    % r^2_{delta} AFW:
    r_sqr_delta_AFW_Ne_pct_CI(1,1) = 1/(3*(r_sqr_delta_AFW_pct_CI(1,2) - correction_factor_har_mean_S));
    r_sqr_delta_AFW_Ne_pct_CI(1,2) = 1/(3*(r_sqr_delta_AFW_pct_CI(1,1) - correction_factor_har_mean_S));

    r_sqr_delta_AFW_Ne_pct_CI(1,3) = 1/(3*(r_sqr_delta_AFW_pct_CI(1,2)*UB_weight_delta - correction_factor_har_mean_S));
    r_sqr_delta_AFW_Ne_pct_CI(1,4) = 1/(3*(r_sqr_delta_AFW_pct_CI(1,1)*UB_weight_delta - correction_factor_har_mean_S));
    
    r_sqr_delta_AFW_Ne_norm_CI(1,1) = 1/(3*(r_sqr_delta_AFW_norm_CI(1,2) - correction_factor_har_mean_S));
    r_sqr_delta_AFW_Ne_norm_CI(1,2) = 1/(3*(r_sqr_delta_AFW_norm_CI(1,1) - correction_factor_har_mean_S));
        
    r_sqr_delta_AFW_Ne_norm_CI(1,3) = 1/(3*(r_sqr_delta_AFW_norm_CI(1,2)*UB_weight_delta - correction_factor_har_mean_S));
    r_sqr_delta_AFW_Ne_norm_CI(1,4) = 1/(3*(r_sqr_delta_AFW_norm_CI(1,1)*UB_weight_delta - correction_factor_har_mean_S));
    
    % r^2_{delta} AFT:
    r_sqr_delta_AFT_Ne_pct_CI(1,1) = 1/(3*(r_sqr_delta_AFT_pct_CI(1,2) - correction_factor_har_mean_S));
    r_sqr_delta_AFT_Ne_pct_CI(1,2) = 1/(3*(r_sqr_delta_AFT_pct_CI(1,1) - correction_factor_har_mean_S));
 
    r_sqr_delta_AFT_Ne_pct_CI(1,3) = 1/(3*(r_sqr_delta_AFT_pct_CI(1,2)*UB_weight_delta - correction_factor_har_mean_S));
    r_sqr_delta_AFT_Ne_pct_CI(1,4) = 1/(3*(r_sqr_delta_AFT_pct_CI(1,1)*UB_weight_delta - correction_factor_har_mean_S));
    
    r_sqr_delta_AFT_Ne_norm_CI(1,1) = 1/(3*(r_sqr_delta_AFT_norm_CI(1,2) - correction_factor_har_mean_S));
    r_sqr_delta_AFT_Ne_norm_CI(1,2) = 1/(3*(r_sqr_delta_AFT_norm_CI(1,1) - correction_factor_har_mean_S));

    r_sqr_delta_AFT_Ne_norm_CI(1,3) = 1/(3*(r_sqr_delta_AFT_norm_CI(1,2)*UB_weight_delta - correction_factor_har_mean_S));
    r_sqr_delta_AFT_Ne_norm_CI(1,4) = 1/(3*(r_sqr_delta_AFT_norm_CI(1,1)*UB_weight_delta - correction_factor_har_mean_S));

    
    if print_output

        fprintf('\n');
        fprintf('========================================================================================\n');
        fprintf('Results of delete-one jackknifing over indivduals.\n');
        fprintf('Ne = 1/(3*(r^2 - r^2 correction factor)).\n');
        fprintf('r^2 correction factor is median r^2 permute for r^2_{c} and harmonic mean S for r^2_{delta}.\n');
        fprintf('Estimates are either allele frequency weighted (AFW) or allele frequency thresholded (AFT).\n');
        fprintf('Estimates are with and without Weir''s (1979) S/(S-1) weighting of r^2 for bias correction.\n');
        fprintf('For unbiased estimates, r^2_{c} correction factor is (median r^2 permute)*harmonic mean S\n');
        fprintf('and is unchanged for r^2_{delta}.\n');
        fprintf('----------------------------------------------------------------------------------------\n');
        fprintf('Estimates without Weir''s (1979) S/(S-1) weighting of r^2 for bias correction.\n\n');
        
        fprintf('r^2_{comp} AFW:\n');
        fprintf('\tpoint estimate & percentile confidence interval: %f (%f to %f)\n',r_squared_c_W,r_sqr_c_AFW_pct_CI(1,1),r_sqr_c_AFW_pct_CI(1,2));
        fprintf('\tNe estimate: %f ',1/(3*(r_squared_c_W - median_r_sqr_c_AFW_perm)) );
        
        if r_sqr_c_AFW_Ne_pct_CI(1,1) >= 0
        	fprintf(' (%8.6f to ',r_sqr_c_AFW_Ne_pct_CI(1,1));
        else
            fprintf('(%i to ',0);
        end

        if r_sqr_c_AFW_Ne_pct_CI(1,2) >= 0 
            fprintf('%f)\n',r_sqr_c_AFW_Ne_pct_CI(1,2));
        else
            fprintf('%s)\n',' infinity');
        end
        fprintf('\n');

        fprintf('\tpoint estimate & normal distribution confidence interval: %f (%f to %f)\n',r_squared_c_W,r_sqr_c_AFW_norm_CI(1,1),r_sqr_c_AFW_norm_CI(1,2));
        fprintf('\tNe estimate: %f ',1/(3*(r_squared_c_W - median_r_sqr_c_AFW_perm)) );
        
        if r_sqr_c_AFW_Ne_pct_CI(1,1) >= 0
        	fprintf(' (%f to ',r_sqr_c_AFW_Ne_norm_CI(1,1));
        else
            fprintf('(%i to ',0);
        end

        if r_sqr_c_AFW_Ne_norm_CI(1,2) >= 0 
            fprintf('%f)\n',r_sqr_c_AFW_Ne_norm_CI(1,2));
        else
            fprintf('%s)\n',' infinity');
        end
        fprintf('\n');

        
        fprintf('r^2_{comp} AFT:\n');
        fprintf('\tpoint estimate & percentile confidence interval: %f (%f to %f)\n',r_squared_c_TH,r_sqr_c_AFT_pct_CI(1,1),r_sqr_c_AFT_pct_CI(1,2));        
        fprintf('\tNe estimate: %f ',1/(3*(r_squared_c_TH - median_r_sqr_c_AFT_perm)) );
        
        if r_sqr_c_AFT_Ne_pct_CI(1,1) >= 0
        	fprintf(' (%f to ',r_sqr_c_AFT_Ne_pct_CI(1,1));
        else
            fprintf('(%i to ',0);
        end

        if r_sqr_c_AFT_Ne_pct_CI(1,2) >= 0 
            fprintf('%f)\n',r_sqr_c_AFT_Ne_pct_CI(1,2));
        else
            fprintf('%s)\n',' infinity');
        end
        fprintf('\n');

        fprintf('\tpoint estimate & normal distribution confidence interval: %f (%f to %f)\n',r_squared_c_TH,r_sqr_c_AFT_norm_CI(1,1),r_sqr_c_AFT_norm_CI(1,2));        
        fprintf('\tNe estimate: %f ',1/(3*(r_squared_c_TH - median_r_sqr_c_AFT_perm)) );
        
        if r_sqr_c_AFT_Ne_norm_CI(1,1) >= 0
        	fprintf(' (%f to ',r_sqr_c_AFT_Ne_norm_CI(1,1));
        else
            fprintf('(%i to ',0);
        end

        if r_sqr_c_AFT_Ne_norm_CI(1,2) >= 0 
            fprintf('%8.6f)\n',r_sqr_c_AFT_Ne_norm_CI(1,2));
        else
            fprintf('%s)\n',' infinity');
        end
        fprintf('\n');

        
        fprintf('r^2_{delta} AFW:\n');
        fprintf('\tpoint estimate & percentile confidence interval: %f (%f to %f)\n',r_squared_delta_W,r_sqr_delta_AFW_pct_CI(1,1),r_sqr_delta_AFW_pct_CI(1,2));
        fprintf('\tNe estimate: %f ',1/(3*(r_squared_delta_W - correction_factor_har_mean_S)) );
        
        if r_sqr_delta_AFW_Ne_pct_CI(1,1) >= 0
        	fprintf(' (%f to ',r_sqr_delta_AFW_Ne_pct_CI(1,1));
        else
            fprintf('(%i to ',0);
        end

        if r_sqr_delta_AFW_Ne_pct_CI(1,2) >= 0 
            fprintf('%f)\n',r_sqr_delta_AFW_Ne_pct_CI(1,2));
        else
            fprintf('%s)\n',' infinity');
        end
        fprintf('\n');

        fprintf('\tpoint estimate & normal distribution confidence interval: %f (%f to %f)\n',r_squared_delta_W,r_sqr_delta_AFW_norm_CI(1,1),r_sqr_delta_AFW_norm_CI(1,2));
        fprintf('\tNe estimate: %f ',1/(3*(r_squared_delta_W - correction_factor_har_mean_S)) );
        
        if r_sqr_delta_AFW_Ne_norm_CI(1,1) >= 0
        	fprintf(' (%f to ',r_sqr_delta_AFW_Ne_norm_CI(1,1));
        else
            fprintf('(%i to ',0);
        end

        if r_sqr_delta_AFW_Ne_norm_CI(1,2) >= 0 
            fprintf('%f)\n',r_sqr_delta_AFW_Ne_norm_CI(1,2));
        else
            fprintf('%s)\n',' infinity');
        end
        fprintf('\n');

        
        fprintf('r^2_{delta} AFT:\n');
        fprintf('\tpoint estimate & percentile confidence interval: %f (%f to %f)\n',r_squared_delta_TH,r_sqr_delta_AFT_pct_CI(1,1),r_sqr_delta_AFT_pct_CI(1,2));
        fprintf('\tNe estimate: %f ',1/(3*(r_squared_delta_TH - correction_factor_har_mean_S)) );
        
        if r_sqr_delta_AFT_Ne_pct_CI(1,1) >= 0
        	fprintf(' (%f to ',r_sqr_delta_AFT_Ne_pct_CI(1,1));
        else
            fprintf('(%i to ',0);
        end

        if r_sqr_delta_AFT_Ne_pct_CI(1,2) >= 0 
            fprintf('%f)\n',r_sqr_delta_AFT_Ne_pct_CI(1,2));
        else
            fprintf('%s)\n',' infinity');
        end
        fprintf('\n');

        fprintf('\tpoint estimate & normal distribution confidence interval: %f (%f to %f)\n',r_squared_delta_TH,r_sqr_delta_AFT_norm_CI(1,1),r_sqr_delta_AFT_norm_CI(1,2));
        fprintf('\tNe estimate: %f ',1/(3*(r_squared_delta_TH - correction_factor_har_mean_S)) );
        
        if r_sqr_delta_AFT_Ne_norm_CI(1,1) >= 0
        	fprintf(' (%f to ',r_sqr_delta_AFT_Ne_norm_CI(1,1));
        else
            fprintf('(%i to ',0);
        end

        if r_sqr_delta_AFT_Ne_norm_CI(1,2) >= 0 
            fprintf('%8.6f)\n',r_sqr_delta_AFT_Ne_norm_CI(1,2));
        else
            fprintf('%s)\n',' infinity');
        end
        fprintf('\n');

        fprintf('-------------------\n');
        
        fprintf('Estimates with Weir''s (1979) S/(S-1) weighting of r^2 for bias correction.\n\n');
        
        r_squared_c_W_UB = r_squared_c_W * UB_weight_comp;
        r_squared_c_TH_UB = r_squared_c_TH * UB_weight_comp;
        r_squared_delta_W_UB = r_squared_delta_W * UB_weight_delta;
        r_squared_delta_TH_UB = r_squared_delta_TH * UB_weight_delta;
        
        fprintf('r^2_{comp} AFW:\n');
        fprintf('\tpoint estimate & percentile confidence interval: %f (%f to %f)\n',r_squared_c_W_UB,r_sqr_c_AFW_pct_CI(1,3),r_sqr_c_AFW_pct_CI(1,4));
        fprintf('\tNe estimate: %f ',1/(3*(r_squared_c_W_UB - median_r_sqr_c_AFW_perm*UB_weight_comp)) );
        
        if r_sqr_c_AFW_Ne_pct_CI(1,3) >= 0
        	fprintf(' (%f to ',r_sqr_c_AFW_Ne_pct_CI(1,3));
        else
            fprintf('(%i to ',0);
        end
        
        if r_sqr_c_AFW_Ne_pct_CI(1,4) >= 0 
            fprintf('%f)\n',r_sqr_c_AFW_Ne_pct_CI(1,4));
        else
            fprintf('%s)\n',' infinity');
        end
        fprintf('\n');

        fprintf('\tpoint estimate & normal distribution confidence interval: %f (%f to %f)\n',r_squared_c_W_UB,r_sqr_c_AFW_norm_CI(1,3),r_sqr_c_AFW_norm_CI(1,4));
        fprintf('\tNe estimate: %f ',1/(3*(r_squared_c_W_UB - median_r_sqr_c_AFW_perm*UB_weight_comp)) );
        
        if r_sqr_c_AFW_Ne_norm_CI(1,3) >= 0
        	fprintf(' (%f to ',r_sqr_c_AFW_Ne_norm_CI(1,3));
        else
            fprintf('(%i to ',0);
        end

        if r_sqr_c_AFW_Ne_norm_CI(1,4) >= 0 
            fprintf('%f)\n',r_sqr_c_AFW_Ne_norm_CI(1,4));
        else
            fprintf('%s)\n',' infinity');
        end
        fprintf('\n');

        
        fprintf('r^2_{comp} AFT:\n');
        fprintf('\tpoint estimate & percentile distribution confidence interval: %f (%f to %f)\n',r_squared_c_TH_UB,r_sqr_c_AFT_pct_CI(1,3),r_sqr_c_AFT_pct_CI(1,4));
        fprintf('\tNe estimate: %f ',1/(3*(r_squared_c_TH_UB - median_r_sqr_c_AFT_perm*UB_weight_comp)) );
        
        if r_sqr_c_AFT_Ne_pct_CI(1,3) >= 0
        	fprintf(' (%f to ',r_sqr_c_AFT_Ne_pct_CI(1,3));
        else
            fprintf('(%i to ',0);
        end

        if r_sqr_c_AFT_Ne_pct_CI(1,4) >= 0 
            fprintf('%f)\n',r_sqr_c_AFT_Ne_pct_CI(1,4));
        else
            fprintf('%s)\n',' infinity');
        end
        fprintf('\n');

        fprintf('\tpoint estimate & normal distribution confidence interval: %f (%f to %f)\n',r_squared_c_TH_UB,r_sqr_c_AFT_norm_CI(1,3),r_sqr_c_AFT_norm_CI(1,4));
        fprintf('\tNe estimate: %f ',1/(3*(r_squared_c_TH_UB - median_r_sqr_c_AFT_perm*UB_weight_comp)) );
        
        if r_sqr_c_AFT_Ne_norm_CI(1,3) >= 0
        	fprintf(' (%f to ',r_sqr_c_AFT_Ne_norm_CI(1,3));
        else
            fprintf('(%i to ',0);
        end

        if r_sqr_c_AFT_Ne_norm_CI(1,4) >= 0 
            fprintf('%f)\n',r_sqr_c_AFT_Ne_norm_CI(1,4));
        else
            fprintf('%s)\n',' infinity');
        end
        fprintf('\n');

        
        fprintf('r^2_{delta} AFW:\n');
        fprintf('\tpoint estimate & percentile confidence interval: %f (%f to %f)\n',r_squared_delta_W_UB,r_sqr_delta_AFW_pct_CI(1,3),r_sqr_delta_AFW_pct_CI(1,4));
        fprintf('\tNe estimate: %f ',1/(3*(r_squared_delta_W_UB - correction_factor_har_mean_S)) );
        
        if r_sqr_delta_AFW_Ne_pct_CI(1,3) >= 0
        	fprintf(' (%f to ',r_sqr_delta_AFW_Ne_pct_CI(1,3));
        else
            fprintf('(%i to ',0);
        end

        if r_sqr_delta_AFW_Ne_pct_CI(1,4) >= 0 
            fprintf('%f)\n',r_sqr_delta_AFW_Ne_pct_CI(1,4));
        else
            fprintf('%s)\n',' infinity');
        end
        fprintf('\n');

        
        fprintf('\tpoint estimate & normal distribution confidence interval: %f (%f to %f)\n',r_squared_delta_W_UB,r_sqr_delta_AFW_norm_CI(1,3),r_sqr_delta_AFW_norm_CI(1,4));
        fprintf('\tNe estimate: %f ',1/(3*(r_squared_delta_W_UB - correction_factor_har_mean_S)) );
        
        if r_sqr_delta_AFW_Ne_norm_CI(1,3) >= 0
        	fprintf(' (%f to ',r_sqr_delta_AFW_Ne_norm_CI(1,3));
        else
            fprintf('(%i to ',0);
        end

        if r_sqr_delta_AFW_Ne_norm_CI(1,4) >= 0 
            fprintf('%f)\n',r_sqr_delta_AFW_Ne_norm_CI(1,4));
        else
            fprintf('%s)\n',' infinity');
        end
        fprintf('\n');

        
        fprintf('r^2_{delta} AFT:\n');
        fprintf('\tpoint estimate & percentile confidence interval: %f (%f to %f)\n',r_squared_delta_TH_UB,r_sqr_delta_AFT_pct_CI(1,3),r_sqr_delta_AFT_pct_CI(1,4));
        fprintf('\tNe estimate: %f ',1/(3*(r_squared_delta_TH_UB - correction_factor_har_mean_S)) );
        
        if r_sqr_delta_AFT_Ne_pct_CI(1,3) >= 0
        	fprintf(' (%f to ',r_sqr_delta_AFT_Ne_pct_CI(1,3));
        else
            fprintf('(%i to ',0);
        end

        if r_sqr_delta_AFT_Ne_pct_CI(1,4) >= 0 
            fprintf('%f)\n',r_sqr_delta_AFT_Ne_pct_CI(1,4));
        else
            fprintf('%s)\n',' infinity');
        end
        fprintf('\n');

        fprintf('\tpoint estimate & normal distribution confidence interval: %f (%f to %f)\n',r_squared_delta_TH_UB,r_sqr_delta_AFT_norm_CI(1,3),r_sqr_delta_AFT_norm_CI(1,4));
        fprintf('\tNe estimate: %f ',1/(3*(r_squared_delta_TH_UB - correction_factor_har_mean_S)) );
        
        if r_sqr_delta_AFT_Ne_norm_CI(1,3) >= 0
        	fprintf(' (%f to ',r_sqr_delta_AFT_Ne_norm_CI(1,3));
        else
            fprintf('(%i to ',0);
        end

        if r_sqr_delta_AFT_Ne_norm_CI(1,4) >= 0 
            fprintf('%f)\n',r_sqr_delta_AFT_Ne_norm_CI(1,4));
        else
            fprintf('%s)\n',' infinity');
        end
        
        fprintf('-------------------\n');
        fprintf('Confidence intervals are %g%% to %g%%\n',100*alpha/2,100*(1-alpha/2));
        fprintf('========================================================================================\n');

    end % if print_output

    
    if make_graphs
        % show distribution of jackknife estimates to evalute distribution
        figure;

        hold on;
        subplot(2,2,1);
        %hist(jack_r_squared_c(:,1)); 
        
        histogram(jack_r_squared_c(:,1));
        xlabel('r_c^2 AFW')
        ylabel('Count')

        subplot(2,2,2);
        %hist(jack_r_squared_c(:,2)); 
        
        histogram(jack_r_squared_c(:,2));
        xlabel('r_c^2 AFT')
        ylabel('Count')

        subplot(2,2,3)
        %hist(jack_r_squared_delta(:,1)); 
        
        histogram(jack_r_squared_delta(:,1));
        xlabel('r_{\Delta}^2 AFW')
        ylabel('Count')

        subplot(2,2,4)
        %hist(jack_r_squared_delta(:,2)); 
        
        histogram(jack_r_squared_delta(:,2));
        xlabel('r_{\Delta}^2 AFT')
        ylabel('Count')

        sgtitle('Distributions of r^2 sub-estimates from delete-one jackknife over individuals.');
        
        hold off;

    end % if make_graphs
