function [mean_r_sqr_c_W_perm,mean_r_sqr_c_TH_perm,mean_r_sqr_delta_W_perm,mean_r_sqr_delta_TH_perm,...
    median_r_sqr_c_W_perm,median_r_sqr_c_TH_perm,median_r_sqr_delta_W_perm,median_r_sqr_delta_TH_perm] = ...
    permute_data(genotypes,n_permutations,num_individuals,num_loci,UB_weight_comp,UB_weight_delta,threshold_allele_freq,r_squared_c_W,r_squared_c_TH,r_squared_delta_W,r_squared_delta_TH,alpha_value,make_table,make_graphs,make_files)
%
%   version 1.4 11 March 2024
%
% Function to permute genotype data based on Sved et al. (2013).
%
% Inputs:
%   genotypes - num_individuals by 2*num_loci matrix of genotype data
%   n_permutations - number of permutations
%   num_individuals - number of rows in genotype matrix
%   num_loci - integer half the number of columns in the genotype data matrix
%   UB_weight_comp - weighting for unbiased r_squared_comp har_mean_S_c/(har_mean_S_c - 1);
%   UB_weight_delta - weighting for unbiased r_squared_delta har_mean_S_delta/(har_mean_S_delta - 1)
%   threshold_allele_freq - minimum frequency for AFT estimates of r^2
%   r_squared_c_W - multilocus estimate of r^2_{comp} with AFW
%   r_squared_c_TH - multilocus estimate of r^2_{comp} with AFT
%   r_squared_delta_W - multilocus estimate of r^2_{delta} with AFW
%   r_squared_delta_TH - multilocus estimate of r^2_{delta} with AFT   
%   alpha_value - probability mass in tail of confidence intervals such that lower
%       CI = (alpha/2) and upper %CI = 1 - (alpha/2).
%   make_table - boolean for output of results table
%   make_graphs - boolean for output of graphs
%
%   Outputs:
%   mean_r_sqr_c_W_perm - mean of r_sqr_c_AFW from the permuted data sets
%   mean_r_sqr_delta_W_perm - mean of r_sqr_delta_AFW from the permuted data sets
%   mean_r_sqr_c_TH_perm - mean of r_sqr_c_AFT from the permuted data sets
%   mean_r_sqr_delta_TH_perm - mean of r_sqr_c_AFT from the permuted data sets
%   median_r_sqr_c_W_perm - median of r_sqr_c_AFW from the permuted data sets
%   median_r_sqr_delta_W_perm - median of r_sqr_delta_AFW from the permuted data sets
%   median_r_sqr_c_TH_perm - median of r_sqr_c_AFT from the permuted data sets
%   median_r_sqr_delta_TH_perm - median of r_sqr_c_AFT from the permuted data sets
%
%**************
%   Copyright 2024 Matthew B Hamilton.
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

    % carry out permutation of genotypes and generate vectors of permute
    % estimates of r^2 measures
    [r_squared_c_W_permute,r_squared_c_TH_permute,r_squared_delta_W_permute,r_squared_delta_TH_permute] = ...
        permute(genotypes,n_permutations,num_individuals,num_loci,threshold_allele_freq,make_files);


    % produce graphs of the distributions of permuted values
    if make_graphs == 1
        figure;        
        hold on;
        subplot(2,2,1);

        histogram(r_squared_c_W_permute); 
        xlabel('r_c^2 AFW permute')
        ylabel('Count')

        subplot(2,2,2);
        histogram(r_squared_c_TH_permute); 
        xlabel('r_c^2 AFT permute')
        ylabel('Count')

        subplot(2,2,3)
        histogram(r_squared_delta_W_permute);
        xlabel('r_{\Delta}^2 AFW permute')
        ylabel('Count')

        subplot(2,2,4)
        histogram(r_squared_delta_TH_permute); 
        xlabel('r_{\Delta}^2 AFT permute')
        ylabel('Count')

        sgtitle('Distributions of r^2 values with permuted genotype data.');

        hold off;

    end % if make_graphs

    % get means
    mean_r_sqr_c_W_perm = mean(r_squared_c_W_permute); 
    mean_r_sqr_c_TH_perm = mean(r_squared_c_TH_permute);
    mean_r_sqr_delta_W_perm = mean(r_squared_delta_W_permute); 
    mean_r_sqr_delta_TH_perm = mean(r_squared_delta_TH_permute);

    % get medians
    median_r_sqr_c_W_perm = median(r_squared_c_W_permute); 
    median_r_sqr_c_TH_perm = median(r_squared_c_TH_permute);
    median_r_sqr_delta_W_perm = median(r_squared_delta_W_permute); 
    median_r_sqr_delta_TH_perm = median(r_squared_delta_TH_permute);
    
    if make_table

        fprintf('\n');
        fprintf('==============================================================================\n');
        fprintf('Mean and median r^2 from genotype data permuted among loci. This value is used as\n');
        fprintf('the correction factor in Ne = 1/(3*(r^2 - r^2 correction factor)).\n');
        fprintf('Estimates are for both allele frequency weighted (AFW) and allele\n');
        fprintf('frequency thresholded (AFT) estimators.\n');
        fprintf('------------------------------------------------------------------------------\n');
        fprintf('Permuted genotypes %d times\n\n', n_permutations);
        fprintf('---------------------\n');
        fprintf('permuted r^2_{comp} allele frequency weighted (AFW):\n');
        fprintf('\tmean = %8.6g\n',mean_r_sqr_c_W_perm);
        fprintf('\tmedian = %8.6g\n\n',median_r_sqr_c_W_perm);

        Ne = 1/(3*(r_squared_c_W - mean_r_sqr_c_W_perm));
        if Ne >= 0
            fprintf('\tNe (mean permuted r^2 as correction factor) = %8.6g\n', Ne);
        else
            fprintf('\tNe (mean permuted r^2 as correction factor) = infinity\n');
        end

        Ne = 1/(3*(r_squared_c_W - median_r_sqr_c_W_perm));
        if Ne >= 0
            fprintf('\tNe (median permuted r^2 as correction factor) = %8.6g\n', Ne);
        else
            fprintf('\tNe (median permuted r^2 as correction factor) = infinity\n');
        end

        
        [lower_CI,upper_CI] = percentile_CI(r_squared_c_W_permute,n_permutations,alpha_value);
        
        upper_bound_Ne = 1/(3*(r_squared_c_W - upper_CI));
        lower_bound_Ne = 1/(3*(r_squared_c_W - lower_CI));
        if upper_bound_Ne >= 0 
            fprintf('\tNe %g%% confidence intervals:\t%8.6g - %8.6g\n',(1-alpha_value)*100,lower_bound_Ne,upper_bound_Ne);
        else
            fprintf('\tNe %g%% confidence intervals:\t%8.6g - %s\n',(1-alpha_value)*100,lower_bound_Ne,'infinity');
        end

        fprintf('---------------------\n');

        fprintf('permuted r^2_{comp} allele frequency thresholded (AFT):\n');
        fprintf('\tmean = %8.6g\n',mean_r_sqr_c_TH_perm);
        fprintf('\tmedian = %8.6g\n\n',median_r_sqr_c_TH_perm);

        Ne = 1/(3*(r_squared_c_TH - mean_r_sqr_c_TH_perm));
        if Ne >= 0
            fprintf('\tNe (mean permuted r^2 as correction factor) = %8.6g\n', Ne);
        else
            fprintf('\tNe (mean permuted r^2 as correction factor) = infinity\n');
        end


        Ne = 1/(3*(r_squared_c_TH - median_r_sqr_c_TH_perm));
        if Ne >= 0
            fprintf('\tNe (median permuted r^2 as correction factor) = %8.6g\n', Ne);
        else
            fprintf('\tNe (median permuted r^2 as correction factor) = infinity\n');
        end

        [lower_CI,upper_CI] = percentile_CI(r_squared_c_TH_permute,n_permutations,alpha_value);
        upper_bound_Ne = 1/(3*(r_squared_c_TH - upper_CI));
        lower_bound_Ne = 1/(3*(r_squared_c_TH - lower_CI));
        if upper_bound_Ne >= 0 
            fprintf('\tNe %g%% confidence intervals:\t%8.6g - %8.6g\n',(1-alpha_value)*100,lower_bound_Ne,upper_bound_Ne);
        else
            fprintf('\tNe %g%% confidence intervals:\t%8.6g - %s\n',(1-alpha_value)*100,lower_bound_Ne,'infinity');
        end

        fprintf('---------------------\n');

        fprintf('permuted r^2_{delta} allele frequency weighted (AFW):\n');
        fprintf('\tmean = %8.6g\n',mean_r_sqr_delta_W_perm);
        fprintf('\tmedian = %8.6g\n\n',median_r_sqr_delta_W_perm);

        Ne = 1/(3*(r_squared_delta_W - mean_r_sqr_delta_W_perm));
        if Ne >= 0
            fprintf('\tNe (mean permuted r^2 as correction factor) = %8.6g\n', Ne);
        else
            fprintf('\tNe (mean permuted r^2 as correction factor) = infinity\n');
        end

        Ne = 1/(3*(r_squared_delta_W - median_r_sqr_delta_W_perm));
        if Ne >= 0
            fprintf('\tNe (median permuted r^2 as correction factor) = %8.6g\n', Ne);
        else
            fprintf('\tNe (median permuted r^2 as correction factor) = infinity\n');
        end

        [lower_CI,upper_CI] = percentile_CI(r_squared_delta_W_permute,n_permutations,alpha_value);
        upper_bound_Ne = 1/(3*(r_squared_delta_W - upper_CI));
        lower_bound_Ne = 1/(3*(r_squared_delta_W - lower_CI));
        if upper_bound_Ne >= 0 
            fprintf('\tNe %g%% confidence intervals:\t%8.6g - %8.6g\n',(1-alpha_value)*100,lower_bound_Ne,upper_bound_Ne);
        else
            fprintf('\tNe %g%% confidence intervals:\t%8.6g - %s\n',(1-alpha_value)*100,lower_bound_Ne,'infinity');
        end

        fprintf('---------------------\n');

        fprintf('permuted r^2_{delta} allele frequency thresholded (AFT):\n');
        fprintf('\tmean = %8.6g\n',mean_r_sqr_delta_TH_perm);
        fprintf('\tmedian = %8.6g\n\n',median_r_sqr_delta_TH_perm); 

        Ne = 1/(3*(r_squared_delta_TH - mean_r_sqr_delta_TH_perm));
        if Ne >= 0
            fprintf('\tNe (mean permuted r^2 as correction factor) = %8.6g\n', Ne);
        else
            fprintf('\tNe (mean permuted r^2 as correction factor) = infinity\n');
        end

        Ne = 1/(3*(r_squared_delta_TH - median_r_sqr_delta_TH_perm));
        if Ne >= 0
            fprintf('\tNe (median permuted r^2 as correction factor) = %8.6g\n', Ne);
        else
            fprintf('\tNe (median permuted r^2 as correction factor) = infinity\n');
        end

        [lower_CI,upper_CI] = percentile_CI(r_squared_delta_TH_permute,n_permutations,alpha_value);
        upper_bound_Ne = 1/(3*(r_squared_delta_TH - upper_CI));
        lower_bound_Ne = 1/(3*(r_squared_delta_TH - lower_CI));
        if upper_bound_Ne >= 0 
            fprintf('\tNe %g%% confidence intervals:\t%8.6g - %8.6g\n',(1-alpha_value)*100,lower_bound_Ne,upper_bound_Ne);
        else
            fprintf('\tNe %g%% confidence intervals:\t%8.6g - %s\n',(1-alpha_value)*100,lower_bound_Ne,'infinity');
        end

        fprintf('==============================================================================\n');
    
    
        % make a table for unbiased estimates (UB) with sample size correction
        fprintf('\n');
        fprintf('==============================================================================\n');
        fprintf('Mean and median r^2 from genotype data permuted among loci. This value is used as\n');
        fprintf('the correction factor in Ne = 1/(3*(r^2*UB - r^2 correction factor)).\n');
        fprintf('Estimates and bounds in this table are corrected for harmonic mean (over loci)\n');
        fprintf('sample sizes of individuals, S, where UB = S/(S-1)) and are therefore unbiased\n');
        fprintf('(Weir, 1979). Estimates are both allele frequency weighted (AFW) and allele\n');
        fprintf('frequency thresholded (AFT).\n');        
        fprintf('------------------------------------------------------------------------------\n');
        fprintf('Permuted genotypes %d times\n\n', n_permutations);
        fprintf('---------------------\n');

        fprintf('permuted r^2_{comp} allele frequency weighted (AFW):\n');
        fprintf('\tmean = %8.6g\n',mean_r_sqr_c_W_perm);
        fprintf('\tmedian = %8.6g\n\n',median_r_sqr_c_W_perm);

        Ne = 1/(3*(r_squared_c_W*UB_weight_comp - mean_r_sqr_c_W_perm)); % Ne estimated with r_squared_c AFW UB mean permute
        if Ne >= 0
            fprintf('\tNe (mean permuted r^2 as correction factor) = %8.6g\n', Ne);
        else
            fprintf('\tNe (mean permuted r^2 as correction factor) = infinity\n');
        end

        Ne = 1/(3*(r_squared_c_W*UB_weight_comp - median_r_sqr_c_W_perm)); % Ne estimated with r_squared_c AFW UB median permute
        if Ne >= 0
            fprintf('\tNe (median permuted r^2 as correction factor) = %8.6g\n', Ne);
        else
            fprintf('\tNe (median permuted r^2 as correction factor) = infinity\n');
        end

        [lower_CI,upper_CI] = percentile_CI(r_squared_c_W_permute,n_permutations,alpha_value);
        
        upper_bound_Ne = 1/(3*(r_squared_c_W*UB_weight_comp - upper_CI));
        lower_bound_Ne = 1/(3*(r_squared_c_W*UB_weight_comp - lower_CI));
        if upper_bound_Ne >= 0 
            fprintf('\tNe %g%% confidence intervals:\t%8.6g - %8.6g\n',(1-alpha_value)*100,lower_bound_Ne,upper_bound_Ne);
        else
            fprintf('\tNe %g%% confidence intervals:\t%8.6g - %s\n',(1-alpha_value)*100,lower_bound_Ne,'infinity');
        end

        fprintf('---------------------\n');

        fprintf('permuted r^2_{comp} allele frequency thresholded (AFT):\n');
        fprintf('\tmean = %8.6g\n',mean_r_sqr_c_TH_perm);
        fprintf('\tmedian = %8.6g\n\n',median_r_sqr_c_TH_perm);

        %Ne = 1/(3*(r_squared_c_TH - mean_r_sqr_c_TH_perm));
        Ne = 1/(3*(r_squared_c_TH*UB_weight_comp - mean_r_sqr_c_TH_perm)); % Ne estimated with r_squared_c AFT UB mean permute

        if Ne >= 0
            fprintf('\tNe (mean permuted r^2 as correction factor) = %8.6g\n', Ne);
        else
            fprintf('\tNe (mean permuted r^2 as correction factor) = infinity\n');
        end


        %Ne = 1/(3*(r_squared_c_TH - median_r_sqr_c_TH_perm));
        Ne = 1/(3*(r_squared_c_TH*UB_weight_comp - median_r_sqr_c_TH_perm)); % Ne estimated with r_squared_c AFT UB median permute

        if Ne >= 0
            fprintf('\tNe (median permuted r^2 as correction factor) = %8.6g\n', Ne);
        else
            fprintf('\tNe (median permuted r^2 as correction factor) = infinity\n');
        end

        [lower_CI,upper_CI] = percentile_CI(r_squared_c_TH_permute,n_permutations,alpha_value);
        upper_bound_Ne = 1/(3*(r_squared_c_TH*UB_weight_comp - upper_CI));
        lower_bound_Ne = 1/(3*(r_squared_c_TH*UB_weight_comp - lower_CI));
        if upper_bound_Ne >= 0 
            fprintf('\tNe %g%% confidence intervals:\t%8.6g - %8.6g\n',(1-alpha_value)*100,lower_bound_Ne,upper_bound_Ne);
        else
            fprintf('\tNe %g%% confidence intervals:\t%8.6g - %s\n',(1-alpha_value)*100,lower_bound_Ne,'infinity');
        end

        fprintf('---------------------\n');

        fprintf('permuted r^2_{delta} allele frequency weighted (AFW):\n');
        fprintf('\tmean = %8.6g\n',mean_r_sqr_delta_W_perm);
        fprintf('\tmedian = %8.6g\n\n',median_r_sqr_delta_W_perm);

        %Ne = 1/(3*(r_squared_delta_W - mean_r_sqr_delta_W_perm));
        Ne = 1/(3*(r_squared_delta_W*UB_weight_delta - mean_r_sqr_c_W_perm)); % Ne estimated with r_squared_delta AFW UB mean permute

        if Ne >= 0
            fprintf('\tNe (mean permuted r^2 as correction factor) = %8.6g\n', Ne);
        else
            fprintf('\tNe (mean permuted r^2 as correction factor) = infinity\n');
        end

        %Ne = 1/(3*(r_squared_delta_W - median_r_sqr_delta_W_perm));
        Ne = 1/(3*(r_squared_delta_W*UB_weight_delta - median_r_sqr_c_W_perm)); % Ne estimated with r_squared_delta AFW UB median permute

        if Ne >= 0
            fprintf('\tNe (median permuted r^2 as correction factor) = %8.6g\n', Ne);
        else
            fprintf('\tNe (median permuted r^2 as correction factor) = infinity\n');
        end

        [lower_CI,upper_CI] = percentile_CI(r_squared_delta_W_permute,n_permutations,alpha_value);
        upper_bound_Ne = 1/(3*(r_squared_delta_W*UB_weight_delta - upper_CI));
        lower_bound_Ne = 1/(3*(r_squared_delta_W*UB_weight_delta - lower_CI));
        if upper_bound_Ne >= 0 
            fprintf('\tNe %g%% confidence intervals:\t%8.6g - %8.6g\n',(1-alpha_value)*100,lower_bound_Ne,upper_bound_Ne);
        else
            fprintf('\tNe %g%% confidence intervals:\t%8.6g - %s\n',(1-alpha_value)*100,lower_bound_Ne,'infinity');
        end

        fprintf('---------------------\n');

        fprintf('permuted r^2_{delta} allele frequency thresholded (AFT):\n');
        fprintf('\tmean = %8.6g\n',mean_r_sqr_delta_TH_perm);
        fprintf('\tmedian = %8.6g\n\n',median_r_sqr_delta_TH_perm); 

        %Ne = 1/(3*(r_squared_delta_TH - mean_r_sqr_delta_TH_perm));
        Ne = 1/(3*(r_squared_delta_TH*UB_weight_delta - mean_r_sqr_c_TH_perm)); % Ne estimated with r_squared_delta AFT UB mean permute

        if Ne >= 0
            fprintf('\tNe (mean permuted r^2 as correction factor) = %8.6g\n', Ne);
        else
            fprintf('\tNe (mean permuted r^2 as correction factor) = infinity\n');
        end

        %Ne = 1/(3*(r_squared_delta_TH - median_r_sqr_delta_TH_perm));
        Ne = 1/(3*(r_squared_delta_TH*UB_weight_delta - median_r_sqr_c_TH_perm)); % Ne estimated with r_squared_delta AFT UB median permute

        if Ne >= 0
            fprintf('\tNe (median permuted r^2 as correction factor) = %8.6g\n', Ne);
        else
            fprintf('\tNe (median permuted r^2 as correction factor) = infinity\n');
        end

        [lower_CI,upper_CI] = percentile_CI(r_squared_delta_TH_permute,n_permutations,alpha_value);
        upper_bound_Ne = 1/(3*(r_squared_delta_TH*UB_weight_delta - upper_CI));
        lower_bound_Ne = 1/(3*(r_squared_delta_TH*UB_weight_delta - lower_CI));
        if upper_bound_Ne >= 0 
            fprintf('\tNe %g%% confidence intervals:\t%8.6g - %8.6g\n',(1-alpha_value)*100,lower_bound_Ne,upper_bound_Ne);
        else
            fprintf('\tNe %g%% confidence intervals:\t%8.6g - %s\n',(1-alpha_value)*100,lower_bound_Ne,'infinity');
        end

        fprintf('==============================================================================\n');
        
    
    end % make_table
    
    % save permuted values to files
    %save('r_squared_c_W_permuted.mat','r_squared_c_W_permute'); % save permuted values
    %save('r_squared_c_TH_permuted.mat','r_squared_c_TH_permute'); % save permuted values
    %save('r_squared_delta_W_permuted.mat','r_squared_delta_W_perm'); % save permuted values
    %save('r_squared_delta_TH_permuted.mat','r_squared_delta_TH_perm'); % save permuted values

end

    %______________________________________
    function [r_squared_c_W_permute,r_squared_c_TH_permute,r_squared_delta_W_permute,r_squared_delta_TH_permute] =... 
            permute(genotypes,n_permutations,num_individuals,num_loci,threshold_allele_freq,make_files)
    
        % Make random permutations of genotypes across all loci in data set, making new random combinations 
        % of genotypes of all individuals. 
        % With each permuted data set, estimate r_squared_comp[permute] and r_squred_delta[permute].
        %
        %   Inputs:
        %   n_permutations - integer number of times to permute data set
        %   genotypes - the genotype data with pairs of columns containing the two integer alleles for a diploid locus
        %   num_individuals - integer number of rows in the genotype data matrix
        %   num_loci - integer half the number of columns in the genotype data matrix
        %   print_notes - boolean switch for sending error fprintf output to console
        %   print_table - boolean switch for output of table of pairwise estimates
        %
        %   Outputs:
        %   r_squared_c_perm - a num_individuals by 1 vector of float r^2_c values,
        %       each value is the estimate from one permuted data set.
        %   r_squared_delta_perm - a num_individuals by 1 vector of float r^2_delta values,
        %       each value is the estimate from one permuted data set.

        % suppress any output from function make_estimate
        print_notes = false; % do not print notes
        print_table = false; % do not show table
        print_locus_table = false; % do not show table        
        
        r_squared_c_W_permute = zeros(n_permutations,1); %allocate space
        r_squared_c_TH_permute = zeros(n_permutations,1); %allocate space
        r_squared_delta_W_permute = zeros(n_permutations,1); %allocate space
        r_squared_delta_TH_permute = zeros(n_permutations,1); %allocate space

        %fprintf('\nPermuting genotypes %d times...\n', n_permutations);

        h = waitbar(0,'permuting genotype data ...');
        
        for i=1:n_permutations

            waitbar(i/n_permutations);

            [perm_genotypes] = permute_genotypes(genotypes,num_individuals,num_loci);   % permute genotype data

            % estimate the two r^2 measures (both with AFW and AFT) with permuted genotypes
            [r_squared_c_W_permute(i,1),r_squared_c_TH_permute(i,1),r_squared_delta_W_permute(i,1),r_squared_delta_TH_permute(i,1),~,~,~,~,~,~] = ...
                make_estimate(print_notes,perm_genotypes,num_individuals,num_loci,print_table,threshold_allele_freq,print_locus_table,make_files);
        
        end %for i=1:n_permutations

        close(h);

    end

    %______________________________________    
    % function [permuted_genotypes] = permute_alleles(genotypes,num_individuals,num_loci)
    % % Generate a random permutation of a genotype matrix.
    % % The permutation is within each locus. This scrambles allele pairings in genotypes 
    % % and makes new random combinations of alleles into genotypes.
    % %
    % % Inputs:   
    % %   genotypes - matrix of genotype data, individuals in rows and loci
    % %               in columns
    % %	num_individuals - integer number of rows in genotypes
    % %	num_loci - integer half the number of columns in genotypes
    % %
    % % Outputs:  
    % %   permuted_genotypes is the random permutation of the genotype matrix
    % 
    %     for i = 1:2*num_loci
    % 
    %         permutation = randperm(num_individuals);    % get random permutation, samples without replacement
    %         index_col = permutation'; % transpose to get column vector
    % 
    %         permuted_genotypes(:,i) = genotypes(index_col,i);   % assign random permutation of genotype data in column i to permute_genotypes
    %     end


    %______________________________________    
    function [permuted_genotypes] = permute_genotypes(genotypes,num_individuals,num_loci)
    % Generate a random permutation of a genotype matrix.
    % The permutation is among pairs of columns. This scrambles multilocus genotypes (keeping each single locus genotype intact)
    % and makes new random combinations of loci into genotypes.
    %
    % Inputs:   
    %   genotypes - matrix of genotype data, individuals in rows and loci
    %               in columns
    %	num_individuals - integer number of rows in genotypes
    %	num_loci - integer half the number of columns in genotypes
    %
    % Outputs:  
    %   permuted_genotypes is the random permutation of the genotype matrix

        permuted_genotypes = zeros(num_individuals,2*num_loci);

        for i = 1:num_loci

            permutation = randperm(num_individuals);    % get random permutation, samples without replacement
            rows = permutation'; % transpose to get column vector

            %paired_rows = repmat(rows,[1,2])

            col = 2*i-1; % the first column of the locus i in the genotype matrix

            permuted_genotypes(:,col:col+1) = genotypes(rows,col:col+1);   % assign random rows of genotype data to rows in permuted_genotypes
        end

    end
