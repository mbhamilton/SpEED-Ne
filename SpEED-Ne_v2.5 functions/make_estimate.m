function [r_squared_c_W,r_squared_c_TH,r_squared_delta_W,r_squared_delta_TH,r_squared_c_values,r_squared_delta_values,...
    S_c_values,S_delta_values,c_locus_pairs,delta_locus_pairs,all_locus_pairs_table,all_allele_pairs_table] = ...
    make_estimate(print_notes,genotypes,num_individuals,num_loci,print_table,threshold_allele_freq,print_locus_table,make_files)
% Estimate r^2_{c} and r^2_{delta} and also tally sample sizes while accounting for missing data and possible cases of div
% zero. 
%
% version 1.0 27 Sept 2016
%
%   Inputs:
%   print_notes - boolean switch for sending error fprintf output to console
%   genotypes - the genetic data with pairs of columns containing the two
%   	integer alleles for a diploid locus
%       -> all missing alleles should be coded as NaN
%   num_individuals - integer number of rows in the genotype data matrix
%   num_loci - integer half the number of columns in the genotype data matrix
%   print_table - boolean switch for output table of r^2 estimates for all pairs of loci
%   threshold_allele_freq - real number threshold of allele frequency to use in
%       estimates; estimates involving alleles with frequency < threshold
%       are not used
%   print_locus_table - boolean switch for output of table of estimates for
%       all pairs of alleles within each locus 
%   make_files - boolean switch for output of results to files - requires
%       prep of results tables
%
%   Outputs:
%   r_squared_c_W - estimate of r^2_{c} with weighting by allele frequency
%   r_squared_c_TH - estimate of r^2_{c} with allele frequency thresholding
%   r_squared_delta_W - estimate of r^2_{delta} with weighting by allele frequency
%   r_squared_delta_TH - estimate of r^2_{delta} with allele frequency thresholding

%   r_squared_c_values - c_pairs by 3 vector of r_squared_c for each locus pair, weighted by  alleles per locus according to Sved et al (2103) eqn. 6 and 7. 
%       First column = numerator for weighted value according to Sved et al 2013 equation 7
%       second column = value for threshold allele frequency value, 
%       third column = denominator.
%   r_squared_delta_values - delta_pairs by 3 vector of r_squared_c for each locus pair, weighted by  alleles per locus according to Sved et al (2103) eqn. 6 and 7. 
%       First column = numerator for weighted value according to Sved et al 2013 equation 7
%       second column = value for threshold allele frequency value, 
%       third column = denominator.
%   S_c_values - c_pairs by 1 vector of integer values of number of
%       individuals for each locus pair used to estimate r^2_{c}
%   S_delta_values - c_pairs by 1 vector of integer values of number of
%       individuals for each locus pair used to estimate r^2_{delta}
%   c_locus_pairs - sample size of pairs of loci for r_squared_c estimates
%   delta_locus_pairs  - sample size of pairs of loci for estimates of r_squared_delta
%   all_locus_pairs_table - matrix of all locus pairs r^2 estimates for later output to file
%   all_allele_pairs_table - matrix of all allele pairs r^2 estimates for later output to file
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


    % initialize variables
    sum_r_squared_c_numerator_W = 0;
    sum_r_squared_c_numerator_TH = 0;
    
    sum_r_squared_delta_numerator_W = 0;
    sum_r_squared_delta_numerator_TH = 0;
    
    sum_r_squared_c_denominator = 0;
    sum_r_squared_delta_denominator = 0;

    c_locus_pairs = 0;
    delta_locus_pairs = 0;
    
    table_row = 0;
    
    locus_pairs = (num_loci*(num_loci-1))/2;
    
    % allocate space
    r_squared_c_W = zeros(num_loci,num_loci); %allocate space
    r_squared_c_TH = zeros(num_loci,num_loci); %allocate space
    
    r_squared_delta_W = zeros(num_loci,num_loci); %allocate space
    r_squared_delta_TH = zeros(num_loci,num_loci); %allocate space
    
    %r_squared_table = zeros(locus_pairs,10); %allocate space, will not work for very large numbers of loci
    
    r_squared_table = zeros(1,10); %allocate space for one pair of loci at a time
    
    %all_allele_pairs_table = zeros(1,10);
    
    %r_squared_table_part = zeros(10,num_loci); %allocate space for part of the pairs of loci, avoids memory limits when there are many loci
        
    % the path of the location to save scratch files for datastore
    %savePath = './Scratchdata/'; 
    
    %stubname = 'scratchfile_'; % stub of file name for scratch files
    
    % determine numbers of alleles per locus, allele frequencies at each locus, and allelic states
    [num_alleles, ~, ~] = alleles(genotypes,num_individuals,num_loci);
    
    
    S_c_values = zeros(locus_pairs,1); %allocate space for all locus_pairs, shrink after loop if needed
    
    r_squared_c_values = zeros(locus_pairs,3); %allocate space for all locus_pairs, shrink after loop if needed
    
    S_delta_values = zeros(locus_pairs,1); %allocate space for all locus_pairs, shrink after loop if needed
    
    r_squared_delta_values = zeros(locus_pairs,3); %allocate space for all locus_pairs, shrink after loop if needed

    
    table_row = 1;
        
    %fprintf('estimating for all pairs of loci ...\n');

    % iterate across two locus haplotypes for a pair of loci
    for i=1:num_loci
        
        %fprintf('estimating all pairs of loci %d\n', i);
        
        %r_squared_table_part(:) = 0; % initialize by assigning zeros to all entries in matrix
        
        %fprintf('in make_estimate: i = %d \n', i); % debugging trace
        
        for j=i+1:num_loci            
            
            %table_row = table_row + 1;            
            
            
            if num_alleles(1,i) <= 1    %locus is monomorphic
%                 if print_notes
%                     fprintf('NOTE: locus %d has only one allele and was skipped!\n', i);
%                 end

                r_squared_table(table_row,1) = i;
                r_squared_table(table_row,2) = j;
                r_squared_table(table_row,3) = NaN;
                r_squared_table(table_row,4) = num_alleles(1,i);
                r_squared_table(table_row,5) = num_alleles(1,j);
            
                continue;  %do not execute this iteration of the for ith loop
            end

            if num_alleles(1,j) <= 1    %locus is monomorphic
%                 if print_notes
%                     fprintf('NOTE: locus %d has only one allele and was skipped!\n', j);
%                 end
                
                r_squared_table(table_row,1) = i;
                r_squared_table(table_row,2) = j;
                r_squared_table(table_row,3) = NaN;
                r_squared_table(table_row,4) = num_alleles(1,i);
                r_squared_table(table_row,5) = num_alleles(1,j);
                
                continue;  %do not execute this iteration of the for jth loop
            end

            col1 = (i*2); %second column of genotype table containing genotypes of locus i
            col2 = (j*2); %second column of genotype table containing genotypes of locus j

            G1 = genotypes(:,col1-1:col1); % copy out two columns of genotype data for locus i
            G2 = genotypes(:,col2-1:col2); % copy out two columns of genotype data for locus j

            % make data set for pair of loci (four columns)
            paired_genotypes = cat(2,G1,G2);

            % determine numbers of alleles per locus, allele frequencies at each locus, and allelic states for pair of loci
            [p_num_alleles, p_allele_freqs, p_allele_states] = alleles(paired_genotypes,num_individuals,2);

            if p_num_alleles(1,1) <= 1 % locus is monomorphic 
                if print_notes
                    fprintf('NOTE: locus %d has only one allele and the locus pair %d-%d was skipped!\n', i,i,j);
                end
                continue;  %do not execute this iteration of the for jth loop
            elseif p_num_alleles(1,2) <= 1 % locus is monomorphic
                if print_notes
                    fprintf('NOTE: locus %d has only one allele and the locus pair %d-%d was skipped!\n', j,i,j);
                end
                continue;  %do not execute this iteration of the for jth loop
            end

            [c_hap_table, S, paa, pbb] = comp_hap_table(paired_genotypes(:,1:2),paired_genotypes(:,3:4),p_num_alleles(1,1),p_num_alleles(1,2),p_allele_states(:,1),p_allele_states(:,2));
            
            %print_locus_table = true;
            
            % estimate  r_squared_c (equation 4) and r_squared_delta for each pair of loci, weighting by allele frequency
            if print_locus_table
                fprintf('********** Locus %d - Locus %d pair **********\n',i,j);
            end
                        
            
            [r_squared_c_W(i,j),r_squared_c_TH(i,j),r_squared_delta_W(i,j),r_squared_delta_TH(i,j),r_delta_div_zero,allele_pairs_table] = estimate_r_squared(c_hap_table,p_num_alleles(1,1),...
                p_num_alleles(1,2),p_allele_freqs(:,1),p_allele_freqs(:,2),S,paa,pbb,i,j,print_notes,threshold_allele_freq,print_locus_table,make_files);            
            % In values returned from estimate_r_squared function "_W" stands for allele-frequency weighted for each pair of alleles,
            % "_TH" stands for allele frequency thresholded where estimates from any allele pair with an allele below the frequency threshold is dropped.
            
            
            %%%%
            % compute components of r_squared_c_lm over all pairs of loci for equation 6
            % sum weighted r_lm_squared_(comp) over all pairs of loci according to equation 6
            r_squared_c_numerator_W = r_squared_c_W(i,j)*(S^2)*(p_num_alleles(1,1) - 1)*(p_num_alleles(1,2) - 1);
            r_squared_c_numerator_TH = r_squared_c_TH(i,j)*(S^2)*(p_num_alleles(1,1) - 1)*(p_num_alleles(1,2) - 1);
            r_squared_c_denominator = (S^2)*(p_num_alleles(1,1) - 1)*(p_num_alleles(1,2) - 1);

            if r_squared_c_denominator == 0
                fprintf('r_squared_c_denominator equal zero for loci i = %d and j = %d in function make_estimate\n',i,j);
            else
            	% sum r_squared_c over all pairs of loci according to equation 6
                sum_r_squared_c_numerator_W = sum_r_squared_c_numerator_W + r_squared_c_numerator_W;
                sum_r_squared_c_numerator_TH = sum_r_squared_c_numerator_TH + r_squared_c_numerator_TH;
                sum_r_squared_c_denominator = sum_r_squared_c_denominator + r_squared_c_denominator;
                
                c_locus_pairs = c_locus_pairs + 1; % counter
                
                % store sample size of individuals for one locus pair
                S_c_values(c_locus_pairs,1) = S;
                
                % store result for one locus pair 
                r_squared_c_values(c_locus_pairs,1) = r_squared_c_numerator_W;
                r_squared_c_values(c_locus_pairs,2) = r_squared_c_numerator_TH;
                r_squared_c_values(c_locus_pairs,3) = r_squared_c_denominator;

            end %if r_squared_c_denominator == 0

            
            if r_delta_div_zero == false
                % Compute components of r_squared_delta_lm over all pairs of loci for equation 6 in Sved et al.
                % Skip cases when denominator of r_squared_delta is zero in estimate_r_squared function;
                % avoids bias of forcing r_squared_delta = 0 or some arbitrary value to prevent div zero.
                
                r_squared_delta_numerator_W = r_squared_delta_W(i,j)*(S^2)*(p_num_alleles(1,1) - 1)*(p_num_alleles(1,2) - 1);
                r_squared_delta_numerator_TH = r_squared_delta_TH(i,j)*(S^2)*(p_num_alleles(1,1) - 1)*(p_num_alleles(1,2) - 1);

                if isnan(r_squared_delta_numerator_TH) || isnan(r_squared_delta_numerator_W)
                    fprintf('r_squared_delta_numerator W or TH equal NaN for loci i = %d and j = %d in function make_estimate\n',i,j);
                end

                % sum r_squared_delta over all pairs of loci according to equation 6
                sum_r_squared_delta_numerator_W = sum_r_squared_delta_numerator_W + r_squared_delta_numerator_W;
                sum_r_squared_delta_numerator_TH = sum_r_squared_delta_numerator_TH + r_squared_delta_numerator_TH;

                r_squared_delta_denominator = (S^2)*(p_num_alleles(1,1) - 1)*(p_num_alleles(1,2) - 1);

                if (r_squared_delta_denominator == 0)
                    fprintf('r_squared_delta_denominator equal zero for loci i = %d and j = %d in function make_estimate\n',i,j);
                end

                sum_r_squared_delta_denominator = sum_r_squared_delta_denominator + r_squared_delta_denominator;
                
                delta_locus_pairs = delta_locus_pairs + 1;
                
                % store sample size of individuals for one locus pair
                S_delta_values(delta_locus_pairs,1) = S;

                % store result for one locus pair 
                r_squared_delta_values(delta_locus_pairs,1) = r_squared_delta_numerator_W;
                r_squared_delta_values(delta_locus_pairs,2) = r_squared_delta_numerator_TH;
                r_squared_delta_values(delta_locus_pairs,3) = r_squared_delta_denominator;
                
                % save values for table
                r_squared_table(table_row,1) = i;
                r_squared_table(table_row,2) = j;
                r_squared_table(table_row,3) = S;
                r_squared_table(table_row,4) = p_num_alleles(1,1);
                r_squared_table(table_row,5) = p_num_alleles(1,2);
                r_squared_table(table_row,6) = r_squared_delta_numerator_TH;
                r_squared_table(table_row,7) = r_squared_delta_denominator;
                
            end %if r_delta_div_zero == false
            
            if i==1 && j==2 % && make_files 
                all_allele_pairs_table = allele_pairs_table; % assign first table as place holder so something is passed back
            end
            
            if make_files && i~=1 && j~=2
                all_allele_pairs_table = vertcat(all_allele_pairs_table,allele_pairs_table); % add new table of allele pairs to tables from other pairs
            end
            
%             if make_files 
%                 if i==1 && j==2
%                     all_allele_pairs_table = allele_pairs_table;
%                 else
%                     all_allele_pairs_table = vertcat(all_allele_pairs_table,allele_pairs_table); % add new table of allele pairs to tables from other pairs
%                 end
%             else
%                 all_allele_pairs_table = allele_pairs_table; % assign first table as place holder so something is passed back
%             end

        end %for j=i+1:num_loci
        
        
    end %for i=1:num_loci
    
    % compute Sved et al. equation 6
    r_squared_c_W = sum_r_squared_c_numerator_W/sum_r_squared_c_denominator;
    r_squared_c_TH = sum_r_squared_c_numerator_TH/sum_r_squared_c_denominator;
    
    r_squared_delta_W = sum_r_squared_delta_numerator_W/sum_r_squared_delta_denominator;
    r_squared_delta_TH = sum_r_squared_delta_numerator_TH/sum_r_squared_delta_denominator;
    
    
    % check on realized dimensions of pre-allocated array
    if c_locus_pairs ~= locus_pairs
        
        %remove trailing rows initially allocated but not used
        S_c_values = S_c_values(1:c_locus_pairs, :);
        
        r_squared_c_values = r_squared_c_values(1:c_locus_pairs, :);
        
    end
    
    % check on realized dimensions of pre-allocated array
    if delta_locus_pairs ~= locus_pairs
        
        %remove trailing rows initially allocated but not used
        S_delta_values = S_delta_values(1:delta_locus_pairs, :);
        
        r_squared_delta_values = r_squared_delta_values(1:delta_locus_pairs, :);
        
    end

    
    
    
    
%     % Cannot output pairwise estimates to file nor make pairwise estimates table since
%     % all_locus_pairs_table does not store estimates in this version!!!
%     
    if make_files
        % make table of all pairwise estimates for file output
        all_locus_pairs_table = zeros(locus_pairs,11); % allocate space
        %output_table(1,:) = ['locus ','locus 2','sample of individuals','k1','k2','r^2{c} weight','r^2{delta} weight','r^2{c} AFW','r^2{delta} AFW','r^2{c} AFT','r^2{delta} AFT'];

        for row=1:locus_pairs
            all_locus_pairs_table(row,1) = r_squared_table(row,1); % locus one
            all_locus_pairs_table(row,2) = r_squared_table(row,2); % locus two

            all_locus_pairs_table(row,3) = r_squared_table(row,3); % sample size of individuals

            all_locus_pairs_table(row,4) = r_squared_table(row,4); % allele number at locus one
            all_locus_pairs_table(row,5) = r_squared_table(row,5); % allele number at locus two

            c_weight = (S^2*(r_squared_table(row,4) - 1)*(r_squared_table(row,5) - 1))/sum_r_squared_c_denominator;
            delta_weight = (S^2*(r_squared_table(row,4) - 1)*(r_squared_table(row,5) - 1))/sum_r_squared_delta_denominator;
            all_locus_pairs_table(row,6) = c_weight; % c_weight
            all_locus_pairs_table(row,7) = delta_weight; % delta_weight

            all_locus_pairs_table(row,8) = (r_squared_c_values(row,1)/r_squared_c_values(row,3))*c_weight; % AFW numerator/denominator * ML weight
            all_locus_pairs_table(row,9) = (r_squared_delta_values(row,1)/r_squared_c_values(row,3))*delta_weight; % AFW numerator/denominator * ML weight

            all_locus_pairs_table(row,10) = (r_squared_c_values(row,2)/r_squared_c_values(row,3))*c_weight; % AFT numerator/denomiator * ML weight
            all_locus_pairs_table(row,11) = (r_squared_delta_values(row,2)/r_squared_c_values(row,3))*delta_weight; % AFT numerator/denomiator * ML weight

        end %for row
    else
        all_locus_pairs_table = 0; % assign zero so something returned
    end % if make_files
    
    
    if print_table
        % make output table of r_squared values for each locus pair
        fprintf('\n');
        fprintf('======================================================================================================================\n');
        fprintf('r_squared values for each locus pair with multilocus weights based on numbers of individuals and numbers \n');
        fprintf('of alleles (eqn. 6 in Sved et al. 2013). For each locus pair estimates are either allele frequency weighted (AFW)\n');
        fprintf('or allele frequency thresholded (AFT).\n\n');        
        fprintf('Locus pair\tS\tk(1)\tk(2)\tmultilocus\tAFW, ML weighted\t\tAFT, ML weighted\n');
        fprintf('\t\t\t\t\tweighting\tr^2_{comp}\tr^2_{delta}\tr^2_{comp}\tr^2_{delta}\n');
        fprintf('----------------------------------------------------------------------------------------------------------------------\n');

        for row=1:locus_pairs

            fprintf('%3d,%d\t',r_squared_table(row,1),r_squared_table(row,2)); % pair of loci
            fprintf('%12.2f\t',r_squared_table(row,3)); % sample size of individuals
            fprintf('%3d\t%3d\t',r_squared_table(row,4),r_squared_table(row,5)); % number of alleles at the loci
           
            c_weight = (S^2*(r_squared_table(row,4) - 1)*(r_squared_table(row,5) - 1))/sum_r_squared_c_denominator;
            delta_weight = (S^2*(r_squared_table(row,4) - 1)*(r_squared_table(row,5) - 1))/sum_r_squared_delta_denominator;
            fprintf('%6.6f\t',c_weight);
            %sum_weights = sum_weights + weight;
            
            fprintf('%6.6f\t',(r_squared_c_values(row,1)/r_squared_c_values(row,3))*c_weight); % W numerator/denomiator * ML weight
            fprintf('%6.6f\t',(r_squared_delta_values(row,1)/r_squared_c_values(row,3))*delta_weight); % W numerator/denomiator * ML weight
            
            
            fprintf('%6.6f\t',(r_squared_c_values(row,2)/r_squared_c_values(row,3))*c_weight); % TH numerator/denomiator * ML weight
            fprintf('%6.6f\t',(r_squared_delta_values(row,2)/r_squared_c_values(row,3))*delta_weight); % TH numerator/denomiator * ML weight

            fprintf('\n');

        end %for row
        
        fprintf('----------------------------------------------------------------------------------------------------------------------\n');
        fprintf('Estimates over all locus pairs:\n\n');
        
        fprintf('\tAllele frequency weighted (AFW) within locus pairs:\n');
        fprintf('\t\tr^2_{comp} = %6.6f\n',r_squared_c_W);
        fprintf('\t\tr^2_{delta} = %6.6f\n',r_squared_delta_W);
        
        fprintf('\n');
        
        fprintf('\tAllele frequency thresholded (AFT) within locus pairs:\n');
        fprintf('\t\tr^2_{comp} = %6.6f\n', r_squared_c_TH);
        fprintf('\t\tr^2_{delta} = %6.6f\n', r_squared_delta_TH);

        fprintf('======================================================================================================================\n');
        fprintf('\n');
        
    end % if print_table
    
