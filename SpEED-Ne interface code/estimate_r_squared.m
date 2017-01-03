function [r_squared_comp_weighted,r_squared_comp_thresholded,r_squared_delta_weighted,r_squared_delta_thresholded,r_delta_div_zero,output_table] = ...
    estimate_r_squared(c_hap_table,k1,k2,freqs1,freqs2,S,paa,pbb,locus_one,locus_two,print_notes,threshold_allele_freq,print_table,make_files)
%   Function to compute r_squared_comp and r_squared_delta from a composite
%   haplotype table. 
%
% version 1.0 27 Sept 2016
%
%
%	The function requires a composite haplotype table (page 2 in Sved et al. 2013) for each pair of loci. 
%	The overall estimates of r_squared_comp and r_squared_delta for each pair of loci are based
%	on allele frequency weighting (AFW) as well as the mean over all pairs
%	of alleles with both allele frequencies above a specificed threshold
%	(thresholding or AFT). 
%
%   Inputs:
%   c_hap_table - k1 by k2 matrix of integers that make the composite haplotype table 
%   k1 - integer number of alleles for locus one
%   k2 - integer number of alleles for locus two
%   freqs1 - num_loci by one vector of allele frequencies for locus one
%   freqs2 - num_loci by one vector of allele frequencies for locus two
%   S - integer sample size of individuals
%   states1 - k1 by 1 vector of allelic states for locus one
%   states2 - k2 by 1 vector of allelic states for locus two
%   paa - k1 by 1 vector of the observed frequency of homozygous genotypes at locus one
%   pbb - k2 by 1 vector of the observed frequency of homozygous genotypes at locus two
%   locus_one - integer index of locus one used for error output
%	locus_two - integer index of locus two used for error output
%   print_notes - boolean switch for sending error fprintf output to console
%   threshold_allele_freq - real number threshold of allele frequency to use in
%       estimates; estimates involving alleles with frequency < threshold
%       are not used
%   print_table - boolean switch for output of table of pairwise estimates
%   make_files - boolean switch for output of results to files - requires
%       prep of results tables
%
%   Outputs:
%   r_squared_comp_weighted - r^2 composite genetic association measure with each allele pair estimate weighted by allele frequencies
%   r_squared_comp_thresholded - r^2 composite genetic association measure with each allele pair estimate subject to allele frequencies threshold
%   r_squared_delta_weighted - r^2 delta genetic association measure adjusted for within locus disequilibrium with each allele pair estimate weighted by allele frequencies
%   r_squared_delta_thresholded - r^2 delta genetic association measure adjusted for within locus disequilibrium with each allele pair estimate subject to allele frequencies threshold
%   r_delta_div_zero - boolean to flag returned true when the demoninator of r_ij_delta was zero.
%   output_table - matrix to store all pairwise estimates for later output to file
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


    Dcomp = zeros(k1,k2); %initialize, rows = number of alleles for locus with fewer alleles

    sum_pipj = 0;   % initialize
    
    counter = 0; % initialize
    
    for loc_one_allele = 1:k1        
        for loc_two_allele = 1:k2
            
            counter = counter + 1; % iterate
            
        	M = c_hap_table(loc_one_allele,loc_two_allele);

            two_na = sum(c_hap_table(loc_one_allele,:)); % sum all columns
            two_nb = sum(c_hap_table(:,loc_two_allele)); % sum all rows
            
            Dcomp(loc_one_allele,loc_two_allele) = (M/(4*S)) - (two_na/(4*S))*(two_nb/(4*S));
            
            % compute r_squared_lm from equation 4
            % Numerator multipled by 4 since r^2_c = 4D^2(comp)/pa(1-pa)pb(1-pb) in equation 4.             
            r_comp_denominator = freqs1(loc_one_allele,1)*(1 - freqs1(loc_one_allele,1))*freqs2(loc_two_allele,1)*(1 - freqs2(loc_two_allele,1));
            r_ij_comp = (4*Dcomp(loc_one_allele,loc_two_allele)^2) / r_comp_denominator;
                        
                        
            % compute r_squared_delta denominator terms for equation 5.
            term_one = freqs1(loc_one_allele,1)*(1 - freqs1(loc_one_allele,1));
            term_two = paa(loc_one_allele,1) - freqs1(loc_one_allele,1)^2;
            term_three = freqs2(loc_two_allele,1)*(1 - freqs2(loc_two_allele,1));            
            term_four = pbb(loc_two_allele,1) - freqs2(loc_two_allele,1)^2;

            delta_denom = ((term_one + term_two)*(term_three + term_four));
            
            if delta_denom == 0
                if print_notes
                    fprintf('r_squared_delta_denominator equal zero in function estimate_r_squared\n\t(loc_one_allele = %d and loc_two_allele = %d of loci %d and %d)\n',loc_one_allele,loc_two_allele,locus_one,locus_two);
                end
                r_ij_delta = 0; % force value of zero, avoid divide by zero
                r_delta_div_zero = true;   %set boolean
            else
                r_ij_delta = (4*Dcomp(loc_one_allele,loc_two_allele)^2) / delta_denom;
                r_delta_div_zero = false;   %set boolean
            end
            
            
            % store results from all pairs of alleles
            estimates(counter,1) = r_ij_comp; % store estimate of r_squared_comp from one pair of alleles
            estimates(counter,2) = r_ij_delta; % store estimate of r_squared_delta from one pair of alleles
            estimates(counter,3) = freqs1(loc_one_allele,1); % store freq of allele at locus 1
            estimates(counter,4) = freqs2(loc_two_allele,1); % store freq of allele at locus 2
            
            estimates(counter,5) = freqs1(loc_one_allele,1)*freqs2(loc_two_allele,1); % store product piqj of allele freqs
            
            
            if (freqs1(loc_one_allele,1) < threshold_allele_freq) || (freqs2(loc_two_allele,1) < threshold_allele_freq)
                estimates(counter,6) = 0; % mark estimate as coming from at least one allele with freq below threshold
            else
                estimates(counter,6) = 1; % mark estimate as coming from both alleles with freqs >= threshold
            end
            
            estimates(counter,7) = loc_one_allele;  % allelic state at locus one
            estimates(counter,8) = loc_two_allele;  % allelic state at locus two
            
        end %for loc_two_allele = 1:k2        
	end %loc_one_allele = 1:k2
    
    %sum_pipj = sum(estimates(:,5));
    %fprintf('debug: %f\n', sum_pipj);

    
    % Determine allele frequency weighted (AFW) estimates
    % r_squared_comp and r-squared_delta weighted by allele frequencies as in equation 7 of Sved et al 2013. 
    % Weighting piqj values are stored in estimates(:,5). 
    
    r_squared_lm_comp_weighted = (estimates(:,5).*estimates(:,1));
    r_squared_comp_weighted = sum(r_squared_lm_comp_weighted); % sum down column
    
    r_squared_lm_delta_weighted = (estimates(:,5).*estimates(:,2));
    r_squared_delta_weighted = sum(r_squared_lm_delta_weighted); % sum down column

    
    % Determine allele frequency thresholded (AFT) estimates
    % Estimates based on any alleles below threshold removed through multication by indicator value in estimates(:,6).
    threshold_n = sum(estimates(:,6)); % number of r_squared estimates above threshold
    
    r_squared_comp_thresholded = sum(estimates(:,6).*estimates(:,1))/threshold_n; % sum down column, take average
    
    r_squared_delta_thresholded = sum(estimates(:,6).*estimates(:,2))/threshold_n; % sum down column, take average
    

    if print_table
        % make output table for estimates from all allele pairings

        fprintf('==========================================================================================\n');
        fprintf('Table of r^2 values for all allele pairings for one pair of loci.\n');
        fprintf('Values employing allele frequency weighting (AFW):\n\n');
        fprintf('k(1)\tk(2)\tallele freq 1\tallele freq 2\tweight\t\tr^2_{comp}\tr^2_{delta}\n');
        fprintf('------------------------------------------------------------------------------------------\n');

        for row=1:counter

            fprintf('%2d\t',estimates(row,7));
            fprintf('%2d\t',estimates(row,8));

            fprintf('%3.6f\t',estimates(row,3));
            fprintf('%3.6f\t',estimates(row,4));

            fprintf('%3.6f\t',estimates(row,5));

            fprintf('%3.6f\t',estimates(row,1));
            fprintf('%3.6f\t',estimates(row,2));
            
            fprintf('\n');

        end %for row
        fprintf('------------------------------------------------------------------------------------------\n');
        fprintf('sum weighted r^2_{comp}  = %6.6f\n', r_squared_comp_weighted);
        fprintf('sum weighted r^2_{delta} = %6.6f\n',r_squared_delta_weighted);
        fprintf('==========================================================================================\n\n\n');
      
        % Print table for estimates using allele frequency threshold in average for all loci.
        fprintf('==========================================================================================\n');
        fprintf('Table of r^2 values for all allele pairings for one pair of loci.\n');
        fprintf('Values for allele frequency threshold (AFT) of: %3.6f\n\n', threshold_allele_freq);
        fprintf('k(1)\tk(2)\tallele freq 1\tallele freq 2\tthreshold\tr^2_{comp}\tr^2_{delta}\n');
        fprintf('------------------------------------------------------------------------------------------\n');

        for row=1:counter
            fprintf('%2d\t',estimates(row,7));
            fprintf('%2d\t',estimates(row,8));

            fprintf('%3.6f\t',estimates(row,3));
            fprintf('%3.6f\t',estimates(row,4));
            
            if estimates(row,6) == 1
                fprintf('Above\t\t');
            else
                fprintf('Below\t\t');
            end
            
            fprintf('%3.6f\t',estimates(row,1));
            fprintf('%3.6f\t',estimates(row,2));
            
            fprintf('\n');

        end %for row
        
        fprintf('------------------------------------------------------------------------------------------\n');
        fprintf('sum thresholded r^2_{comp}  = %6.6f\n', r_squared_comp_thresholded);
        fprintf('sum thresholded r^2_{delta} = %6.6f\n',r_squared_delta_thresholded);
        fprintf('==========================================================================================\n\n\n');

    end % if print_table
    
    
    if make_files
        % make table of all pairwise estimates for file output
        output_table = zeros(counter,10); % allocate space
        %output_table(1,:) = ['locus 1','locus 2','allele 1','allele 2','freq allele 1','freq allele 2','above threhold','product allele freqs','r^2{c}','r^2{delta}'];

        for row=1:counter
            output_table(row,1) = locus_one; % locus one
            output_table(row,2) = locus_two; % locus two

            output_table(row,3) = estimates(row,7); % allele one
            output_table(row,4) = estimates(row,8); % allele two
            
            output_table(row,5) = estimates(row,3); % freq of allele at locus 1/allele 1
            output_table(row,6) = estimates(row,4); % freq of allele at locus 2/allele 2
            
            output_table(row,7) = estimates(row,6); %  0 = at least one allele with freq below threshold, 1 = both above threshold

            output_table(row,8) = estimates(row,5);  % product piqj of allele freqs

            output_table(row,9) = estimates(row,1); % r_squared_comp
            output_table(row,10) = estimates(row,2); % r_squared_delta
            
        end; %for 1:counter
    else
        output_table = 0; % assign zero so something returned
    end; % if make_files

    