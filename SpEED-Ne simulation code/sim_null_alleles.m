function [mod_genotypes] = sim_null_alleles(genotypes,num_individs,num_loci,max_freq_null,num_alleles,orig_allele_freqs,allele_states)
% Function to simulate null alleles in diploid genotype data. Given a data set of genotypes, the function
% will sample one allele at random with a frequency ? max_freq_null at each
% locus. All heterozygous genotypes containing the null allele are altered to be homozygous for 
% the non-null allele and all genotypes homozygous for the null allele were re-coded as missing data.
% If a locus has no alleles ? pnull the locus is skipped.
%
% This version simulates genotype data for each iteration using generate_pop.m. 
%
% version 1.0 10 Oct 2016
%
% Inputs:
%   genotypes - number of individuals by 2*num_loci matrix of integers of allelic states
%   num_individs - number of rows in genotypes
%   num_loci - integer half the number of columns in the genotype data matrix
%   max_freq_null - maximum alle frequency of allele to convert to a null allele
%   num_alleles - one by num_loci vector of the integer number of alleles per locus
%   allele_freqs - max num_alleles by num_loci vector of allele frequencies
%                  (entries without an allele contain zeros)
%   allele_states - max num_alleles by num_loci vector of integer allelic
%                   states (entries without an allele contain zeros)
% Outputs:
%   mod_genotypes - number of individuals by 2*num_loci matrix of genotypes
%       with some heterozygous genotypes modified to exhibit null alleles
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
%   A copy of the GNU General Public License is available to http://www.gnu.org/licenses/.
%**************


    mod_genotypes = genotypes; % initialize
    
    % determine alleles at each locus <= max_freq_null
    less_than_max = (orig_allele_freqs <= max_freq_null);
    greater_than_zero = (0 < orig_allele_freqs);
    
    poss_nulls = less_than_max.*greater_than_zero;
    
    total_poss_nulls = sum(poss_nulls,1); % sum columns to determine number of alleles at each locus
    
    for i=1:num_loci
        if total_poss_nulls(1,i) > 0
            
            allele_row = randi([1,num_alleles(1,i)],1,1); % pick random allele
            freq = orig_allele_freqs(allele_row,i); % get frequency of random allele
            
            while freq > max_freq_null
            	allele_row = randi([1,num_alleles(1,i)],1,1); % pick random allele
                freq = orig_allele_freqs(allele_row,i); % get frequency of random allele
            end; % while

            %allele_state = allele_states(allele_row,i);
            
            for row=1:num_individs
            	% homozygote for allele, set both alleles in genotype to NaN
                if genotypes(row,(2*i)-1) == allele_states(allele_row,i) && genotypes(row,2*i) == allele_states(allele_row,i)
                    mod_genotypes(row,(2*i)-1) = NaN;
                    mod_genotypes(row,2*i) = NaN;
                end

                % heterozygote for allele, make into homozygote for second allele
                if genotypes(row,(2*i)-1) == allele_states(allele_row,i) && genotypes(row,2*i) ~= allele_states(allele_row,i)
                    mod_genotypes(row,(2*i)-1) = genotypes(row,2*i);
                end
                
                % heterozygote for allele, make into homozygote for first allele
                if genotypes(row,(2*i)-1) ~= allele_states(allele_row,i) && genotypes(row,2*i) == allele_states(allele_row,i)
                    mod_genotypes(row,2*i) = genotypes(row,(2*i)-1);
                end
                
            end % for row=1:num_individs
            
        end; % if total_poss_nulls(1,i) > 0
        
    end; % for loop
    
    
    