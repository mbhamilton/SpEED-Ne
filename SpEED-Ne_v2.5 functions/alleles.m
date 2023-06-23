function [num_alleles, allele_freqs, allele_states] = alleles(genotypes,num_individuals,num_loci);
% Function determine numbers of alleles per locus, allele frequencies at each locus,
% and allelic states.
%
%   Inputs:
%   genotypes - the genetype data with pairs of columns containing the two
%               integer alleles for a diploid locus, rows represent individuals
%   num_individuals - integer number of rows in the genotype data matrix
%   num_loci - integer half the number of columns in the genotype data matrix
%   print_table - boolean switch for output of table of pairwise estimates
%
%   Outputs:
%   num_alleles - one by num_loci vector of the integer number of alleles per locus
%   allele_freqs - max num_alleles by num_loci vector of allele frequencies
%                  (entries without an allele are filled with zeros)
%   allele_states - max num_alleles by num_loci vector of integer allelic
%                   states (entries without an allele are filled with zeros)
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


    % organize data into alleles 
    alleles = zeros(2*num_individuals,num_loci); %allocate space for column vectors of all alleles at each locus

    cols = 2*num_loci;
    
    j = 1;
    for i=1:2:cols
        col_one = genotypes(:,i);
        col_two = genotypes(:,i+1);

        alleles(:,j) = cat(1,col_one, col_two);

        j = j+1;

    end

    % determine alleles present at each locus, get number of alleles per locus
    for i=1:num_loci
        table = tabulate(alleles(:,i)); % get allele states (col 1), counts (col 2), and freqs (row 3) of alleles
        indices = table(:,2) > 0;   % get indices of rows that have a count > zero
        num_alleles(1,i) = sum(indices); %count up indices that are one to get number of alleles
    end

    max_alleles = max(num_alleles); %determine the maximum number of alleles

    allele_states = zeros(max_alleles, num_loci);   %allocate space
    allele_counts = zeros(max_alleles, num_loci);   %allocate space
    allele_freqs = zeros(max_alleles, num_loci);   %allocate space

    % record allelic states and allele frequencies for each state
    for i=1:num_loci
        table = tabulate(alleles(:,i)); % get allele states (col 1), counts (col 2), and freqs (row 3) of alleles
        rows_in_table = size(table,1);   % get number of rows in table produced by tabulate function

        state = 1;  %initialize counter
        for j=1:rows_in_table
            if table(j,2) > 0
                allele_states(state, i) = table(j,1);    %record allele state
                allele_counts(state, i) = table(j,2);    %record allele count
                allele_freqs(state, i) = table(j,3)/100;    %record allele frequency, must divide by 100 since tabulate provides percent
                state = state + 1;  %increment counter
            end

        end%for j
    end    