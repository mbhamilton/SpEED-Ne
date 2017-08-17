function [miss_allele_genotypes,miss_geno_genotypes] = sim_missing_data(genotypes,num_individs,num_loci,freq_loci_missing,freq_missing_individs)

% version 1.0 07 July 2017
%
% Inputs:
%   genotypes - number of individuals by 2*num_loci matrix of integers of allelic states
%   num_individs - number of rows in genotypes
%   num_loci - integer half the number of columns in the genotype data matrix
%   freq_loci_missing - frequency of loci with missing data
%   freq_missing_individs - expected frequency of individuals (rows) with
%       missing data per locus
%
% Outputs:
%   miss_allele_genotypes - number of individuals by 2*num_loci matrix of genotypes
%       with missing data for individual alleles
%   miss_geno_genotypes - number of individuals by 2*num_loci matrix of genotypes
%       with missing data for pairs of alleles (one locus genotypes)
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

    rng('shuffle'); %initialize random number generator with a seed based on the current time
    
    miss_allele_genotypes = genotypes; % initialize
    miss_geno_genotypes = genotypes; % initialize
    
    for i=1:num_loci
        if rand <= freq_loci_missing % this locus has missing data
            
            for row=1:num_individs
            
                % determine if row has missing data
                if rand <= freq_missing_individs
                    if rand <= 0.5
                        miss_allele_genotypes(row,(2*i)-1) = NaN;   % missing allele in first column of locus data
                    else
                        miss_allele_genotypes(row,2*i) = NaN;   % missing allele in second column of locus data
                    end
                
                    %set both alleles to missing in miss_geno_genotypes version of data
                    miss_geno_genotypes(row,(2*i)-1) = NaN;
                    miss_geno_genotypes(row,2*i) = NaN;
                end
                
            end % end for row
        end % if
    end %for i
    