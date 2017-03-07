function [table, S, paa, pbb] = comp_hap_table(G1,G2,k1,k2,states1,states2)
% Function to compute composite haplotype table from genotypes observed at
% two loci. Also tally homozygotes observed at each locus.
%
% This version tallies haplotypes rather than genotypes.
%
%   Inputs:
%   G1 - num_individuals by 2 matrix of integer genotype data for locus one
%   G2 - num_individuals by 2 matrix of integer genotype data for locus two
%       -> funciton assumes that genotype data have already been checked
%       -> for missing data and that only full rows are provided
%
%   k1 - integer number of alleles for locus one
%   k2 - integer number of alleles for locus two
%   states1 - k1 by 1 vector of allelic states for locus one
%   states2 - k2 by 1 vector of allelic states for locus two
%
%   Outputs:
%   table - k1 by k2 matrix of integers that make the composite haplotype table 
%   S - integer sample size of individuals
%   paa - k1 by 1 vector of the observed frequency of homozygous genotypes at locus one
%   pbb - k2 by 1 vector of the observed frequency of homozygous genotypes at locus two
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

    
	if(size(G1,1)~= size(G2,1)),
		error('unequal lengths of G1 and G2 in function comp_hap_table');
	end;

    tally = 0; %initialize
    count = 0;
    
    paa = zeros(k1,1); %initialize
    pbb = zeros(k2,1); %initialize
    
    table = zeros(k1, k2); %initialize
    
    index = zeros(size(G1,1),4); %initialize
    
    for loc_one_allele = 1:k1
        
        % count up locus one homozygote genotypes
        index(:,1) = (states1(loc_one_allele,1) == G1(:,1)) & (states1(loc_one_allele,1) == G1(:,2));
        paa(loc_one_allele,1) = sum(index(:,1));
        
        for loc_two_allele = 1:k2

            count = 0; %initialize
            index(:,:) = 0; %initialize

            %count x_ by x_ haplotype
            index(:,1) = (states1(loc_one_allele,1) == G1(:,1)) & (states2(loc_two_allele,1) == G2(:,1));
            count = sum(index(:,1));
            table(loc_one_allele, loc_two_allele) = table(loc_one_allele, loc_two_allele) + count;
            tally = tally + count;
                        
            %count _x by x_ haplotype
            index(:,2) = (states1(loc_one_allele,1) == G1(:,2)) & (states2(loc_two_allele,1) == G2(:,1));
            count = sum(index(:,2));
            table(loc_one_allele, loc_two_allele) = table(loc_one_allele, loc_two_allele) + count;
            tally = tally + count;

            %count x_ by _x haplotype
            index(:,3) = (states1(loc_one_allele,1) == G1(:,1)) & (states2(loc_two_allele,1) == G2(:,2));
            count = sum(index(:,3));
            table(loc_one_allele, loc_two_allele) = table(loc_one_allele, loc_two_allele) + count;
            tally = tally + count;

             %count _x by _x haplotype
            index(:,4) = (states1(loc_one_allele,1) == G1(:,2)) & (states2(loc_two_allele,1) == G2(:,2));
            count = sum(index(:,4));
            table(loc_one_allele, loc_two_allele) = table(loc_one_allele, loc_two_allele) + count;
            tally = tally + count;

        end %for loc_two_allele = 1:k2        
	end %loc_one_allele = 1:k1
    
    
    for loc_two_allele = 1:k2
        % count up locus two homozygote genotypes
        % do it outside of nested loop above to be more efficient
        index(:,:) = 0;

        index(:,1) = (states2(loc_two_allele,1) == G2(:,1)) & (states2(loc_two_allele,1) == G2(:,2));
        pbb(loc_two_allele,1) = sum(index(:,1));
    end
    
    S = tally/4;
    
    paa = paa./S;   %divide counts by sample size to make frequency
    pbb = pbb./S;   %divide counts by sample size to make frequency
