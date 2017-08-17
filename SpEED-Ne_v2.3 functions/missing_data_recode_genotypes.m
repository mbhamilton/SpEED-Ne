function [genotypes] = missing_data_recode_genotypes(genotypes, num_loci, num_individuals,missing_data_value);
% Scan genotype data for missing values (integer <= missing_data_value) and create version of data where missing
% alleles are re-coded to NaN AND the entire genotype is recoded to NaN. (NaN will be recognized as missing data in other functions). 
%
%   Inputs:
%   genotypes - the genotype data with pairs of columns containing the two
%               integer alleles for a diploid locus, rows represent individuals
%   num_loci - integer half the number of columns in the genotype data matrix
%   num_individuals - integer number of rows in the genotype data matrix
%   missing_data_value - integer that defines upper limit of missing data
%                       values
%
%   Outputs:
%   genotypes - genotype data with all genotypes with missing data entries recoded as NaN

    cols = 2*num_loci;
    
    % scan genotypes for one allele with missing data
    for i=1:2:cols
        for j=1:num_individuals

            if (genotypes(j,i) <= missing_data_value) || (genotypes(j,i+1) <= missing_data_value) || (isnan(genotypes(j,i))) || (isnan(genotypes(j,i+1)))
                genotypes(j,i) = NaN;   %set first allele of genotype to NaN
                genotypes(j,i+1) = NaN;   %set second allele of genotype to NaN
            end

        end
    end
