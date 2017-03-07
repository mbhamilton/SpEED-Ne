function allele_freq_table(num_loci, num_alleles, allele_freqs, allele_states)
% Function to print out a table of allelic states and frequencies for all
% loci.
%
%   Inputs:
%   num_loci - integer number of loci
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

    % make allele state and frequency output table
    fprintf('================================================\n');
    fprintf('Table of allele states and frequencies\n\n');
    fprintf('Locus\tk\tallele\tfrequency\n');
    fprintf('------------------------------------------------\n');

    for locus=1:num_loci

        fprintf('%3d\t',locus); % locus
        fprintf('%3d\n',num_alleles(1,locus)); % num alleles for this locus

        for state=1:num_alleles(1,locus)
            fprintf('\t\t%3d\t',allele_states(state,locus)); % allelic state
            fprintf('%4.4f\n',allele_freqs(state,locus)); % frequency of state
        end % for state

    end %for row
    fprintf('================================================\n\n');
