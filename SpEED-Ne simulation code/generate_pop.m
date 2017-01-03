function [output] = generate_pop(nloci,k,Ne,gens,write_file,outfilename)

% Function to generate a sample of genotypes from a population with known effective population size.
% All loci have the same number of alleles, alleles have equal initial frequencies, and that allelic
% states are adjacent integers.
%
%   version 1.2  02 January 2017
%
% Inputs:
%   nloci - number of independent loci
%   k - number of alleles at each locus
%   Ne - number of individuals sampled at final time point
%   gens - number of generations to simulate before sample taken
%   write_file - boolean that dictates whether or not to output a file
%   outfilename - name of output file if one is written
%
% Outputs
%   genotypes - n_individs by 2*nloci matrix of allelic states 
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
%   A copy of the GNU General Public License is available at http://www.gnu.org/licenses/.
%**************


    rng('shuffle'); %initialize random number generator
    
    states = [1:k]; % define allelic states
    n_individs = 100000; % number of individuals in founding population

    init_freqs = zeros(1,k);

    for i=1:k
        init_freqs(1,i) = 1/k;
    end

    founders = zeros(n_individs, 2*nloci, 1); % genotypes of large founding population
    
    genotypes = zeros(Ne, 2*nloci, gens); % genotypes of population of size Ne individuals
    
    
    % populate 1st generation by sampling alleles using initial freqs as expected values
    for i=1:nloci
        for j=1:n_individs
            
            founders(j,2*i-1,1) = drawstate(init_freqs,k,states);
            founders(j,2*i,1) = drawstate(init_freqs,k,states);

        end %for j
    end %for i

    
    % sample Ne individuals from founder population under randon mating 
    % to establish 2nd generation
    for i=1:Ne

        % sample parents (with replacement so selfing permitted) for
        % each progeny in the next generation
        parent_one = randi(n_individs);
        parent_two = randi(n_individs);

        for locus=1:nloci

            % sample allele from parent one
            if rand < 0.5 %determine which allele to sample from parent
                genotypes(i,2*locus-1,1) = founders(parent_one,2*locus-1);
            else
                genotypes(i,2*locus-1,1) = founders(parent_one,2*locus-1);
            end

            % sample allele from parent two
            if rand < 0.5 %determine which allele to sample from parent
                genotypes(i,2*locus,1) = founders(parent_two,2*locus-1);
            else
                genotypes(i,2*locus,1) = founders(parent_two,2*locus);
            end

        end %for locus
    end %for i

    
    % have population reproduce by random mating
    for g = 2:gens
        
        % sample one locus gametes for each parent and combine to make
        % progeny genotype at each locus
        for i=1:Ne
            
            % sample parents (with replacement so selfing permitted) for
            % each progeny in the next generation
            parent_one = randi(Ne);
            parent_two = randi(Ne);

            for locus=1:nloci

                % sample allele from parent one
                if rand < 0.5 %determine which allele to sample from parent
                    genotypes(i,2*locus-1,g) = genotypes(parent_one,2*locus-1,g-1);
                else
                    genotypes(i,2*locus-1,g) = genotypes(parent_one,2*locus,g-1);
                end

                % sample allele from parent two
                if rand < 0.5 %determine which allele to sample from parent
                    genotypes(i,2*locus,g) = genotypes(parent_two,2*locus-1,g-1);
                else
                    genotypes(i,2*locus,g) = genotypes(parent_two,2*locus,g-1);
                end

            end %for locus
        end %for i
                
    end %for gens
    
    %genotypes
    
    % sample individuals from base population and export to Excel file
    if Ne == n_individs
        output = genotypes(:,:,gens);
        
        if write_file
            xlswrite(outfilename,output);
        end % write_file
    else
        % sample first Ne individuals from population of n_individs
        output = genotypes(1:Ne,:,gens);
        
        if write_file
            xlswrite(outfilename,output);
        end % write_file
    end
    
    
    function [allele] = drawstate(freqs,k,states)
        % function to randomly sample an allele given an allele freq
        % distribution and a vector of allelic states
        
        num = rand; % draw random number
        
        cumulative_allele_freq = freqs(1,1);
        for i=1:k
            
            if num < cumulative_allele_freq
                allele = states(1,i);
                break; %exit for loop since we have our allele
            else
                cumulative_allele_freq = cumulative_allele_freq + freqs(1,i);
            end
            
        end
        
