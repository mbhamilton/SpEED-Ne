function [output] = generate_pop_r(mating_model,r,nloci,k,Ne,gens,write_file,outfilename)

% Function to generate a sample of genotypes from a population with known effective population size.
% All loci have the same number of alleles, alleles have equal initial frequencies, and that allelic
% states are adjacent integers.
%
% This version of the function allows for the user to input an array of
% recombination rates between each pair of loci, for the number of loci specified.
% Recombination is simulated by checking if a randomly generated number
% is higher or lower than the given recombination rate, with the latter
% producing recombination. one_is_top and two_is_top keep track of whether
% the allele from the first and second parent chosen was the "top" option,
% as opposed to the other "bottom" option. If the loci experience
% recombination, the opposite allele position is chosen at the next locus, and if
% there is no recombination the same allele position is chosen. 
%
%   version 1.3  03 January 2017
%
% Inputs:
%   mating_model - integer to specify random mating with possibility of selfing (=1)
%       or random mating without the possibility of selfing (=2)
%   r - an array containing nloci-1 recombination rate values between each
%   pair of loci
%   nloci - integer number of independent loci
%   k - integer number of alleles at each locus
%   Ne - integer number of individuals sampled at final time point
%   gens - integer number of generations to simulate before sample taken
%   write_file - boolean that dictates whether or not to output a file
%   outfilename - strong name of output file if one is written
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
% to establish 2nd generation, ignore possibility of selfing since
% founder population is very large
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

if mating_model == 1
    
    % random mating with the possibility of self-fertilization
    for g = 2:gens
        
        % sample one locus gametes for each parent and combine to make
        % progeny genotype at each locus
        for i=1:Ne
            
            % sample parents (with replacement so selfing permitted) for
            % each progeny in the next generation
            parent_one = randi(Ne);
            parent_two = randi(Ne);
    
            % boolean values that keep track of whether the allele chosen
            % from each parent is the "top" or "bottom" allele
            one_is_top = true;
            two_is_top = true;
            
            % sample initial allele from parent one randomly
            if rand < 0.5
                genotypes(i,1,g) = genotypes(parent_one,1,g-1);
            else
                genotypes(i,1,g) = genotypes(parent_one,2,g-1);
                one_is_top = false; % means bottom allele was chosen
            end
            % sample initial allele from parent two randomly
            if rand < 0.5
                genotypes(i,2,g) = genotypes(parent_two,1,g-1);
            else
                genotypes(i,2,g) = genotypes(parent_two,2,g-1);
                two_is_top = false; % means bottom allele was chosen
            end
            
            
            % sample the rest of the loci based on r values provided by user
            for locus=2:nloci
                
                % sample allele from parent one
                if rand < r(locus-1) % loci undergo recombination
                    if one_is_top == true % if top was previously chosen, bottom is chosen in next loci
                        genotypes(i,2*locus-1,g) = genotypes(parent_one,2*locus,g-1);
                        one_is_top = false; %change allele being tracked
                    else % bottom was prevously chosen
                        genotypes(i,2*locus-1,g) = genotypes(parent_one,2*locus-1,g-1); % top is chosen 
                        one_is_top = true; % change allele being tracked
                    end
                else % loci do not undergo recombination
                    if one_is_top == true % if top was previously chosen, top is chosen again 
                        genotypes(i,2*locus-1,g) = genotypes(parent_one,2*locus-1,g-1);
                    else   % if bottom was previously chosen, bottom is chosen again
                        genotypes(i,2*locus-1,g) = genotypes(parent_one,2*locus,g-1);
                    end
                end
                
                
                % sample allele from parent two, same process as for parent
                % one
                if rand < r(locus-1)
                    if two_is_top == true
                        genotypes(i,2*locus,g) = genotypes(parent_two,2*locus,g-1);
                        two_is_top = false;
                    else
                        genotypes(i,2*locus,g) = genotypes(parent_two,2*locus-1,g-1);
                        two_is_top = true;
                    end
                    
                else
                    if two_is_top == true
                        genotypes(i,2*locus,g) = genotypes(parent_two,2*locus-1,g-1);
                    else
                        genotypes(i,2*locus,g) = genotypes(parent_two,2*locus,g-1);
                    end
                    
                end
                
                
            end %for locus
            
        end %for i
        
    end %for gens
    
end; % if mating_model == 1


if mating_model == 2
    % random mating without the possibility of self-fertilization
    for g = 2:gens
        
        % sample one locus gametes for each parent and combine to make
        % progeny genotype at each locus
        for i=1:Ne
            
            % sample parents (with replacement so selfing permitted) for
            % each progeny in the next generation
            parent_one = randi(Ne);
            
            parent_two = randi(Ne);
            while parent_one == parent_two
                parent_two = randi(Ne); % keep drawing random parent two
            end; % while
            
            % boolean values that keep track of whether the allele chosen
            % from each parent is the "top" or "bottom" allele
            one_is_top = true;
            two_is_top = true;
            
             % sample initial allele from parent one randomly
            if rand < 0.5
                genotypes(i,1,g) = genotypes(parent_one,1,g-1);
            else
                genotypes(i,1,g) = genotypes(parent_one,2,g-1);
                one_is_top = false; % means bottom allele was chosen
            end
            % sample initial allele from parent two randomly
            if rand < 0.5
                genotypes(i,2,g) = genotypes(parent_two,1,g-1);
            else
                genotypes(i,2,g) = genotypes(parent_two,2,g-1);
                two_is_top = false; % means bottom allele was chosen
            end
            
             % sample the rest of the loci based on r values provided by user
            for locus=2:nloci
                
                % sample allele from parent one
                if rand < r(locus-1) % loci undergo recombination
                    if one_is_top == true % if top was previously chosen, bottom is chosen in next loci
                        genotypes(i,2*locus-1,g) = genotypes(parent_one,2*locus,g-1);
                        one_is_top = false; %change allele being tracked
                    else % bottom was prevously chosen
                        genotypes(i,2*locus-1,g) = genotypes(parent_one,2*locus-1,g-1); % top is chosen 
                        one_is_top = true; % change allele being tracked
                    end
                else % loci do not undergo recombination
                    if one_is_top == true % if top was previously chosen, top is chosen again 
                        genotypes(i,2*locus-1,g) = genotypes(parent_one,2*locus-1,g-1);
                    else   % if bottom was previously chosen, bottom is chosen again
                        genotypes(i,2*locus-1,g) = genotypes(parent_one,2*locus,g-1);
                    end
                end
                
                
                % sample allele from parent two, same process as for parent
                % one
                if rand < r(locus-1)
                    if two_is_top == true
                        genotypes(i,2*locus,g) = genotypes(parent_two,2*locus,g-1);
                        two_is_top = false;
                    else
                        genotypes(i,2*locus,g) = genotypes(parent_two,2*locus-1,g-1);
                        two_is_top = true;
                    end
                    
                else
                    if two_is_top == true
                        genotypes(i,2*locus,g) = genotypes(parent_two,2*locus-1,g-1);
                    else
                        genotypes(i,2*locus,g) = genotypes(parent_two,2*locus,g-1);
                    end
                    
                end
                
            end %for locus
        end %for i
        
    end %for gens
    
end; % if mating_model == 2


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

