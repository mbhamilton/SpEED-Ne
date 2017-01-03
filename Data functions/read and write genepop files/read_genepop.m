function [genotypes,num_pops,total_ind,ind_per_pop,sample_names] = read_genepop(text_file_name)
% Function to obtain genotype data from text file containing GenePop format genotype data. 

% version 1 - 03 Aug 2015

% Inputs:
% 
% text_file_name: the file name or path name of a text input file
% 
% input: an text file with:
% 1st line: string describing data set or with comments 
% 2nd line: the name of the first locus (locus names all must contain at least one character)
% 3rd: the name of the second locus (if needed)
%
% Alternate form are all locus names on one line with each name separated
% by a comma (i.e. locus_1, locus_2, locus_3)
%
% Line number of loci + 1: "POP" or "Pop", or "pop" to indicate the first population
% Line number of loci + 2: genotype data in two or three digit per allele format 
%
%
% Outputs:
% 
% genotypes: a total number of individuals by 2xnumber of loci matrix of integer allelic state data; pairs of columns define the two alleles in a genotype
% num_pops: integer number of populations
% total_ind: total number of individuals with genotype data summed over all populations
% ind_per_pop: num_pops by one vector of number of individuals with genotype data for each population
% sample_names: total number of individuals by one vector of sample names
% 
%**************
%   Copyright 2017 Matthew B Hamilton.
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


% inFile = 'example genepop file two.txt'; %name of the input text file containing Genepop format genotype data
[rawdata] = importdata(text_file_name, '\n'); %open the file and read in data, generates cell array of file contents

% determine how many lines until first occurance of POP or Pop or pop
line = 2;
found_pop = false;

while found_pop == false
    string = char(rawdata(line,1)); %convert line to character array
    line = line + 1; %increment counter
    
    result = strfind(string, 'POP');
    if result == 1
        result1 = true;
    else
        result1 = false;
    end
        
    result = strfind(string, 'Pop');
    if result == 1
        result2 = true;
    else
        result2 = false;
    end
    
    result = strfind(string, 'pop');
    if result == 1
        result3 = true;
    else
        result3 = false;
    end
    
    if result1 == true || result2 == true || result3 == true
        found_pop = true;
        line = line - 1; %decrement counter since we found it
    end
    
end % while

% determine how many loci
if (line - 2) == 1 % locus name or names takes up only one row
    % therefore, there must be only one locus or the locus names are comma separated
    lines_to_genotypes = 3;

    string = char(rawdata(2,1)); % convert 2nd line to character array
    raw_loc_names = strsplit(string, ','); % split string into separate names separated by commas
    num_loci = numel(raw_loc_names);
    %loc_names = cell(num_loci,1); % allocate space
    
    % remove white spaces from locus names
    for i=1:num_loci
        trimmed_string = strtrim(raw_loc_names(1,i)); 
        loc_names(i,1) = trimmed_string; % store locus name
    end
    
    lines_to_genotypes = 3;

else % locus names are on multiple lines
    
    num_loci = line - 2; % subtract one for line with POP and first line
    loc_names = cell(num_loci,1); % allocate space in cell array
    
    for i=1:num_loci
        loc_names(i,1) = rawdata(i+1,1);
    end
    
    lines_to_genotypes = num_loci + 2;
end

%num_loci
%loc_names

% determine how many pops and how many individuals within each pop
total_lines = size(rawdata, 1); % total lines in file
num_pops = 0; % initialize counter
counter = 0; % initialize counter

for i = lines_to_genotypes:total_lines
    string = char(rawdata(i,1)); %convert line to character array
    counter = counter + 1; % counter of individuals within pops
    
    result = strfind(string, 'POP'); % strfind returns [] if string not found
    if isempty(result)
        result1 = false;
    else
        result1 = true;
    end
        
    result = strfind(string, 'Pop'); % strfind returns [] if string not found
    if isempty(result)
        result2 = false;
    else
        result2 = true;
    end
    
    result = strfind(string, 'pop'); % strfind returns [] if string not found
    if isempty(result)
        result3 = false;
    else
        result3 = true;
    end
    
    if result1 == true || result2 == true || result3 == true
        found_pop = true;
        num_pops = num_pops + 1; % increment counter
        pop_start_rows(num_pops,1) = i+1; % store line number of line after where POP found
        counter = 0; % reset counter of individuals within pops
    else
        ind_per_pop(num_pops,1) = counter; % increment counter of individuals within pop
    end
    
end % for i

total_ind = sum(ind_per_pop); % total individuals with genotype data in file


% determine the number of digits used to code alleles
string = char(rawdata(lines_to_genotypes+2,1)); % read first line of genotype data, convert line to character array
raw_line_data = strsplit(string, ','); % split string into sections separated by comma

genotype_string = char(raw_line_data(1,2)); % convert genotype cell contents to array of char
genotype_trim_string = strtrim(genotype_string); % remove leading and trailing spaces
genotype_data = strsplit(genotype_trim_string); % part of split string with genotypes into sections separated by white space
genotype = genotype_data{1,1}; % use curly brackets toaccess the content of the cell
allele_length = size(genotype,2)/2; % half the number of characters in first genotype


% parse lines into sample names and numerical allelic states
sample_names = cell(total_ind,1); % allocate space
genotypes = zeros(total_ind, 2*num_loci); % allocate space

counter = 0; % initialize
for pop = 1:num_pops
    for line = 1:ind_per_pop(pop,1)

        counter = counter + 1; % lines in genotypes matrix
        index = (pop_start_rows(pop,1) - 1) + line;

        string = char(rawdata(index,1)); % read line, convert line to character array
        raw_line_data = strsplit(string, ','); % split string into sections separated by comma

        sample_names(counter,1) = strtrim(raw_line_data(1,1)); % save sample name

        genotype_string = char(raw_line_data(1,2)); % convert cell contents to array of char
        genotype_trim_string = strtrim(genotype_string); % remove leading and trailing spaces
        genotype_data = strsplit(genotype_trim_string); % part of split string with genotypes into sections separated by white space

        for g = 1:num_loci
            col = (2*g)-1;

            if allele_length == 2
                alleles = sscanf(genotype_data{1,g}, '%2d%2d');
            elseif allele_length == 3
                alleles = sscanf(genotype_data{1,g}, '%3d%3d');
            else
                fprintf('allele length not 2 or 3');
            end

            genotypes(counter,col) = alleles(1,1); 
            genotypes(counter,col+1) = alleles(2,1);

        end % for g

        end % for line = 1:ind_per_pop(pop,1)
end % for pop = 1:num_pops

