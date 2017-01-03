function [genotypes,num_individuals,num_loci,result] = read_arp(text_file_name)
% Function to obtain two column genotype data from text file containing Arlequin format genotype data. 
%
% version 1.1 - 02 January 2017
%
% Inputs:
%   text_file_name: the file name (if in same directory as code) or full path name of a text input file with .arp extension
%
% Outputs:
%	genotypes: number of individuals by 2xnumber of loci matrix of integer allelic state data
%	num_individuals: number of individuals and rows of genotype data
%	num_loci: integer number of loci in genotype data
%	result: 1 if end of genotype data '}' found and zero if not found
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


% inFile = 'example genepop file two.txt'; %name of the input text file containing Genepop format genotype data
[rawdata] = importdata(text_file_name, '\n'); %open the file and read in data, generates cell array of file contents

% determine how many lines until first occurance of POP or Pop or pop
line = 2;
found_data = false;

% read down until SampleData= { line
while found_data == false
    string = char(rawdata(line,1)); %convert line to character array
    line = line + 1; %increment counter
    
    result = strfind(string, '		SampleData= {');
    if result == 1
        found_data = true;
    else
        found_data = false;
    end
        
end % while


% get sample size of individuals give by SampleSize=xx on line
% before SampleData= { 
string = char(rawdata(line - 2,1)); %convert line to character array
trim_string = strtrim(string); % remove leading and trailing spaces
cell_values = strsplit(trim_string,'='); % split string into separate values separated by equals sign
num_individuals = str2double(cell_values{1,2}); % convert cell value to double

% get number of loci by determining number of alleles on second line of
% data after SampleData= {.
% Use second line of data since first line has individual name.
string2 = char(rawdata(line+1,1)); % convert line to cell array
trim_string2 = strtrim(string2); % remove leading and trailing spaces
raw_alleles2 = strsplit(trim_string2); % split string into separate values separated by any whitespace

num_loci = numel(raw_alleles2);

genotypes = zeros(num_individuals,2*num_loci); % allocate space
counter = 1; % initialize for genotype rows
% read in genotype data and write to num_individuals by num_loci genotype matrix
for row=line:2:((line-1) + 2*num_individuals)

    string1 = char(rawdata(row,1)); % convert line to cell array
    trim_string1 = strtrim(string1); % remove leading and trailing spaces
    raw_alleles1 = strsplit(trim_string1); % split string into separate values separated by any whitespace

    string2 = char(rawdata(row+1,1)); % convert line to cell array
    trim_string2 = strtrim(string2); % remove leading and trailing spaces
    raw_alleles2 = strsplit(trim_string2); % split string into separate values separated by any whitespace

    for col=1:num_loci
    
        genotypes(counter,2*col-1) = str2double(raw_alleles1{1,col+2}); % convert cell value to double, col+2 to skip name fields
        genotypes(counter,2*col) = str2double(raw_alleles2{1,col}); % convert cell value to double
    
    end; % for col=1:num_loci

    counter = counter + 1; % increment

end; % for row


% test for } at end of genotype data
string = char(rawdata(row+2,1)); %convert line to character array

result = strfind(string, '}'); % returns 1 if found and zero if not found


