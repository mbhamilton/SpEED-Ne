function excel_to_genepop
% Code to take genotype data in Excel spread sheet and create a text file in
% GenePop format. 

% version 4 - 28 July 2016

% input: an Excel file with:
% 1st row: number of loci; comment
% 2nd row: number of populations; comment
% 3rd row: title of dataset
% 4th row: locus names
% 5th row: name of 1st population
% numeric data starts in 6th row, 5 rows down
% each population is separated by one line with the population name in the left-most cell
% Population names must contain at least one character so that they are NaN.
% Each population must contain at least one individual with non-missing data at each locus.
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



inFile = 'Msax 2012 Ches Hud Del 12 loci v1 genotypes only.xlsx';    %name of the input Excel file
outFile = 'Msax 2012 Ches Hud Del 12 loci v1 genepop.txt';   %name of tab delimited text file for output

[data, labels] = xlsread(inFile);    %open the file and read in data, get numeric data and text data separately

%parse input file for genotypes and alphanumeric data
[nloci,npops,genotypes,popsizes,popindices,total_N,cols] = getdata(data,labels);

fid = fopen(outFile, 'w+');

counter = 1;
for pop=1:npops
    
    fprintf(fid,'POP\n');
    fprintf('POP\n')
    
    if pop == 1
        pop_name_row = 3;
    else
        pop_name_row = pop_name_row + popsizes(pop-1,1) + 1;
    end

    
    for indiv=1:popsizes(pop,1)
    
        pop_name = char( labels(pop_name_row,1) );
        fprintf(fid,'%s, ', pop_name);
        fprintf('%s, ', pop_name)
        
        for locus=1:nloci
            
            fprintf(fid,'%03.0f%03.0f ', genotypes(counter, (locus*2)-1), genotypes(counter, locus*2));

            fprintf('%03.0f%03.0f ', genotypes(counter, (locus*2)-1), genotypes(counter, locus*2))
        end %for loci
        
        fprintf(fid,'\n');
        fprintf('\n')
        
        counter = counter +1;
        
    end %for indiv
    
end %for pops


%_______________________________________
function[nloci,npops,genotypes,popsizes,popindices,rows,cols] = getdata(data,labels)    %parse file for genotypes and other data
% input: an Excel file with:
% 1st row: number of loci; comment
% 2nd row: number of populations; comment
% 3rd row: title of dataset
% 4th row: locus names
% 5th row: name of 1st population
% numeric data starts in 6th row, 5 rows down
% each population is separated by one line with the population name in the left-most cell
% Population names must contain at least one character so that they are NaN.
% Each population must contain at least one individual with non-missing data at each locus.

nloci = data(1,1);   %number of loci, top left cell
npops = data(2,1);   %number of populations, one cell below top left cell

%calculate individual population sizes
start = 6;                  %numerical values always start on 6th row
popsizes = zeros(npops,1);  %allocate space
popindices = zeros(npops,2);    %allocate space
vec = isnan(data(:,1));     %get vector of 0 = non-nan, 1 = nan, p. 160 of Hanselman
maxindex = size(vec,1);

vec(maxindex+1,1)=1;    %cram an extra 1 for nan into last cell since it is lacking
maxindex = maxindex+1;  %increment maxindex by one

counter=0;
pop=1;
for i=6:maxindex
    % count rows until encountering a value of 1 (equal to nan).  
    % NANs are population name fields that separate population data blocks.
    if (vec(i,1) == 1)
        popsizes(pop,1) = counter;
        popindices(pop,2)=i-1;  %last row of data for this pop
        popindices(pop,1)=i-counter;  %first row of data for this pop
        counter=0;
        pop=pop+1;
    else
        counter=counter+1;
    end
end
% popsizes
% popindices

% now need to get JUST numerical data into a matrix
rows = sum(popsizes,1); %sum down column
cols = 2*nloci;  %total number of allele entries, 2 per locus
genotypes = zeros(rows,cols);   %allocate space

% counter=1;
% for i=1:npops
%     for j=popindices(i,1):1:popindices(i,2)
%         genotypes(counter,:) = data(j,:);
%         counter=counter+1
%     end
% end
% %genotypes

for i=1:npops
    
    start_row = popindices(i,1);
    stop_row = popindices(i,2);
    
    geno_start = (start_row - 5) - (i-1);
    geno_stop = (stop_row - 5) - (i-1);
    
    genotypes(geno_start:geno_stop,:) = data(start_row:stop_row,:);
    
end



