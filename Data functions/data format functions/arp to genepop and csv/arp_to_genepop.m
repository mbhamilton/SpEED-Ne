function arp_to_genepop
% Function to read text file containing arlequin format genotype data and write that genotype data 
% out in a genepop format file and a comma separated value (.csv) format file. 
% The program assumes a single population.
%
% version 1 - 08 Nov. 2016

% Inputs: 
%   none - inputs coded into function at top of file
%
% Outputs:
%   function writes new file named outFile in genepop format
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


    arp_file_name = 'fsc_msat_u0p1_r0p5_Ne100_1_1.arp';    %name or path of the input arlequin format data file

    gp_file_name = 'fsc_msat_u0p1_r0p5_Ne100_1_1.arp_genepop.txt';   %name of tab delimited text file for output in genepop format

    csv_file_name = 'fsc_msat_u0p1_r0p5_Ne100_1_1.arp.csv';    %name or path of the output .csv format data file

    
    % read in .arp file
    [genotypes,num_individuals,num_loci,result] = read_arp(arp_file_name);

    if result == 0
        fprintf('*** Genotype .arp file not converted currectly in call to function read_arp. ***\n');
    end;


    % write genepop format file, assuming a single population

    fid = fopen(gp_file_name, 'w+'); % open new file for writing

    fprintf(fid,'Genepop format conversion of arlequin file %s\n',arp_file_name); % comment on first line

    % write arbitary locus names, one per line
    for i=1:num_loci
        fprintf(fid,'locus_%d\n',i); % locus name
    end; % for i=1:num_loci
    
    fprintf(fid,'POP\n');

    for indiv=1:num_individuals

        fprintf(fid,'individual_%d, ',indiv);

        for locus=1:num_loci

            fprintf(fid,'%03.0f%03.0f ', genotypes(indiv, (locus*2)-1), genotypes(indiv, locus*2));

        end %for loci

        fprintf(fid,'\n');

    end %for indiv

    % write Excel format file
    dlmwrite(csv_file_name,genotypes);
