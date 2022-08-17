function varargout = SpEED_Ne(varargin)
%   SpEED-Ne: A program to estimate gametic (linkage) disequilibrium and the genetic 
%   effective population size (Ne) using a multilocus diploid genotype dataset. 
%   The Statistics Toolbox of Matlab is required. Uses SpEED_Ne.fig. 
%  
%   version 2.3  17 Aug 2022
%
%**************
%   Copyright 2022 Matthew B Hamilton.
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
%
%	Last Modified by GUIDE v2.5 18-Aug-2016 07:42:15


dbstop if error % debug if program crashes


% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SpEED_Ne_OpeningFcn, ...
                   'gui_OutputFcn',  @SpEED_Ne_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before SpEED_Ne is made visible.
function SpEED_Ne_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SpEED_Ne (see VARARGIN)

% Choose default command line output for SpEED_Ne
handles.output = hObject;

set(handles.missing_data_button_group,'SelectionChangeFcn',@missing_data_button_group_SelectionChangeFcn);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes SpEED_Ne wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = SpEED_Ne_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in run_button.
function run_button_Callback(hObject, eventdata, handles)

    % test for required toolbox
    v = ver;
    has_stats = any(strcmp({v.Name}, 'Statistics and Machine Learning Toolbox'));
    if has_stats == 0
        beep;
        fprintf('Sorry. Statistics and Machine Learning Toolbox required to run this app.\nExiting program.\n');
        return;  % no toolbox present, exit function
    end


    rng('shuffle'); %initialize random number generator with a seed based on the current time
    
    %get values from interface
    threshold_allele_freq = get(handles.min_allele_freq,'UserData');
    n_permutations = get(handles.data_permute_reps,'UserData');
    alpha = get(handles.alpha,'UserData');   % total width of confidence interval tails in jackknife and bootstrap
    in_file_path = get(handles.data_file_string,'String'); % full path name of input file
    missing_data = get(handles.missing_data_button_group,'UserData');
    
    % test check boxes for jackknifing
    jack_individuals = get(handles.jackknife_individs_checkbox,'Value'); % get state of checkbox, 1=checked, 0=unchecked
    jack_loci = get(handles.jackknife_loci_checkbox,'Value'); % get state of checkbox, 1=checked, 0=unchecked
    
    % test check box for permuting genotype data
    permute = get(handles.permute_checkbox,'Value'); % get state of checkbox, 1=checked, 0=unchecked
    
    % test check box for print table of allele states and frequencies
    print_allele_table = get(handles.allele_freq_table_checkbox,'Value'); % get state of checkbox, 1=checked, 0=unchecked
    
    % test check box for print table of r^2 estimates for each pair of alleles within all loci
    print_locus_table = get(handles.paired_loci_tables_checkbox,'Value'); % get state of checkbox, 1=checked, 0=unchecked
    
    % test check box for output of graphs
    make_graphs = get(handles.output_graphs_checkbox,'Value'); % get state of checkbox, 1=checked, 0=unchecked
    
    % test check box for output of files containing estimates and values
    make_files = get(handles.file_output_checkbox,'Value'); % get state of checkbox, 1=checked, 0=unchecked

    
    % test that input file has been specified
    if strcmpi(in_file_path,'no file specified') %compare strings, case insensitive
        beep;
        fprintf('No data file specified.\n');
        return;  % no file specified, exit function
    end
    
    genotypes = importdata(in_file_path); % read in gentype data from file        

    fprintf('----------------------------------------------------\n');
    fprintf('SpEED-Ne is free software under the terms of the GNU General Public License\n\n');
    
    t = datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss Z'); % get date and time
    date_string = datestr(t); % convert date and time to string
    fprintf('%s\n\n',date_string);
    
    fprintf('Input genotype data file: %s\n\n',in_file_path);

    [num_individuals, cols] = size(genotypes);   %get dimensions of genotype table
    num_loci = cols/2;
    
    fprintf('Number of loci in data set = %d \n',num_loci);
    fprintf('\n');

    % recode missing data
    missing_data_value = 0;
    
    if missing_data == 2
        % recode missing alleles as NaN
        [genotypes] = missing_data_recode_alleles(genotypes,num_loci,num_individuals,missing_data_value);
        fprintf('Single non-missing alleles used in estimate when other allele is missing.\n\n',threshold_allele_freq);
    else % missing data == 1 or default
        % recode entire genotypes with missing alleles as NaN
        [genotypes] = missing_data_recode_genotypes(genotypes,num_loci,num_individuals,missing_data_value); 
        fprintf('Entire genotypes treated as missing if one allele is missing.\n\n',threshold_allele_freq);
    end
    
    [orig_num_alleles, orig_allele_freqs, orig_allele_states] = alleles(genotypes,num_individuals,num_loci); % get allele frequencies and states
    
    if print_allele_table == 1
        fprintf('Allele frequencies in genotype data:\n');
        allele_freq_table(num_loci, orig_num_alleles, orig_allele_freqs, orig_allele_states); % print table of allele states and frequencies for all loci
    end
    
    % exclude alleles less than a given threshold frequency
    fprintf('Minimum allele frequency used for allele frequency threshold (AFT) estimates = %4.4f\n\n',threshold_allele_freq);

    %%%%%%%
    % make estimates from data
    print_notes = true;  % print error output strings to console
    
    if print_locus_table == 1
        print_table = true;  % print table of r^2 estimates for each pair of loci
    else
        print_table = false;  % do not print table of r^2 estimates for each pair of loci
    end


    fprintf('estimating for all pairs of loci ...\n');

    [r_squared_c_W,r_squared_c_TH,r_squared_delta_W,r_squared_delta_TH,r_squared_c_values,r_squared_delta_values,...
        S_c_values,S_delta_values,c_locus_pairs,delta_locus_pairs,all_locus_pairs_table,~] = ...
        make_estimate(print_notes,genotypes,num_individuals,num_loci,print_table,threshold_allele_freq,print_locus_table,make_files);

    
    % if checked, output estimates for all allele pairs and all locus pairs
    if make_files
        % prepare file names for output files
        [file_path,file_name,~] = fileparts(in_file_path); % get file name and path
        full_path = strcat(file_path,'/');
        path_plus_filename = strcat(full_path,file_name);
        pairwise_alleles_filename = strcat(path_plus_filename,'_all_allele_pair_estimates.csv');
        pairwise_loci_filename = strcat(path_plus_filename,'_all_locus_pairs_estimates.csv');
        
        cell_data = num2cell(all_locus_pairs_table);
        T = cell2table(cell_data,'VariableNames',{'locus1','locus2','S','k1','k2','r_sqr_c_weight','r_sqr_delta_weight','r_sqr_c_AFW','r_sqr_delta_AFW','r_sqr_c_AFT','r_sqr_delta_AFT'});
        writetable(T,pairwise_loci_filename);

        cell_data = num2cell(all_allele_pairs_table);
        T = cell2table(cell_data,'VariableNames',{'locus1','locus2','allele1','allele2','freq_allele_1','freq_allele_2','above_threhold','allele_freqs_product','r_sqr_c','r_sqr_delta'});
        writetable(T,pairwise_alleles_filename);
        
    end % if make_files
    
        
    % compute an average for S, the number of individuals for each estimate of r. There is one version for r^2_c and another for
    % r^2_delta because of div zero instances that may impact one but not the other.
    % Waples and Do 2008 used "weighted harmonic mean with weights proportional to the nij"

    sum_S = sum(S_c_values);
    mean_S_c = sum_S/c_locus_pairs;

    inv_S_c = 1./S_c_values;
    sum_inv_S_c = sum(inv_S_c);
    har_mean_S_c = 1/((1/c_locus_pairs)*sum_inv_S_c);

    sum_S_delta = sum(S_delta_values);
    mean_S_delta = sum_S_delta/delta_locus_pairs;

    inv_S_delta = 1./S_delta_values;
    sum_inv_S_delta = sum(inv_S_delta);
    har_mean_S_delta = 1/((1/delta_locus_pairs)*sum_inv_S_delta);
   
    % Estimate genetic association caused by finite sampling, see equations 13 and 14 in Sved et al. 2013. 
    correction_factor = (1/num_individuals)*(1-(1/((2*num_individuals - 1)^2))); %due to finite sampling of individuals
    correction_factor_avg_S = (1/mean_S_c)*(1-(1/((2*mean_S_c - 1)^2))); % arithmetic mean of sample sizes for locus pairs
    correction_factor_har_mean_S = (1/har_mean_S_c)*(1-(1/((2*har_mean_S_c - 1)^2))); % harmonic mean of sample sizes for locus pairs
   
    r_delta_correction_factor = (1/num_individuals)*(1-(1/((2*num_individuals - 1)^2))); %due to finite sampling of individuals
    r_delta_correction_factor_avg_S = (1/mean_S_delta)*(1-(1/((2*mean_S_delta - 1)^2))); % arithmetic mean of sample sizes for locus pairs
    r_delta_correction_factor_har_mean_S = (1/har_mean_S_delta)*(1-(1/((2*har_mean_S_delta - 1)^2))); % harmonic mean of sample sizes for locus pairs
    
   
    % make output table for sample sizes and correction factors
    fprintf('\n');
    fprintf('========================================================================================\n');
    fprintf('Table of sample sizes of diploid individuals (S) and\n');
    fprintf('E(r^2) for a finite sample of diploid individuals.\n\n');
    fprintf('Estimator\t\tS\tArithmetic mean S\tHarmonic mean S\n');
    fprintf('----------------------------------------------------------------------------------------\n');
    
    fprintf('r^2_{comp}\n');
    
    fprintf('\tSample size\t');
    fprintf('%8.2f\t',num_individuals);  
    fprintf('%8.2f\t',mean_S_c);
    fprintf('%8.2f\n',har_mean_S_c);
    
    fprintf('\tr^2 correction\t');
    fprintf('%8.6g\t',correction_factor);  
    fprintf('%8.6g\t',correction_factor_avg_S);
    fprintf('%8.6g\n',correction_factor_har_mean_S);
    
    fprintf('\n');
    
    fprintf('r^2_{delta}\n');

    fprintf('\tSample size\t');
    fprintf('%8.2f\t',num_individuals);  
    fprintf('%8.2f\t',mean_S_delta);
    fprintf('%8.2f\n',har_mean_S_delta);

    fprintf('\tr^2 correction\t');
    fprintf('%8.6g\t',r_delta_correction_factor);  
    fprintf('%8.6g\t',r_delta_correction_factor_avg_S);
    fprintf('%8.6g\n',r_delta_correction_factor_har_mean_S);

    fprintf('========================================================================================\n\n');
    
    % compute weighting factor for Weir's (1979) unbiased estimates of r_squared
    UB_weight_comp = har_mean_S_c/(har_mean_S_c - 1);
    UB_weight_delta = har_mean_S_delta/(har_mean_S_delta - 1);
    
    % make output table for r_squared values
    fprintf('\n');
    fprintf('=================================================================================================\n');
    fprintf('Overall r_squared values with multilocus weights based on numbers of individuals and numbers \n');
    fprintf('of alleles (eqn. 6 in Sved et al. 2013). For each locus pair estimates are either allele\n');
    fprintf('frequency weighted (AFW) or allele frequency thresholded (AFT).\n');
    fprintf('Unbiased estimates are multiplied by S/(S-1) where S is the harmonic mean sample size as in\n');
    fprintf('Weir (1979).\n');
    fprintf('-------------------------------------------------------------------------------------------------\n');
    fprintf('Estimates over all locus pairs:\n\n');

    fprintf('\tAllele frequency weighted (AFW) within locus pairs:\n');
    fprintf('\t\tr^2_{comp} = %6.6f\n',r_squared_c_W);
    fprintf('\t\tr^2_{delta} = %6.6f\n',r_squared_delta_W);
    
    ub_r_squared_c_W = UB_weight_comp*r_squared_c_W;
    ub_r_squared_delta_W = UB_weight_delta*r_squared_delta_W;
    
    fprintf('\n');
    fprintf('\t\tunbaised r^2_{comp} = %6.6f\n',ub_r_squared_c_W);
    fprintf('\t\tunbaised r^2_{delta} = %6.6f\n',ub_r_squared_delta_W);

    fprintf('\n');

    fprintf('\tAllele frequency thresholded (AFT) within locus pairs:\n');
    fprintf('\t\tr^2_{comp} = %6.6f\n', r_squared_c_TH);
    fprintf('\t\tr^2_{delta} = %6.6f\n', r_squared_delta_TH);
    
    ub_r_squared_c_TH = UB_weight_comp*r_squared_c_TH;
    ub_r_squared_delta_TH = UB_weight_delta*r_squared_delta_TH;
    
    fprintf('\n');
    fprintf('\t\tunbaised r^2_{comp} = %6.6f\n',ub_r_squared_c_TH);
    fprintf('\t\tunbaised r^2_{delta} = %6.6f\n',ub_r_squared_delta_TH);

    fprintf('=================================================================================================\n');
    fprintf('\n');
    
    
    % Provide some Ne estimates
    fprintf('=================================================================================================\n');
    fprintf('Estimates of Ne based on 1/(3 x r^2 - (r^2 correction factor)) with and without S/(S-1)\n');
    fprintf('bias correction to r^2. The r^2 correction factor is based on harmonic mean sample size of individuals, S.\n');
    fprintf('-------------------------------------------------------------------------------------------------\n');
    fprintf('Ne = %5.2f based on r^2_{c} AFW\n', 1/(3*(r_squared_c_W - correction_factor_har_mean_S)) );
    fprintf('Ne = %5.2f based on r^2_{c} AFT\n', 1/(3*(r_squared_c_TH - correction_factor_har_mean_S)) );
    fprintf('Ne = %5.2f based on r^2_{c} AFW unbiased\n', 1/(3*(ub_r_squared_c_W - correction_factor_har_mean_S)) );
    fprintf('Ne = %5.2f based on r^2_{c} AFT unbiased\n', 1/(3*(ub_r_squared_c_TH - correction_factor_har_mean_S)) );
    fprintf('\n');
    fprintf('Ne = %5.2f based on r^2_{delta} AFW\n', 1/(3*(r_squared_delta_W - r_delta_correction_factor_har_mean_S)) );
    fprintf('Ne = %5.2f based on r^2_{delta} AFT\n', 1/(3*(r_squared_delta_TH - r_delta_correction_factor_har_mean_S)) );
    fprintf('Ne = %5.2f based on r^2_{delta} AFW unbiased\n', 1/(3*(ub_r_squared_delta_W - r_delta_correction_factor_har_mean_S)) );
    fprintf('Ne = %5.2f based on r^2_{delta} AFT unbiased\n', 1/(3*(ub_r_squared_delta_TH - r_delta_correction_factor_har_mean_S)) );
    fprintf('=================================================================================================\n');
    
    
    % Compute Waples (2006) statistical-fit corrected Ne estimate.  
    % Equations based on Table 1 in Waples and Do (2007). Utilizes estimates of r_squared delta only. 
    fprintf('\n');
    fprintf('=================================================================================================\n');
    fprintf('Ne estimate based on statistical fit correction of Waples (2006) as implemented\n');
    fprintf('in LDNe and NeEstimator [see Waples & Do (2007) Table 1].\n');
    fprintf('Allele frequency thresholded (AFT) and allele frequency weighted (AFW) and well as \n');
    fprintf('estimates with and without S/(S-1) bias correction given for comparison. \n');
    fprintf('All estimates use harmonic mean sample size of individuals (S),\n');
    fprintf('A value of NaN indicates a denominator of zero.\n');
    fprintf('-------------------------------------------------------------------------------------------------\n');

    if har_mean_S_delta < 30
        fprintf('Estimate based on statistical fit for sample size < 30:\n');

        E_r2_delta_sample = 0.0018 + 0.907/har_mean_S_delta + 4.44/(har_mean_S_delta^2);

        r_2_prime_AFW = r_squared_delta_W - E_r2_delta_sample;
        
        if r_2_prime_AFW == 0
            % prevent div zero
            waples_AFW_Ne = NaN;   % assign NaN
        else
            waples_AFW_Ne = (0.308 + sqrt(0.308^2 + 2.08*r_2_prime_AFW^2) )/(2*r_2_prime_AFW); % Waples r_squared_delta AFW
        end
        fprintf('AFW Ne = %8.2f\n', waples_AFW_Ne);
        
        r_2_prime_AFT = r_squared_delta_TH - E_r2_delta_sample;
        if r_2_prime_AFT == 0
            % prevent div zero
            waples_AFT_Ne = NaN;   % assign NaN
        else
            waples_AFT_Ne = (0.308 + sqrt(0.308^2 + 2.08*r_2_prime_AFT^2) )/(2*r_2_prime_AFT); % Waples r_squared_delta AFT 
        end
        fprintf('AFT Ne = %8.2f\n', waples_AFT_Ne);

        r_2_prime_UB_AFW = ub_r_squared_delta_W - E_r2_delta_sample;
        if r_2_prime_UB_AFW == 0
            % prevent div zero
            waples_ub_AFW_Ne = NaN;   % assign NaN
        else
            waples_ub_AFW_Ne = (0.308 + sqrt(0.308^2 + 2.08*r_2_prime_UB_AFW^2) )/ (2*r_2_prime_UB_AFW); % Waples r_squared_delta AFW UB
        end
        fprintf('UB AFW Ne = %8.2f\n', waples_ub_AFW_Ne);

        r_2_prime_UB_AFT = ub_r_squared_delta_TH - E_r2_delta_sample;
        if r_2_prime_UB_AFT == 0
            % prevent div zero
            waples_ub_AFT_Ne = NaN;   % assign NaN
        else
            waples_ub_AFT_Ne = (0.308 + sqrt(0.308^2 + 2.08*r_2_prime_UB_AFT^2) )/(2*r_2_prime_UB_AFT); % Waples r_squared_delta AFT UB
        end
        fprintf('UB AFT Ne = %8.2f\n', waples_ub_AFT_Ne);

    else % har_mean_S_delta >= 30
        fprintf('Estimate based on statistical fit for sample size >= 30:\n');

        E_r2_delta_sample = 1/har_mean_S_delta + 3.19/(har_mean_S_delta^2);

        r_2_prime_AFW = r_squared_delta_W - E_r2_delta_sample;
        if r_2_prime_AFW == 0
            % prevent div zero
            waples_AFW_Ne = NaN;   % assign NaN
        else
            waples_AFW_Ne = (1/3 + sqrt(1/9 + 2.76*r_2_prime_AFW^2) )/(2*r_2_prime_AFW); % Waples r_squared_delta AFW
        end
        fprintf('AFW Ne = %8.2f\n', waples_AFW_Ne);
        
        r_2_prime_AFT = r_squared_delta_TH - E_r2_delta_sample;
        if r_2_prime_AFT == 0
            % prevent div zero
            waples_AFT_Ne = NaN;   % assign NaN
        else
            waples_AFT_Ne = (1/3 + sqrt(1/9 + 2.76*r_2_prime_AFT^2) )/(2*r_2_prime_AFT); % Waples r_squared_delta AFT 
        end
        fprintf('AFT Ne = %8.2f\n', waples_AFT_Ne);
        
        r_2_prime_UB_AFW = ub_r_squared_delta_W - E_r2_delta_sample;
        if r_2_prime_UB_AFW == 0
            % prevent div zero
            waples_ub_AFW_Ne = NaN;   % assign NaN
        else
            waples_ub_AFW_Ne = (1/3 + sqrt(1/9 + 2.76*r_2_prime_UB_AFW^2) )/ (2*r_2_prime_UB_AFW); % Waples r_squared_delta AFW UB
        end
        fprintf('UB AFW Ne = %8.2f\n', waples_ub_AFW_Ne);
        
        r_2_prime_UB_AFT = ub_r_squared_delta_TH - E_r2_delta_sample;
        if r_2_prime_UB_AFT == 0
            % prevent div zero
            waples_ub_AFT_Ne = NaN;   % assign NaN
        else
            waples_ub_AFT_Ne = (1/3 + sqrt(1/9 + 2.76*r_2_prime_UB_AFT^2) )/ (2*r_2_prime_UB_AFT); % Waples r_squared_delta AFT UB
        end
        fprintf('UB AFT Ne = %8.2f\n', waples_ub_AFT_Ne);
        
    end
    fprintf('=================================================================================================\n');

        
    % iteratively permute genotype data to estimate correction factor and associated Ne
    if permute == 1
        fprintf('\nPermuting genotype data %g times ...\n',n_permutations);
        make_table = true;
        make_files = false; % don't need to prepare output tables when permuting
        [~,~,~,~,median_r_sqr_c_AFW_perm,median_r_sqr_c_AFT_perm,~,~] = ...
            permute_data(genotypes,n_permutations,num_individuals,num_loci,threshold_allele_freq,r_squared_c_W,r_squared_c_TH,r_squared_delta_W,r_squared_delta_TH,alpha,make_table,make_graphs,make_files);

    end

    
    % Compute Waples & Do chi-square distribution CIs, "jackknife" CIs, and estimate of effective sample size n' and associated CIs.  
    if num_loci > 2
        print_output = true; % print results to console 
        
        [~,~,~,~,~] = wd_confidence_intervals_v2(num_loci,r_squared_c_values,r_squared_delta_values,c_locus_pairs,delta_locus_pairs,...
            r_squared_c_W,r_squared_c_TH,r_squared_delta_W,r_squared_delta_TH,correction_factor_har_mean_S,median_r_sqr_c_AFW_perm,median_r_sqr_c_AFT_perm,alpha,print_output,make_graphs);
    else
        fprintf('Too few loci to compute Waples & Do chi-square distribution CIs, "jackknife" CIs, and effective sample size n''.\n');
    end 
    
    
    if jack_individuals == 1
        % carry out jackknife across individuals
        if num_individuals > 2
            print_output = true; % print results to console
            make_files = false; % don't need to prepare output tables when jackknifing

            [~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~] = jackknife_individuals(r_squared_c_W,r_squared_delta_W,r_squared_c_TH,r_squared_delta_TH,median_r_sqr_c_AFW_perm,median_r_sqr_c_AFT_perm,...
                r_delta_correction_factor_har_mean_S,genotypes,num_individuals,num_loci,UB_weight_comp,UB_weight_delta,threshold_allele_freq,alpha,print_output,make_graphs,make_files);            
        else
            fprintf('\n');
            fprintf('*** Jackknife over individuals skipped. ***\n');
            fprintf('Three or more individuals required to carry out jackknife estimates of confidence intervals.\n');
        end

    end % if jack_individuals
        
    
    if jack_loci == 1
        % jackknife over loci
        
        if num_loci > 2
            print_output = true; % print results to console
            make_files = false; % don't need to prepare output tables when jackknifing

            [~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~] = jackknife_loci(r_squared_c_W,r_squared_delta_W,r_squared_c_TH,r_squared_delta_TH,median_r_sqr_c_AFW_perm,median_r_sqr_c_AFT_perm,...
                r_delta_correction_factor_har_mean_S,genotypes,num_individuals,num_loci,UB_weight_comp,UB_weight_delta,threshold_allele_freq,alpha,print_output,make_graphs,make_files);            
        else
            fprintf('\n');
            fprintf('*** Jackknife over loci skipped. ***\n');
            fprintf('Three or more loci required to carry out jackknife estimates of confidence intervals.\n');
        end
        
    end % if jack_loci
        
    fprintf('\n');
    fprintf('Program completed.\n');
    
    

    
% --- Executes on button press in data_file_button.
function data_file_button_Callback(hObject, eventdata, handles)

    [filename,filepath] = uigetfile({'*.*','All Files'},'Select input file of genotype data');
    
    full_path = strcat(filepath,filename);
    
    %assign value to UserData field
    set(handles.data_file_string,'String',full_path); 

    % Update handles structure to show text
    guidata(hObject, handles);


function min_allele_freq_Callback(hObject, eventdata, handles)

    threshold = str2double(get(hObject,'String'));
    
    threshold_min = 0.0;
    threshold_max = 1.0;
    
    if threshold < threshold_min
        threshold = threshold_min;  %set lower limit value
        textval = num2str(threshold_min);
        set(handles.min_allele_freq,'String',textval);   % change text to give feedback
    end
    
    if threshold > threshold_max
        threshold = threshold_max;  %set lower limit value
        textval = num2str(threshold_max);
        set(handles.min_allele_freq,'String',textval);   % change text to give feedback
    end
    
    %assign value to UserData field
    set(handles.min_allele_freq,'UserData',threshold); %UserData = threshold

    % Update handles structure to show text
    guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function min_allele_freq_CreateFcn(hObject, eventdata, handles)
% hObject    handle to min_allele_freq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function data_permute_reps_Callback(hObject, eventdata, handles)

    perm_reps = str2double(get(hObject,'String'));
    
    perm_reps_min = 10;
    perm_reps_max = 100000;
    
    if perm_reps < perm_reps_min
        perm_reps = perm_reps_min;  %set lower limit value
        textval = num2str(perm_reps_min);
        set(handles.data_permute_reps,'String',textval);   % change text to give feedback
    end
    
    if perm_reps > perm_reps_max
        perm_reps = perm_reps_max;  %set lower limit value
        textval = num2str(perm_reps_max);
        set(handles.data_permute_reps,'String',textval);   % change text to give feedback
    end
    
    %assign value to UserData field
    set(handles.data_permute_reps,'UserData',perm_reps); %UserData = threshold

    % Update handles structure to show text
    guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function data_permute_reps_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end


function missing_data_button_group_SelectionChangeFcn(hObject, eventdata)
    %UserData = missing
        
    %retrieve GUI data, i.e. the handles structure
    handles = guidata(hObject); 

    switch get(eventdata.NewValue,'Tag')   % Get Tag of selected object
        case 'both_missing'
          %execute this code when radiobutton1 is selected
          missing = 1; % when one allele is missing code genotype as missing

        case 'use_one_allele'
          %execute this code when radiobutton2 is selected
          missing = 2; % use the one scored allele and treat other as missing

        otherwise
           % Code for when there is no match
           missing = 1; % when one allele is missing code genotype as missing

    end
    
    %assign to UserData field
    set(handles.missing_data_button_group,'UserData',missing); %UserData = missing


function alpha_Callback(hObject, eventdata, handles)

    alpha = str2double(get(hObject,'String'));
    
    alpha_min = 0.0;
    alpha_max = 1.0;
    
    if alpha < alpha_min
        alpha = alpha_min;  %set lower limit value
        textval = num2str(alpha_min);
        set(handles.alpha,'String',textval);   % change text to give feedback
    end
    
    if alpha > alpha_max
        alpha = alpha_max;  %set lower limit value
        textval = num2str(alpha_max);
        set(handles.alpha,'String',textval);   % change text to give feedback
    end
    
    %assign value to UserData field
    set(handles.alpha,'UserData',alpha); %UserData = alpha

    % Update handles structure to show text
    guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function alpha_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end

    
% --- Executes on button press in jackknife_individs_checkbox.
function jackknife_individs_checkbox_Callback(hObject, eventdata, handles)

    if (get(hObject,'Value') == get(hObject,'Max'))
        % box is checked
        temp = get(handles.jackknife_individs_checkbox,'Value');
        
        %assign one to UserData field
        set(handles.jackknife_individs_checkbox,'UserData',1);
        
    else
        % box is unchecked
        temp = get(handles.jackknife_individs_checkbox,'Value');
        
        %assign zero to UserData field
        set(handles.jackknife_individs_checkbox,'UserData',0); 
        
    end
    
    
    % --- Executes during object creation, after setting all properties.
function jackknife_individs_checkbox_CreateFcn(hObject, eventdata, handles)

    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end


    
    % --- Executes on button press in jackknife_individs_checkbox.
function jackknife_loci_checkbox_Callback(hObject, eventdata, handles)

    if (get(hObject,'Value') == get(hObject,'Max'))
        % box is checked
        temp = get(handles.jackknife_loci_checkbox,'Value');
        
        %assign one to UserData field
        set(handles.jackknife_loci_checkbox,'UserData',1);
        
    else
        % box is unchecked
        temp = get(handles.jackknife_loci_checkbox,'Value');
        
        %assign zero to UserData field
        set(handles.jackknife_loci_checkbox,'UserData',0); 
        
    end
    
    
    % --- Executes during object creation, after setting all properties.
function jackknife_loci_checkbox_CreateFcn(hObject, eventdata, handles)

    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end

    
    
% --- Executes on button press in permute_checkbox.
function permute_checkbox_Callback(hObject, eventdata, handles)

    if (get(hObject,'Value') == get(hObject,'Max'))
        % box is checked
        temp = get(handles.permute_checkbox,'Value');
        
        %assign one to UserData field
        set(handles.permute_checkbox,'UserData',1);
        
        % make permute_reps text box active
        set(handles.data_permute_reps, 'enable', 'on');
    else
        % box is unchecked
        temp = get(handles.permute_checkbox,'Value');
        
        %assign zero to UserData field
        set(handles.permute_checkbox,'UserData',0); 
        
        % make permute_reps text box inactive
        set(handles.data_permute_reps, 'enable', 'off');
    end


% --- Executes on button press in paired_loci_tables_checkbox.
function paired_loci_tables_checkbox_Callback(hObject, eventdata, handles)

    if (get(hObject,'Value') == get(hObject,'Max'))
        % box is checked
        temp = get(handles.paired_loci_tables_checkbox,'Value');
        
        %assign one to UserData field
        set(handles.paired_loci_tables_checkbox,'UserData',1);
    else
        % box is unchecked
        temp = get(handles.paired_loci_tables_checkbox,'Value');
        
        %assign zero to UserData field
        set(handles.paired_loci_tables_checkbox,'UserData',0); 
    end


% --- Executes on button press in allele_freq_table_checkbox.
function allele_freq_table_checkbox_Callback(hObject, eventdata, handles)

    if (get(hObject,'Value') == get(hObject,'Max'))
        % box is checked
        temp = get(handles.allele_freq_table_checkbox,'Value');
        
        %assign one to UserData field
        set(handles.allele_freq_table_checkbox,'UserData',1);
    else
        % box is unchecked
        temp = get(handles.allele_freq_table_checkbox,'Value');
        
        %assign zero to UserData field
        set(handles.allele_freq_table_checkbox,'UserData',0); 
    end
    
    
    
% --- Executes on button press in output_graphs_checkbox.
function output_graphs_checkbox_Callback(hObject, eventdata, handles)

    if (get(hObject,'Value') == get(hObject,'Max'))
        % box is checked
        temp = get(handles.output_graphs_checkbox,'Value');
        
        %assign one to UserData field
        set(handles.output_graphs_checkbox,'UserData',1);
    else
        % box is unchecked
        temp = get(handles.output_graphs_checkbox,'Value');
        
        %assign zero to UserData field
        set(handles.output_graphs_checkbox,'UserData',0); 
    end


    
% --- Executes on button press in file_output_checkbox.
function file_output_checkbox_Callback(hObject, eventdata, handles)

    if (get(hObject,'Value') == get(hObject,'Max'))
        % box is checked
        temp = get(handles.file_output_checkbox,'Value');
        
        %assign one to UserData field
        set(handles.file_output_checkbox,'UserData',1);
    else
        % box is unchecked
        temp = get(handles.file_output_checkbox,'Value');
        
        %assign zero to UserData field
        set(handles.file_output_checkbox,'UserData',0); 
    end
