%% Connectivity Analysis
%  ------------------------------------------------------------------------
%
% Script Description:
% This MATLAB script performs connectivity analysis of preprocessed EEG data.
%
% Usage:
%   - This script requires the following MATLAB toolboxes:
%       1. Image Processing Toolbox
%       2. Statistics and Machine Learning Toolbox
%       3. Bioinformatics Toolbox
%       4. FieldTrip Toolbox
%   - Note: check which toolboxes are required by [1] running script [2] executing: [fList,pList] = matlab.codetools.requiredFilesAndProducts('run_analysis.m'); disp({pList.Name}');
%
% Inputs:
%   - Raw data is stored in '/EEG-analysis/input' and includes:
%       1. The subdirectory '/EEG_data', which contains the unprocessed EEG data for each subject (.eeg, .vhdr and .vmrk files)
%       2. The subdirectory '/event_codes', which contains the condition codes (1-12) for each trial and for each subject
%
% Outputs:
%   - PLI matrices are saved to '/EEG-analysis/output'
%   - Mean PLIs of interest are saved to '/EEG-analysis/statistics' for statistical analysis 
%
% Author:
%   Jane Chesley
%
% Affiliation:
%   Faculty of Psychology & Neuroscience, Maastricht University
%   Maastricht, the Netherlands
%
% Date:
%   22 February, 2024 
%
% Outline of script content:
%   1. Script setup
%   2. Load data and define parameters  
%   3. Simulations
% 		3.1 Construct simulated data 
%		3.2 Compute PLIs for all channel pairs 
% 		3.3 Compute within-region PLIs 
% 	4. Real data 
% 		4.1 Real data: Compute PLIs for all channel pairs 
% 		4.2 Real data: Compute within-region PLIs for pre-defined ROIs
% 		4.3 Real data: Compute between-region PLIs for pre-defined ROIs
%       4.4 Real data: PLI data-driven ROI selection
%   
% ------------------------------------------------------------------------



%% ------------------------------------------------------------------------
%  Part 1 - Script setup
%  ------------------------------------------------------------------------

% clean working environment 
clear, clc, close all;
restoredefaultpath;

% set up relevant directories 
[dir_parent, ~, ~] = fileparts(pwd);
cd(dir_parent);

% add custom functions to path
addpath('functions');

% specify the path to the Fieldtrip Toolbox on your system
dir_FT = '/Users/jane_chesley/Library/Application Support/MathWorks/MATLAB Add-Ons/Collections/FieldTrip';

% initialize Fieldtrip toolbox
addpath(dir_FT); ft_defaults;



%% ------------------------------------------------------------------------
%  Part 2 - Load data and define parameters 
%  ------------------------------------------------------------------------

% identify all input data to analyze
% data to analyze are preprocessed EEG data, stratified by condition
input = dir(fullfile('output', '*data_clean*'));

% load all files containing 'data_clean' in their names 
disp('loading data...')
for i = 1: length(input)
    disp(strcat('loading data...', num2str(i),' of ', num2str(length(input))));
    load(strcat(input(i).folder, '/',input(i).name));
end

% define time-window of interest (TOI) in seconds
min_T1 = 0; % TOI minimum 
max_T1 = 1; % TOI maximum 

% label all frequency bands of interest 
% the stimuli used for the present experiment contain motion artifacts at 30 Hz
% to handle this, gamma is analyzed twice: once at 35-45 Hz and once at 24-45 Hz with a bandstop at 30 Hz 
all_freq = {'delta';'theta';'alpha';'beta';'gamma_A_35_45';'gamma_B_25_45_BS30';'broadband_1_45'};

% identify all channel labels 
% channel labels are constant across all measurements
all_channels = data_clean_hum_body_norm{1}.label; 



%% ------------------------------------------------------------------------
%  Part 3 - Simulations
%  ------------------------------------------------------------------------
% Create simulated data and run analyses on it to confirm the calculations
% and scripts are correct 

%  ------------------------------------------------------------------------
%  Part 3.1 Simulations: Construct simulated data 
%  ------------------------------------------------------------------------

% generate three sine waves with very small phase offset
% expect high phase coherence (PLI = 1) for these three signals
frequency = 10; % frequency of the sine wave
amplitude = 5; % amplitude of the sine wave
timepts = linspace(0,575,575);
sin_wave1 = amplitude * sin(2*pi*frequency*timepts);
sin_wave2 = amplitude * sin(2*pi*frequency*timepts + pi/256);
sin_wave3 = amplitude * sin(2*pi*frequency*timepts + pi/128);



% create simulated data for each condition, which contains real data from three example subjects
sim_data_clean_hum_body_norm = data_clean_hum_body_norm(1:3);
sim_data_clean_hum_body_scr = data_clean_hum_body_scr(1:3);
sim_data_clean_hum_face_norm = data_clean_hum_face_norm(1:3);
sim_data_clean_hum_face_scr = data_clean_hum_face_scr(1:3);
sim_data_clean_hum_obj_norm = data_clean_hum_obj_norm(1:3);
sim_data_clean_hum_obj_scr = data_clean_hum_obj_scr(1:3);

sim_data_clean_monk_body_norm = data_clean_monk_body_norm(1:3);
sim_data_clean_monk_body_scr = data_clean_monk_body_scr(1:3);
sim_data_clean_monk_face_norm = data_clean_monk_face_norm(1:3);
sim_data_clean_monk_face_scr = data_clean_monk_face_scr(1:3);
sim_data_clean_monk_obj_norm = data_clean_monk_obj_norm(1:3);
sim_data_clean_monk_obj_scr = data_clean_monk_obj_scr(1:3);



% inject simulated data into the real data 
% assign three sine waves to channels 1-3, respectively, for all trial 
% and all subject data for one condition ('human body normal')
total_subjects = 3;
for subject = 1:total_subjects

    total_trials = length(sim_data_clean_hum_body_norm{subject}.trial);

    for trial = 1:total_trials
        sim_data_clean_hum_body_norm{subject}.trial{trial}(1,:) = sin_wave1; % channel 1, all timepts = sin_wave1
        sim_data_clean_hum_body_norm{subject}.trial{trial}(2,:) = sin_wave2; % channel 2, all timepts = sin_wave2
        sim_data_clean_hum_body_norm{subject}.trial{trial}(3,:) = sin_wave3; % channel 3, all timepts = sin_wave3
    end

end



%% ------------------------------------------------------------------------
%  Part 3.2 Simulations: Compute PLIs for all channel pairs
%  ------------------------------------------------------------------------

sim_PLI_hum_body_norm = PLI_all_subjects(sim_data_clean_hum_body_norm, min_T1, max_T1);
sim_PLI_hum_face_norm = PLI_all_subjects(sim_data_clean_hum_face_norm, min_T1, max_T1);
sim_PLI_hum_obj_norm = PLI_all_subjects(sim_data_clean_hum_obj_norm, min_T1, max_T1);

sim_PLI_monk_body_norm = PLI_all_subjects(sim_data_clean_monk_body_norm, min_T1, max_T1);
sim_PLI_monk_face_norm = PLI_all_subjects(sim_data_clean_monk_face_norm, min_T1, max_T1);
sim_PLI_monk_obj_norm = PLI_all_subjects(sim_data_clean_monk_obj_norm, min_T1, max_T1);

sim_PLI_hum_body_scr = PLI_all_subjects(sim_data_clean_hum_body_scr, min_T1, max_T1);
sim_PLI_hum_face_scr = PLI_all_subjects(sim_data_clean_hum_face_scr, min_T1, max_T1);
sim_PLI_hum_obj_scr = PLI_all_subjects(sim_data_clean_hum_obj_scr, min_T1, max_T1);

sim_PLI_monk_body_scr = PLI_all_subjects(sim_data_clean_monk_body_scr, min_T1, max_T1);
sim_PLI_monk_face_scr = PLI_all_subjects(sim_data_clean_monk_face_scr, min_T1, max_T1);
sim_PLI_monk_obj_scr = PLI_all_subjects(sim_data_clean_monk_obj_scr, min_T1, max_T1);



%% ------------------------------------------------------------------------
%  Part 3.3 Simulations: Compute within-region PLIs 
%  ------------------------------------------------------------------------

% pre-allocate variables
sim_PLI_normalized = zeros(3,6); % 3 subjects x 6 conditions
sim_data_long = cell(1,3); % 3 subjects 

% compute within-region PLIs for a region of interest (ROI), which contains simulated data 
roi = {'AFz','Fz','FCz'}; % define channel labels
roi_string = 'simRoi'; % define cluster label 
roi_idx = find(ismember(all_channels, roi)); % extract indices of ROI


% compute PLIs separately for each frequency band
for freq = 1:length(all_freq)

    freq_string = all_freq{freq};

    % compute within-region connectivity for each subject and condition 
    % sim_PLI_allconditions = subject (3) x condition (12) array 
    % 12 conditions:
    %   Human body normal; Human face normal; Human object normal 
    %   Monkey body normal; Monkey face normal; Monkey object normal 
    %   Human body scramble; Human face scramble; Human object scramble 
    %   Monkey body scramble; Monkey face scramble; Monkey object scramble 

    sim_PLI_allconditions(:,1) = PLI_within(sim_PLI_hum_body_norm, roi_idx, freq);
    sim_PLI_allconditions(:,2) = PLI_within(sim_PLI_hum_face_norm, roi_idx, freq);
    sim_PLI_allconditions(:,3) = PLI_within(sim_PLI_hum_obj_norm, roi_idx, freq);

    sim_PLI_allconditions(:,4) = PLI_within(sim_PLI_monk_body_norm, roi_idx, freq);
    sim_PLI_allconditions(:,5) = PLI_within(sim_PLI_monk_face_norm, roi_idx, freq);
    sim_PLI_allconditions(:,6) = PLI_within(sim_PLI_monk_obj_norm, roi_idx, freq);

    sim_PLI_allconditions(:,7) = PLI_within(sim_PLI_hum_body_scr, roi_idx, freq);
    sim_PLI_allconditions(:,8) = PLI_within(sim_PLI_hum_face_scr, roi_idx, freq);
    sim_PLI_allconditions(:,9) = PLI_within(sim_PLI_hum_obj_scr, roi_idx, freq);

    sim_PLI_allconditions(:,10) = PLI_within(sim_PLI_monk_body_scr, roi_idx, freq);
    sim_PLI_allconditions(:,11) = PLI_within(sim_PLI_monk_face_scr, roi_idx, freq);
    sim_PLI_allconditions(:,12) = PLI_within(sim_PLI_monk_obj_scr, roi_idx, freq);

    % normalize PLIs
    % the design for the present experiment contains scrambled conditions to control for low-level visual features 
    % to handle this, the the arithmetic subtraction of normal - scramble conditions normalizes the data
    % sim_PLI_normalized = subject (3) x condition (6) array
    % 6 conditions:
    %   Human body normal-scramble; Human face normal-scramble; Human object normal-scramble 
    %   Monkey body normal-scramble; Monkey face normal-scramble; Monkey object normal-scramble 

    for c = 1:6
        sim_PLI_normalized(:,c) = sim_PLI_allconditions(:,c)-sim_PLI_allconditions(:,c+6);
    end

    % transform data from wide to long format for statistical analysis in R 
    % wide format is organized as: subjects (rows) x conditions (columns)
    % long format contains one experimental observation per row, organized as follows:
    % data_long{s}(:,subject | IV1(Species: 1.Human 2.Monkey) | IV2(Category: 1.Body 2.Face 3.Object) | Measurement(normalized PLI))
    % IV = Independent Variable

    for s = 1:total_subjects
        
        % extract measurements for each subject and condition 
        sim_data_long{s}(1,4) = sim_PLI_normalized(s,1); % human body condition 
        sim_data_long{s}(2,4) = sim_PLI_normalized(s,2); % human face condition 
        sim_data_long{s}(3,4) = sim_PLI_normalized(s,3); % human object condition 

        sim_data_long{s}(4,4) = sim_PLI_normalized(s,4); % monkey body condition 
        sim_data_long{s}(5,4) = sim_PLI_normalized(s,5); % monkey face condition 
        sim_data_long{s}(6,4) = sim_PLI_normalized(s,6); % monkey object condition 
        
        % label subject number
        sim_data_long{s}(:,1) = s; % subject number
        
        % label condition levels for IV1
        sim_data_long{s}(1:3,2) = 1; % IV1(Species): human(1)
        sim_data_long{s}(4:6,2) = 2; % IV1(Species): monkey(2)

        % label condition levels for IV2
        sim_data_long{s}([1,4],3) = 1; % IV2(Category): body(1)
        sim_data_long{s}([2,5],3) = 2; % IV2(Category): face(2)
        sim_data_long{s}([3,6],3) = 3; % IV2(Category): object(3)

    end

    % create a table of subject-level data in *long format*
    sim_data_allsubj = cat(1, sim_data_long{:}); % concatenate all subject-level data
    sim_variable_names = {'Subject' 'Species','Category','Measurement'}; % assign variable names 
    sim_T = array2table(sim_data_allsubj, 'VariableNames', sim_variable_names); % create a table 

    % save data as a table (.xlsx) in the 'statistics/R’ directory, with its file name and folder name representing the analysis
    % data in *long format* are prepared for statistical analysis in *R*

    analysis = strcat('sim_PLI_normalized_',roi_string,'_Within-PLI_',freq_string); % define present analysis
    folder_path = strcat('statistics/R/',analysis);
    file_name = strcat(analysis,'.xlsx');
    full_file_path = strcat(folder_path,'/',file_name);
    
    % if the folder does not exist, then create it
    if ~isfolder(folder_path)
        mkdir(folder_path);
    end

    % if the file name already exists, then delete it 
    if exist(full_file_path, 'file') 
        delete(full_file_path); 
    end

    writetable(sim_T,full_file_path); % save

end



%% ------------------------------------------------------------------------
%  Part 4 - Real data
%  ------------------------------------------------------------------------
% Statistical analyses were run in R and confirmed the expected results of the simulations
% Now, the same PLI pipeline is applied to real data 
% Real data are preprocessed EEG data, stratified by condition



%% ------------------------------------------------------------------------
%  Part 4.1 - Real data: Compute PLIs for all channel pairs 
%  ------------------------------------------------------------------------

% create a log file to document script progress, as these calculations can take several hours 
logFileName = 'logfile.txt';
logfile = fopen(logFileName, 'a'); % open or create a log file for writing (append mode)
fprintf(logfile, 'Script started\n'); % log a message indicating the start of the script


% calculate subject-level PLIs separately for each frequency band, averaged across trials
% input for 'PLI_all_subjects.m' is subject-level clean data for one condition
% output for 'PLI_all_subjects.m' are 4-D PLI data (freq x channel x channel x subject)
% E.g. PLI_hum_body_norm{freq}(channel,channel,subject);
% PLI values are represented in cells of each channel x channel matrix 
% PLI values are between 0-1 (0 = no connectivity; 1 = perfect connectivity)
% frequency bands are segmented as follows: {freq1: delta 0.5-3.9 Hz} {freq2: theta 4-7 Hz} {freq3: alpha 8-12 Hz}
% {freq4: beta 13-30 Hz} {freq5: gamma_A 35-45 Hz} {freq6: gamma_B 25-45 Hz + bandstop 29:31 Hz}

PLI_hum_body_norm = PLI_all_subjects(data_clean_hum_body_norm, min_T1, max_T1);
fprintf(logfile, 'PLI 1/14 completed\n');
save('output/PLI_hum_body_norm.mat','PLI_hum_body_norm')

PLI_hum_face_norm = PLI_all_subjects(data_clean_hum_face_norm, min_T1, max_T1);
fprintf(logfile, 'PLI 2/14 completed\n');
save('output/PLI_hum_face_norm.mat','PLI_hum_face_norm')

PLI_hum_obj_norm = PLI_all_subjects(data_clean_hum_obj_norm, min_T1, max_T1);
fprintf(logfile, 'PLI 3/14 completed\n');
save('output/PLI_hum_obj_norm.mat','PLI_hum_obj_norm')

PLI_monk_body_norm = PLI_all_subjects(data_clean_monk_body_norm, min_T1, max_T1);
fprintf(logfile, 'PLI 4/14 completed\n');
save('output/PLI_monk_body_norm.mat','PLI_monk_body_norm')

PLI_monk_face_norm = PLI_all_subjects(data_clean_monk_face_norm, min_T1, max_T1);
fprintf(logfile, 'PLI 5/14 completed\n');
save('output/PLI_monk_face_norm.mat','PLI_monk_face_norm')

PLI_monk_obj_norm = PLI_all_subjects(data_clean_monk_obj_norm, min_T1, max_T1);
fprintf(logfile, 'PLI 6/14 completed\n');
save('output/PLI_monk_obj_norm.mat','PLI_monk_obj_norm')

PLI_hum_body_scr = PLI_all_subjects(data_clean_hum_body_scr, min_T1, max_T1);
fprintf(logfile, 'PLI 7/14 completed\n');
save('output/PLI_hum_body_scr.mat','PLI_hum_body_scr')

PLI_hum_face_scr = PLI_all_subjects(data_clean_hum_face_scr, min_T1, max_T1);
fprintf(logfile, 'PLI 8/14 completed\n');
save('output/PLI_hum_face_scr.mat','PLI_hum_face_scr')

PLI_hum_obj_scr = PLI_all_subjects(data_clean_hum_obj_scr, min_T1, max_T1);
fprintf(logfile, 'PLI 9/14 completed\n');
save('output/PLI_hum_obj_scr.mat','PLI_hum_obj_scr')

PLI_monk_body_scr = PLI_all_subjects(data_clean_monk_body_scr, min_T1, max_T1);
fprintf(logfile, 'PLI 10/14 completed\n');
save('output/PLI_monk_body_scr.mat','PLI_monk_body_scr')

PLI_monk_face_scr = PLI_all_subjects(data_clean_monk_face_scr, min_T1, max_T1);
fprintf(logfile, 'PLI 11/14 completed\n');
save('output/PLI_monk_face_scr.mat','PLI_monk_face_scr')

PLI_monk_obj_scr = PLI_all_subjects(data_clean_monk_obj_scr, min_T1, max_T1);
fprintf(logfile, 'PLI 12/14 completed\n');
save('output/PLI_monk_obj_scr.mat','PLI_monk_obj_scr')

PLI_pooled_normal = PLI_all_subjects(data_clean_pooled_normal, min_T1, max_T1);
fprintf(logfile, 'PLI 13/14 completed\n');
save('output/PLI_pooled_normal.mat','PLI_pooled_normal')

PLI_pooled_scramble = PLI_all_subjects(data_clean_pooled_scramble, min_T1, max_T1);
fprintf(logfile, 'PLI 14/14 completed\n');
save('output/PLI_pooled_scramble.mat','PLI_pooled_scramble')

fclose(logfile); % close the log file



%% ------------------------------------------------------------------------
%  Part 4.2 - Real data: Compute within-region PLIs 
%  ------------------------------------------------------------------------
% Within-region PLI refers to the connectivity within a specified cluster of channels 

% load input data
% input data are 4-D PLI data (freq x channel x channel x subject) for each condition, computed above in step 4.1
input = dir(fullfile('output','*PLI*'));
for i = 1: length(input)
    load(strcat(input(i).folder, '/',input(i).name)); % load all files containing 'PLI' in their name
end



% pre-allocate variables
total_subjects = 29; 
PLI_normalized = zeros(total_subjects,6); % subjects (29) x conditions (6)
data_long = cell(1,total_subjects);
data_long2 = cell(1,total_subjects);



% compute within-region PLIs, separately for two ROIs 
% these two ROIs (clusters) are selected based off of previous research (see Chesley et al., 2024, doi: https://doi.org/10.1162/imag_a_00150)
for cluster = 1:2
    
    if cluster == 1
        roi = {'C3', 'CP3', 'P3'}; % define channel labels
        roi_string = 'roiHumanBody'; % 'human body selective cluster'

    elseif cluster == 2
        roi = {'AFz', 'FCz', 'Cz', 'CPz', 'Pz', 'Fp1', 'Fp2', 'F3', 'F4', 'F7', 'F8', 'FC3', 'FC4', 'FT7', 'FT8', 'C3', 'C4', 'CP3', 'CP4', 'TP10', 'P3', 'P4', 'P8'} ; % define channel labels 
        roi_string = 'roiObjectLvl'; % 'object-level processing cluster'
    end

    % extract indices of ROI
    roi_idx = find(ismember(all_channels, roi));



    % compute PLIs separately for each frequency band
    for freq = 1:length(all_freq)

        freq_string = all_freq{freq};

        % compute within-region connectivity for each subject and condition
        % sim_PLI_allconditions = subject (3) x condition (12) array
        % 12 conditions:
        %   Human body normal; Human face normal; Human object normal
        %   Monkey body normal; Monkey face normal; Monkey object normal
        %   Human body scramble; Human face scramble; Human object scramble
        %   Monkey body scramble; Monkey face scramble; Monkey object scramble

        PLI_allconditions(:,1) = PLI_within(PLI_hum_body_norm, roi_idx, freq);
        PLI_allconditions(:,2) = PLI_within(PLI_hum_face_norm, roi_idx, freq);
        PLI_allconditions(:,3) = PLI_within(PLI_hum_obj_norm, roi_idx, freq);

        PLI_allconditions(:,4) = PLI_within(PLI_monk_body_norm, roi_idx, freq);
        PLI_allconditions(:,5) = PLI_within(PLI_monk_face_norm, roi_idx, freq);
        PLI_allconditions(:,6) = PLI_within(PLI_monk_obj_norm, roi_idx, freq);

        PLI_allconditions(:,7) = PLI_within(PLI_hum_body_scr, roi_idx, freq);
        PLI_allconditions(:,8) = PLI_within(PLI_hum_face_scr, roi_idx, freq);
        PLI_allconditions(:,9) = PLI_within(PLI_hum_obj_scr, roi_idx, freq);

        PLI_allconditions(:,10) = PLI_within(PLI_monk_body_scr, roi_idx, freq);
        PLI_allconditions(:,11) = PLI_within(PLI_monk_face_scr, roi_idx, freq);
        PLI_allconditions(:,12) = PLI_within(PLI_monk_obj_scr, roi_idx, freq);


        % create a table of subject-level data in *wide format* 
        variable_names = {'hum_body_norm','hum_face_norm','hum_obj_norm','monk_body_norm','monk_face_norm','monk_obj_norm', 'hum_body_scr','hum_face_scr','hum_obj_scr','monk_body_scr','monk_face_scr','monk_obj_scr'};
        T = array2table(PLI_allconditions, 'VariableNames', variable_names);
        
        % save data as a table (.xlsx) in the 'statistics/SPSS' directory, with its file name and folder name representing the analysis
        % data in *wide format* are prepared for statistical analysis in *SPSS*
        analysis = strcat('PLI_allcond_',roi_string,'_Within-PLI_',freq_string);
        folder_path = strcat('statistics/SPSS/',analysis);
        file_name = strcat(analysis,'.xlsx');
        full_file_path = strcat(folder_path,'/',file_name);

        % if the folder does not exist, then create it
        if ~isfolder(folder_path)
            mkdir(folder_path);
        end
        
        % if the file name already exists, delete it
        if exist(full_file_path, 'file') 
            delete(full_file_path); 
        end
        
        writetable(T,full_file_path); % save



        % normalize PLIs
        % the design for the present experiment contains scrambled conditions to control for low-level visual features
        % to handle this, the the arithmetic subtraction of normal - scramble conditions normalizes the data
        % sim_PLI_normalized = subject (3) x condition (6) array
        % 6 conditions:
        %   Human body normal-scramble; Human face normal-scramble; Human object normal-scramble
        %   Monkey body normal-scramble; Monkey face normal-scramble; Monkey object normal-scramble
        for c = 1:6
            PLI_normalized(:,c) = PLI_allconditions(:,c)-PLI_allconditions(:,c+6);
        end



        % create a table of subject-level data in *wide format* 
        variable_names = {'hum_body','hum_face','hum_obj','monk_body','monk_face','monk_obj'};
        T = array2table(PLI_normalized, 'VariableNames', variable_names);
        
        % save data as a table (.xlsx) in the 'statistics/SPSS' directory, with its file name and folder name representing the analysis
        % data in *wide format* are prepared for statistical analysis in *SPSS*
        analysis = strcat('PLI_normalized_',roi_string,'_Within-PLI_',freq_string);
        folder_path = strcat('statistics/SPSS/',analysis);
        file_name = strcat(analysis,'.xlsx');
        full_file_path = strcat(folder_path,'/',file_name);

        % if the folder does not exist, then create it
        if ~isfolder(folder_path)
            mkdir(folder_path);
        end

        % if the file name already exists, delete it
        if exist(full_file_path, 'file') 
            delete(full_file_path); 
        end
        
        writetable(T,full_file_path); % save



        % transform non-normalized data from wide to long format for statistical analysis in R
        % wide format is organized as: subjects (rows) x conditions (columns)
        % long format contains one experimental observation per row, organized as follows:
        % data_long{s}(:,subject | IV1(Species: 1.Human 2.Monkey) | IV2(Category: 1.Body 2.Face 3.Object) | IV3(Configuration: 1.Normal 2.Scramble | Measurement(non-normalized PLI))

        for s = 1:total_subjects
            
            % extract measurements for each subject and condition 
            data_long{s}(1,5) = PLI_allconditions(s,1); % human body normal 
            data_long{s}(2,5) = PLI_allconditions(s,2); % human face normal 
            data_long{s}(3,5) = PLI_allconditions(s,3); % human object normal

            data_long{s}(4,5) = PLI_allconditions(s,4); % monkey body normal 
            data_long{s}(5,5) = PLI_allconditions(s,5); % monkey face normal 
            data_long{s}(6,5) = PLI_allconditions(s,6); % monkey object normal 

            data_long{s}(7,5) = PLI_allconditions(s,7); % human body scramble
            data_long{s}(8,5) = PLI_allconditions(s,8); % human face scramble 
            data_long{s}(9,5) = PLI_allconditions(s,9); % human object scramble

            data_long{s}(10,5) = PLI_allconditions(s,10); % monkey body scramble 
            data_long{s}(11,5) = PLI_allconditions(s,11); % monkey face scramble 
            data_long{s}(12,5) = PLI_allconditions(s,12); % monkey object scramble 


            % label subject number
            data_long{s}(:,1) = s; 
            
            % label condition levels for IV1
            data_long{s}(1:3,2) = 1; % IV1(Species): human(1)
            data_long{s}(7:9,2) = 1; % IV1(Species): human(1)

            data_long{s}(4:6,2) = 2; % IV1(Species): monkey(2)
            data_long{s}(10:12,2) = 2; % IV1(Species): monkey(2)
            
            % label condition levels for IV2
            data_long{s}(1,3) = 1; % IV2(Category): body(1)
            data_long{s}(2,3) = 2; % IV2(Category): face(2)
            data_long{s}(3,3) = 3; % IV2(Category): object(3)

            data_long{s}(4,3) = 1; % IV2(Category): body(1)
            data_long{s}(5,3) = 2; % IV2(Category): face(2)
            data_long{s}(6,3) = 3; % IV2(Category): object(3)

            data_long{s}(7,3) = 1; % IV2(Category): body(1)
            data_long{s}(8,3) = 2; % IV2(Category): face(2)
            data_long{s}(9,3) = 3; % IV2(Category): object(3)

            data_long{s}(10,3) = 1; % IV2(Category): body(1)
            data_long{s}(11,3) = 2; % IV2(Category): face(2)
            data_long{s}(12,3) = 3; % IV2(Category): object(3)

            % label condition levels for IV3
            data_long{s}(1:6,4) = 1; % IV3(Configuration): normal(1)
            data_long{s}(7:12,4) = 2; % IV3(Configuration): scramble(2)

        end

        % create a table of subject-level data in *long format*
        data_allsubj = cat(1, data_long{:}); % concatenate all subject-level data
        variable_names = {'Subject' 'Species','Category','Configuration','Measurement'}; % assign variable names 
        T = array2table(data_allsubj, 'VariableNames', variable_names); % create table 

        % save data as a table (.xlsx) in the 'statistics/R' directory, with its file name and folder name representing the analysis
        % data in *long format* are prepared for statistical analysis in *R*
        analysis = strcat('PLI_allcond_',roi_string,'_Within-PLI_',freq_string); % define present analysis
        folder_path = strcat('statistics/R/',analysis);
        file_name = strcat(analysis,'.xlsx');
        full_file_path = strcat(folder_path,'/',file_name);

        % if the folder does not exist, then create it
        if ~isfolder(folder_path)
            mkdir(folder_path);
        end

        % if the file name already exists, delete it
        if exist(full_file_path, 'file') 
            delete(full_file_path); 
        end
        
        writetable(T,full_file_path); % save



        % transform normalized data from wide to long format for statistical analysis in R
        % wide format is organized as: subjects (rows) x conditions (columns)
        % long format contains one experimental observation per row, organized as follows:
        % data_long{s}(:,subject | IV1(Species: 1.Human 2.Monkey) | IV2(Category: 1.Body 2.Face 3.Object) | Measurement(normalized PLI))

        for s = 1:29
            
            % extract measurements for each subject and condition 
            data_long2{s}(1,4) = PLI_normalized(s,1); % human body condition
            data_long2{s}(2,4) = PLI_normalized(s,2); % human face condition
            data_long2{s}(3,4) = PLI_normalized(s,3); % human object condition

            data_long2{s}(4,4) = PLI_normalized(s,4); % monkey body condition
            data_long2{s}(5,4) = PLI_normalized(s,5); % monkey face condition
            data_long2{s}(6,4) = PLI_normalized(s,6); % monkey object condition
            
            % label subject number
            data_long2{s}(:,1) = s; 

            % label condition levels for IV1
            data_long2{s}(1:3,2) = 1; % IV1(Species): human(1)
            data_long2{s}(4:6,2) = 2; % IV1(Species): monkey(2)

            % label condition levels for IV2
            data_long2{s}(1,3) = 1; % IV2(Category): body(1)
            data_long2{s}(2,3) = 2; % IV2(Category): face(2)
            data_long2{s}(3,3) = 3; % IV2(Category): object(3)

            data_long2{s}(4,3) = 1; % IV2(Category): body(1)
            data_long2{s}(5,3) = 2; % IV2(Category): face(2)
            data_long2{s}(6,3) = 3; % IV2(Category): object(3)

        end

        % create a table of subject-level data in *long format*
        data_allsubj2 = cat(1, data_long2{:}); % concatenate all subject-level data
        variable_names2 = {'Subject' 'Species','Category','Measurement'}; % assign variable names
        T = array2table(data_allsubj2, 'VariableNames', variable_names2); % create table 

        % save data as a table (.xlsx) in the 'statistics/R' directory, with its file name and folder name representing the analysis
        % data in *long format* are prepared for statistical analysis in *R*
        analysis = strcat('PLI_normalized_',roi_string,'_Within-PLI_',freq_string); % define present
        folder_path = strcat('statistics/R/',analysis);
        file_name = strcat(analysis,'.xlsx');
        full_file_path = strcat(folder_path,'/',file_name);

        % if the folder does not exist, then create it
        if ~isfolder(folder_path)
            mkdir(folder_path);
        end

        % if the file name already exists, delete it
        if exist(full_file_path, 'file')
            delete(full_file_path); 
        end
        
        writetable(T,full_file_path); % save



    end

end



%% ------------------------------------------------------------------------
%  Part 4.3 - Real data: Compute between-region PLIs 
%  ------------------------------------------------------------------------
% Between-region PLI refers to the connectivity between two specified clusters 
% In the present analysis, we investigate the connectivity between two clusters and the rest of the scalp 

% pre-allocate variables
PLI_normalized = zeros(total_subjects,6); % subjects (29) x conditions (6)
data_long = cell(1,total_subjects);
data_long2 = cell(1,total_subjects);


% compute between-region PLIs, separately for two ROIs 
for cluster = 1:2

    if cluster == 1

        % connectivity of human body cluster to rest of scalp
        roi_string = 'roiHumanBody';

        % roi1 = 'human body selective cluster'
        roi1 = {'C3', 'CP3', 'P3'};
        roi1_idx = find(ismember(all_channels, roi1));

        % roi2 = rest of scalp
        roi2_idx = find(~ismember(all_channels,roi1));

    elseif cluster == 2

        % connectivity of object-level cluster to rest of scalp
        roi_string = 'roiObjectLvl';

        % roi1 = 'object-level processing cluster'
        roi1 = {'AFz', 'FCz', 'Cz', 'CPz', 'Pz', 'Fp1', 'Fp2', 'F3', 'F4', 'F7', 'F8', 'FC3', 'FC4', 'FT7', 'FT8', 'C3', 'C4', 'CP3', 'CP4', 'TP10', 'P3', 'P4', 'P8'} ;
        roi1_idx = find(ismember(all_channels, roi1));

        % roi2 = rest of scalp
        roi2_idx = find(~ismember(all_channels,roi1));

    end



    % compute PLIs separately for each frequency band
    for freq = 1:length(all_freq)

        freq_string = all_freq{freq};
        
        % compute between-region connectivity for each subject and condition
        % PLI_allconditions = subject (3) x condition (12) array

        PLI_allconditions(:,1) = PLI_between(PLI_hum_body_norm, roi1_idx, roi2_idx,freq);
        PLI_allconditions(:,2) = PLI_between(PLI_hum_face_norm, roi1_idx, roi2_idx,freq);
        PLI_allconditions(:,3) = PLI_between(PLI_hum_obj_norm, roi1_idx, roi2_idx,freq);

        PLI_allconditions(:,4) = PLI_between(PLI_monk_body_norm, roi1_idx, roi2_idx,freq);
        PLI_allconditions(:,5) = PLI_between(PLI_monk_face_norm, roi1_idx, roi2_idx,freq);
        PLI_allconditions(:,6) = PLI_between(PLI_monk_obj_norm, roi1_idx, roi2_idx,freq);

        PLI_allconditions(:,7) = PLI_between(PLI_hum_body_scr, roi1_idx, roi2_idx,freq);
        PLI_allconditions(:,8) = PLI_between(PLI_hum_face_scr, roi1_idx, roi2_idx,freq);
        PLI_allconditions(:,9) = PLI_between(PLI_hum_obj_scr, roi1_idx, roi2_idx,freq);

        PLI_allconditions(:,10) = PLI_between(PLI_monk_body_scr, roi1_idx, roi2_idx,freq);
        PLI_allconditions(:,11) = PLI_between(PLI_monk_face_scr, roi1_idx, roi2_idx,freq);
        PLI_allconditions(:,12) = PLI_between(PLI_monk_obj_scr, roi1_idx, roi2_idx,freq);



        % create a table of subject-level data in *wide format* 
        variable_names = {'hum_body_norm','hum_face_norm','hum_obj_norm','monk_body_norm','monk_face_norm','monk_obj_norm', 'hum_body_scr','hum_face_scr','hum_obj_scr','monk_body_scr','monk_face_scr','monk_obj_scr'};
        T = array2table(PLI_allconditions, 'VariableNames', variable_names);

        % save data as a table (.xlsx) in the 'statistics/SPSS' directory, with its file name and folder name representing the analysis
        % data in *wide format* are prepared for statistical analysis in *SPSS*
        analysis = strcat('PLI_allcond_',roi_string,'_Between-PLI_',freq_string);
        folder_path = strcat('statistics/SPSS/',analysis);
        file_name = strcat(analysis,'.xlsx');
        full_file_path = strcat(folder_path,'/',file_name);

        % if the folder does not exist, then create it
        if ~isfolder(folder_path)
            mkdir(folder_path);
        end
        
        % if the file name already exists, then delete it 
        if exist(full_file_path, 'file') 
            delete(full_file_path); 
        end

        writetable(T,full_file_path); % save



        % normalize PLIs
        for c = 1:6
            PLI_normalized(:,c) = PLI_allconditions(:,c)-PLI_allconditions(:,c+6);
        end



        % create a table of subject-level data in *wide format* 
        variable_names = {'hum_body','hum_face','hum_obj','monk_body','monk_face','monk_obj'};
        T = array2table(PLI_normalized, 'VariableNames', variable_names);

        % save data as a table (.xlsx) in the 'statistics/SPSS' directory, with its file name and folder name representing the analysis
        % data in *wide format* are prepared for statistical analysis in *SPSS*
        analysis = strcat('PLI_normalized_',roi_string,'_Between-PLI_',freq_string);
        folder_path = strcat('statistics/SPSS/',analysis);
        file_name = strcat(analysis,'.xlsx');
        full_file_path = strcat(folder_path,'/',file_name);

        % if the folder does not exist, then create it
        if ~isfolder(folder_path)
            mkdir(folder_path);
        end

        % if the file name already exists, delete it 
        if exist(full_file_path, 'file') 
            delete(full_file_path); 
        end

        writetable(T,full_file_path); % save 



        % transform non-normalized data from wide to long format for statistical analysis in R

        for s = 1:29
            
            % extract measurements for each subject and condition 
            data_long{s}(1,5) = PLI_allconditions(s,1); % human body normal
            data_long{s}(2,5) = PLI_allconditions(s,2); % human face normal 
            data_long{s}(3,5) = PLI_allconditions(s,3); % human object normal

            data_long{s}(4,5) = PLI_allconditions(s,4); % monkey body normal 
            data_long{s}(5,5) = PLI_allconditions(s,5); % monkey face normal 
            data_long{s}(6,5) = PLI_allconditions(s,6); % monkey object normal 

            data_long{s}(7,5) = PLI_allconditions(s,7); % human body scramble
            data_long{s}(8,5) = PLI_allconditions(s,8); % human face scramble 
            data_long{s}(9,5) = PLI_allconditions(s,9); % human object scramble

            data_long{s}(10,5) = PLI_allconditions(s,10); % monkey body scramble 
            data_long{s}(11,5) = PLI_allconditions(s,11); % monkey face scramble 
            data_long{s}(12,5) = PLI_allconditions(s,12); % monkey object scramble 


            % label subject number
            data_long{s}(:,1) = s; % subject number
            
            % label condition levels for IV1
            data_long{s}(1:3,2) = 1; % IV1(Species): human(1)
            data_long{s}(7:9,2) = 1; % IV1(Species): human(1)

            data_long{s}(4:6,2) = 2; % IV1(Species): monkey(2)
            data_long{s}(10:12,2) = 2; % IV1(Species): monkey(2)
            
            % label condition levels for IV2
            data_long{s}(1,3) = 1; % IV2(Category): body(1)
            data_long{s}(2,3) = 2; % IV2(Category): face(2)
            data_long{s}(3,3) = 3; % IV2(Category): object(3)

            data_long{s}(4,3) = 1; % IV2(Category): body(1)
            data_long{s}(5,3) = 2; % IV2(Category): face(2)
            data_long{s}(6,3) = 3; % IV2(Category): object(3)

            data_long{s}(7,3) = 1; % IV2(Category): body(1)
            data_long{s}(8,3) = 2; % IV2(Category): face(2)
            data_long{s}(9,3) = 3; % IV2(Category): object(3)

            data_long{s}(10,3) = 1; % IV2(Category): body(1)
            data_long{s}(11,3) = 2; % IV2(Category): face(2)
            data_long{s}(12,3) = 3; % IV2(Category): object(3)

            % label condition levels for IV3
            data_long{s}(1:6,4) = 1; % IV3(Configuration): normal(1)
            data_long{s}(7:12,4) = 2; % IV3(Configuration): scramble(2)

        end

        % create a table of subject-level data in *long format*
        data_allsubj = cat(1, data_long{:}); % concatenate all subject-level data
        variable_names = {'Subject' 'Species','Category','Configuration','Measurement'}; % assign variable names
        T = array2table(data_allsubj, 'VariableNames', variable_names); % create table

        % save data as a table (.xlsx) in the 'statistics/R' directory, with its file name and folder name representing the analysis
        % data in *long format* are prepared for statistical analysis in *R*        
        analysis = strcat('PLI_allcond_',roi_string,'_Between-PLI_',freq_string); % define present analysis
        folder_path = strcat('statistics/R/',analysis);
        file_name = strcat(analysis,'.xlsx');
        full_file_path = strcat(folder_path,'/',file_name);

        % if the folder does not exist, then create it
        if ~isfolder(folder_path)
            mkdir(folder_path);
        end

        % if the file name already exists, delete it
        if exist(full_file_path, 'file') 
            delete(full_file_path); 
        end

        writetable(T,full_file_path); % save new file



        % transform normalized data from wide to long format for statistical analysis in R
        for s = 1:29

            % extract measurements for each subject and condition 
            data_long2{s}(1,4) = PLI_normalized(s,1); % human body condition
            data_long2{s}(2,4) = PLI_normalized(s,2); % human face condition
            data_long2{s}(3,4) = PLI_normalized(s,3); % human object condition

            data_long2{s}(4,4) = PLI_normalized(s,4); % monkey body condition
            data_long2{s}(5,4) = PLI_normalized(s,5); % monkey face condition
            data_long2{s}(6,4) = PLI_normalized(s,6); % monkey object condition
        
            % label subject number
            data_long2{s}(:,1) = s; 
            
            % label condition levels for IV1
            data_long2{s}(1:3,2) = 1; % IV1(Species): human(1)
            data_long2{s}(4:6,2) = 2; % IV1(Species): monkey(2)

            % label condition levels for IV2
            data_long2{s}(1,3) = 1; % IV2(Category): body(1)
            data_long2{s}(2,3) = 2; % IV2(Category): face(2)
            data_long2{s}(3,3) = 3; % IV2(Category): object(3)

            data_long2{s}(4,3) = 1; % IV2(Category): body(1)
            data_long2{s}(5,3) = 2; % IV2(Category): face(2)
            data_long2{s}(6,3) = 3; % IV2(Category): object(3)

        end

        % create a table of subject-level data in *long format*
        data_allsubj2 = cat(1, data_long2{:}); % concatenate all subject-level data
        variable_names2 = {'Subject' 'Species','Category','Measurement'}; % assign variable names
        T = array2table(data_allsubj2, 'VariableNames', variable_names2); % create table 

        % save data as a table (.xlsx) in the 'statistics/R' directory, with its file name and folder name representing the analysis
        % data in *long format* are prepared for statistical analysis in *R*
         analysis = strcat('PLI_normalized_',roi_string,'_Between-PLI_',freq_string);
        folder_path = strcat('statistics/R/',analysis);
        file_name = strcat(analysis,'.xlsx');
        full_file_path = strcat(folder_path,'/',file_name);

        % if the folder does not exist, then create it
        if ~isfolder(folder_path)
            mkdir(folder_path);
        end

        % if the file name already exists, delete it
        if exist(full_file_path, 'file') 
            delete(full_file_path); 
        end

        writetable(T,full_file_path); % save 


    end

end



%% ------------------------------------------------------------------------
%  Part 4.4 - PLI data-driven ROI selection 
%  ------------------------------------------------------------------------
% Select ROI based on statistical comparison of PLI data of normal versus scramble conditions
% Compute this on a new time-window of 500-1000ms post-stimulus 

%% Part 4.3.1 Compute PLIs on new time window 

% define time-window of interest (TOI) in seconds
min_T1 = 0.5; % TOI minimum 
max_T1 = 1; % TOI maximum 

% load input data
load(fullfile('output','data_clean_pooled_normal.mat'))
load(fullfile('output','data_clean_pooled_scramble.mat'))
load(fullfile('output','all_channels.mat'))
all_freq = {'delta';'theta';'alpha';'beta';'gamma_A_35_45';'gamma_B_25_45_BS30';'broadband_1_45'};

PLI_pooled_normal = PLI_all_subjects(data_clean_pooled_normal, min_T1, max_T1);
save('output/PLI_pooled_normal_500_1000ms.mat','PLI_pooled_normal')

PLI_pooled_scramble = PLI_all_subjects(data_clean_pooled_scramble, min_T1, max_T1);
save('output/PLI_pooled_scramble_500_1000ms.mat','PLI_pooled_scramble')


%% Statistically compare Normal vs Scramble conditions 

% load input data
% input data are 4-D PLI data (freq x channel x channel x subject) for normal and scramble conditions, computed above 
load(fullfile('output','PLI_pooled_normal_500_1000ms.mat'));
load(fullfile('output','PLI_pooled_scramble_500_1000ms.mat'));

% pre-allocate variables for speed
NS_500_1000ms = cell(1,length(all_freq));

% begin loop for all frequencies 
for freq = 1:length(all_freq)
    
    % pre-allocate variables for speed 
    normal_condition = cell(length(all_channels),length(all_channels));
    scramble_condition = cell(length(all_channels),length(all_channels));
    results_normality = cell(length(all_channels),length(all_channels));
    results_ttest = zeros(length(all_channels));
    results_wilrank = zeros(length(all_channels));

    % restructure data for normal and scramble conditions
    % E.g.  original data = PLI_pooled_normal{freq}(channel,channel,subject)
    %       restructured data = normal_condition{channel,channel}(subject)
    for i = 1:length(all_channels) % loop for all channel pairs 
        for j = 1:length(all_channels)
            % create a channel x channel cell with all subject data in each cell,
            % separately for each condition 
            normal_condition{i,j} = squeeze(PLI_pooled_normal{freq}(i,j,:)); 
            scramble_condition{i,j} = squeeze(PLI_pooled_scramble{freq}(i,j,:));
        end
    end

    % statistically compare normal and scramble conditions
    for i = 1:length(all_channels)

        for j = 1:length(all_channels)


            % For each channel pair, calculate paired differences  
            % We are interested in paired differences because this
            % experiment has a repeated measures design (each subject undergoes all conditions)
            % Thus, we check normality and perform statistical tests on the paired differences 
            differences = normal_condition{i,j} - scramble_condition{i,j};


            % Check normality assumption: Perform Shapiro-Wilk test on paired differences 
            [H, pValue, W] = swtest(differences);
            % Store results in a channel x channel matrix 
            if H == 0
                results_normality{i,j} = 'normal';
            else
                results_normality{i,j} = 'ATTN: normality assumption violated';
            end


            % Perform a one-tailed paired samples t-test 
            % *but only interpret this if the assumption of normality is met*
            % 'Tail', 'right' = right-tailed, x > y
            [h, p, ci, stats] = ttest(normal_condition{i,j}, scramble_condition{i,j}, 'Tail', 'right');
            % Store results in a channel x channel matrix 
            results_ttest(i,j) = p; 


            % Perform non-parametric Wilcoxon-rank test on the paired differences
            % *but only intepret if the assumption of normality is violated*
            [p2, h, stats] = signrank(differences);
            % Store results in a channel x channel matrix 
            results_wilrank(i,j) = p2;

        end
    end
    
    % Extract p-values for all unique channel pairs, excluding same-channel pairs 
    dummy = ones(33);
    tri_idx = ~tril(dummy); % get indices of the upper triangle, excluding the diagonal 
    pairs_unique = results_wilrank(tri_idx);

    % Perform FDR-correction for multiple comparisons 
    raw_pValues = pairs_unique; 
    q = 0.05;
    [h, crit_p, adj_ci_cvrg, adj_p] = fdr_bh(raw_pValues,q,'pdep','yes');
    adj_p = round(adj_p,4);

    % Store results in a channel x channel matrix
    results_wilrank_adj_p = zeros(33); % create a channel x channel matrix of zeros 
    results_wilrank_adj_p(tri_idx) = adj_p; % assign WR adj p values to upper triangle, excluding diagonal
    results_wilrank_adj_p = results_wilrank_adj_p + results_wilrank_adj_p'; % make the matrix symmetrical wrt the diagonal
    %     % Check if the matrix is symmetric
    %     isequal(results_wilrank_adj_p, results_wilrank_adj_p');

    % Consolidate results
    NS_500_1000ms{freq}.frequency = all_freq{freq};
    NS_500_1000ms{freq}.results_normality = results_normality;
    NS_500_1000ms{freq}.results_ttest = results_ttest;
    NS_500_1000ms{freq}.results_wilrank = results_wilrank;
    NS_500_1000ms{freq}.results_wilrank_adj_p = results_wilrank_adj_p;

end

% Save results as .mat 
save('output/NS_500_1000ms.mat','NS_500_1000ms')

% Save connectivity matrices as tables for visualizations in R
for i = 1:length(all_freq)
    
    % 1. Create a table of Channel x Channel Wilcoxon Rank p-values (raw) 
    T = array2table(NS_500_1000ms{i}.results_wilrank);
    full_file_path = strcat('output/','connectivity_matrix_WR_p_',all_freq{i},'.xlsx');
    % if the file name already exists, delete it
    if exist(full_file_path, 'file')
        delete(full_file_path);
    end
    writetable(T,full_file_path,"WriteVariableNames",0);

    % 2. Create a table of Channel x Channel Wilcoxon Rank adjusted p-values (FDR) 
    T2 = array2table(NS_500_1000ms{i}.results_wilrank_adj_p);
    full_file_path = strcat('output/','connectivity_matrix_WR_adj_p_',all_freq{i},'.xlsx');
    % if the file name already exists, delete it
    if exist(full_file_path, 'file')
        delete(full_file_path);
    end
    writetable(T2,full_file_path,"WriteVariableNames",0);

end















%% Notes 

% 
% saving all files with a specific name
% %  Save %%%%%% ATTN - do this individually in lines above (running PLIs takes super long)
% variables = who('*PLI*');  % Get all variables containing 'PLI' in their name
% for i = 1:length(variables)
%     save(fullfile('output', [variables{i} '.mat']), variables{i}); % save all to /output
% end

% % Create a 33x33 matrix
% matrix = zeros(5)
% 
% % Define the values you want to assign to the upper or lower triangle
% values = [1:10]; % Example array of values
% matrix(triu(true(5), 1)) = values; % assign values to upper triangle
% matrix(tril(true(5), -1)) = values; % assign values to lower triangle 
% 
% % Display the modified matrix
% disp(matrix)



% 
% % confirm data restructuring was done correctly
% n = 1;
% for channel1 = 1:33
%     for channel2 = 1:33
%         for subject = 1:29
%             x = PLI_pooled_normal{freq}(channel1,channel2,subject);
%             y = normal_condition{channel1,channel2}(subject);
%             list(n) = ~isequal(x,y);
%             n = n+1;
%         end
%     end
% end
% nnz(list)
% 


