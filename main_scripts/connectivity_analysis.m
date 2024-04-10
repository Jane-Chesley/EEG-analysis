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
%   2. PLI computation
%   3. Extract PLIs of interest (within-region and between-region) 
%   
% ------------------------------------------------------------------------



tic; 

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
% dir_FT = '';
dir_FT = '/Users/jane_chesley/Library/Application Support/MathWorks/MATLAB Add-Ons/Collections/FieldTrip';

% initialize Fieldtrip toolbox
addpath(dir_FT); ft_defaults;



%% ------------------------------------------------------------------------
%  Part 2 - PLI computation for all subjects, frequencies, and channels
%  ------------------------------------------------------------------------

% identify all input data to analyze
% data to analyze are preprocessed EEG data, stratified by condition
input = dir(fullfile('output', '*data_clean*'));

% load all (this is slow because files are large)
disp('loading data...')
for i = 1: length(input)
    disp(strcat('loading data...', num2str(i),' of ', num2str(length(input))));
    load(strcat(input(i).folder, '/',input(i).name));
end

% define toi range in SECONDS
min_T1 = 0;
max_T1 = 1;

% calculate subject-level PLIs for each frequency band, averaged across trials
% input is subject-level clean data for ONE condition
% output for each condition includes freq x channel x channel x subject data;
% organized as follows: PLI_hum_body_norm{freq}(channel,channel,subject);
% PLI values are represented in cells of each (channel,channel) matrix; (0 = no connectivity; 1 = perfect connectivity)
% frequency bands are segmented as follows: {freq1: delta 0.5-3.9 Hz} {freq2: theta 4-7 Hz} {freq3: alpha 8-12 Hz}
% {freq4: beta 13-30 Hz} {freq5: gamma_A 35-45 Hz} {freq6: gamma_B 25-45 Hz + bandstop 29:31 Hz}
PLI_hum_body_norm = PLI_all_subjects(data_clean_hum_body_norm, min_T1, max_T1);
PLI_hum_face_norm = PLI_all_subjects(data_clean_hum_face_norm, min_T1, max_T1);
PLI_hum_obj_norm = PLI_all_subjects(data_clean_hum_obj_norm, min_T1, max_T1);

PLI_monk_body_norm = PLI_all_subjects(data_clean_monk_body_norm, min_T1, max_T1);
PLI_monk_face_norm = PLI_all_subjects(data_clean_monk_face_norm, min_T1, max_T1);
PLI_monk_obj_norm = PLI_all_subjects(data_clean_monk_obj_norm, min_T1, max_T1);

PLI_hum_body_scr = PLI_all_subjects(data_clean_hum_body_scr, min_T1, max_T1);
PLI_hum_face_scr = PLI_all_subjects(data_clean_hum_face_scr, min_T1, max_T1);
PLI_hum_obj_scr = PLI_all_subjects(data_clean_hum_obj_scr, min_T1, max_T1);

PLI_monk_body_scr = PLI_all_subjects(data_clean_monk_body_scr, min_T1, max_T1);
PLI_monk_face_scr = PLI_all_subjects(data_clean_monk_face_scr, min_T1, max_T1);
PLI_monk_obj_scr = PLI_all_subjects(data_clean_monk_obj_scr, min_T1, max_T1);

PLI_pooled_normal = PLI_all_subjects(data_clean_pooled_normal, min_T1, max_T1);
PLI_pooled_scramble = PLI_all_subjects(data_clean_pooled_scramble, min_T1, max_T1);



%  Save
variables = who('*PLI*');  % Get all variables containing 'PLI' in their name
for i = 1:length(variables)
    save(fullfile('output', [variables{i} '.mat']), variables{i}); % save all to /output
end



%% ------------------------------------------------------------------------
%  Part 3 - Extract relevant PLIs 
%  ------------------------------------------------------------------------

% identify all input data to analyze
% data to analyze are PLIs for each subject*condition*frequency (computed above in Step 2)
input = dir(fullfile('output','*PLI*'));

% load all
for i = 1: length(input)
    load(strcat(input(i).folder, '/',input(i).name));
end

% define all frequency bands
all_freq = {'delta';'theta';'alpha';'beta';'gamma_A';'gamma_B'};

% get all channel labels for roi extraction 
% channel labels are constant across all measurements
all_channels = data_clean_hum_body_norm{1}.label; 



%% ------------------------------------------------------------------------
%  Within-region PLI
%  ------------------------------------------------------------------------

% pre-allocate vars
PLI_normalized = zeros(29,6);
data_long = cell(1,29);

% compute within-region PLIs, separately for two pre-defined clusters
for cluster = 1:2

    if cluster == 1
        roi = {'C3', 'CP3', 'P3'};
        roi_string = 'roiHumanBody'; % 'human body selective cluster'

    elseif cluster == 2
        roi = {'AFz', 'FCz', 'Cz', 'CPz', 'Pz', 'Fp1', 'Fp2', 'F3', 'F4', 'F7', 'F8', 'FC3', 'FC4', 'FT7', 'FT8', 'C3', 'C4', 'CP3', 'CP4', 'TP10', 'P3', 'P4', 'P8'} ;
        roi_string = 'roiObjectLvl'; % 'object-level processing cluster'
    end

    % indices of roi
    roi_idx = find(ismember(all_channels, roi));



    % compute separately for each frequency band
    for freq = 1:length(all_freq)

        freq_string = all_freq{freq};

        % extract relevant information:
        % Within-region connectivity: subject- and condition-level
        % PLI_allconditions(subject,condition)

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



        % compute normalization (normal - scramble)
        % PLI_normalized(subject,condition)
        % conditions = 'hum_body','hum_face','hum_obj','monk_body','monk_face','monk_obj';
        for c = 1:6
            PLI_normalized(:,c) = PLI_allconditions(:,c)-PLI_allconditions(:,c+6);
        end



        % save data in wide format (for statistical analysis in SPSS)
        variable_names = {'hum_body','hum_face','hum_obj','monk_body','monk_face','monk_obj'};
        T = array2table(PLI_normalized, 'VariableNames', variable_names);
        
        analysis = strcat(roi_string,'_Within-PLI_',freq_string);
        folder_path = strcat('statistics/SPSS/',analysis); 
        file_name = strcat(analysis,'.xlsx');
        full_file_path = fullfile(folder_path, file_name);  % combine folder path and file name

        if exist(full_file_path, 'file') % check if the file name already exists
            delete(full_file_path); % if it exists, delete it
        end

        writetable(T,full_file_path); % save new file



        % transform data from wide format to long format (for statistical analysis in R)
        % format is as follows:
        % data_long{s}(:,subject | IV1(Species: 1.Human 2.Monkey) | IV2(Category: 1.Body 2.Face 3.Object) | Measurement(normalized PLI))

        for s = 1:29

            data_long{s}(1,4) = PLI_normalized(s,1); % human body
            data_long{s}(2,4) = PLI_normalized(s,2); % human face
            data_long{s}(3,4) = PLI_normalized(s,3); % human object

            data_long{s}(4,4) = PLI_normalized(s,4); % monkey body
            data_long{s}(5,4) = PLI_normalized(s,5); % monkey face
            data_long{s}(6,4) = PLI_normalized(s,6); % monkey object

            data_long{s}(:,1) = s; % subject number

            data_long{s}(1:3,2) = 1; % human
            data_long{s}(4:6,2) = 2; % monkey

            data_long{s}(1,3) = 1; % body
            data_long{s}(2,3) = 2; % face
            data_long{s}(3,3) = 3; % object

            data_long{s}(4,3) = 1; % body
            data_long{s}(5,3) = 2; % face
            data_long{s}(6,3) = 3; % object

        end

        % concatenate all subject-level data
        % data_allsubj(:,subject | IV1(Species: 1.Human 2.Monkey) | IV2(Category: 1.Body 2.Face 3.Object) | Measurement(normalized PLI))
        data_allsubj = cat(1, data_long{:});
        variable_names = {'Subject' 'Species','Category','Measurement'};

        T = array2table(data_allsubj, 'VariableNames', variable_names);

        % save
        analysis = strcat(roi_string,'_Within-PLI_',freq_string);
        folder_path = strcat('statistics/R/',analysis); 
        file_name = strcat(analysis,'.xlsx');
        full_file_path = fullfile(folder_path, file_name);  % combine folder path and file name

        if exist(full_file_path, 'file') % check if the file name already exists
            delete(full_file_path); % if it exists, delete it
        end

        writetable(T,full_file_path); % save new file



    end

end




%% ------------------------------------------------------------------------
%  Between-region PLI
%  ------------------------------------------------------------------------

% pre-allocate vars 
PLI_normalized = zeros(29,6);
data_long = cell(1,29);

% compute between-region PLIs, separately for two pre-defined clusters
for cluster = 1:2

    if cluster == 1

        % connectivity of human body cluster to rest of scalp
        roi_string = 'roiHumanBody';

        % roi1 = 'human body selective cluster'
        roi1 = {'C3', 'CP3', 'P3'};
        roi1_idx = find(ismember(all_channels, roi1));

        % roi2 = rest of scalp
        roi2 = {'AFz','Fz','FCz','Cz','CPz','Pz','Oz','Fp1','Fp2','F3','F4','F7','F8','FC3','FC4','FT7','FT8','C4','T7','T8','CP4','TP7','TP8','TP9','TP10','P4','P7','P8','O1','O2'};
        roi2_idx = find(ismember(all_channels, roi2));

    elseif cluster == 2

        % connectivity of object-level cluster to rest of scalp
        roi_string = 'roiObjectLvl';

        % roi1 = 'object-level processing cluster'
        roi1 = {'AFz', 'FCz', 'Cz', 'CPz', 'Pz', 'Fp1', 'Fp2', 'F3', 'F4', 'F7', 'F8', 'FC3', 'FC4', 'FT7', 'FT8', 'C3', 'C4', 'CP3', 'CP4', 'TP10', 'P3', 'P4', 'P8'} ;
        roi1_idx = find(ismember(all_channels, roi1));

        % roi2 = rest of scalp
        roi2 =     {'Fz','O1','O2','Oz','P7','T7','T8','TP7','TP8','TP9'};
        roi2_idx = find(ismember(all_channels, roi2));


    end



    % compute separately for each frequency band
    for freq = 1:6

        freq_string = all_freq{freq};

        % extract relevant information:
        % Between-region connectivity: subject- and condition-level
        % PLI_allconditions(subject,condition)

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



        % compute normalization (normal - scramble)
        % PLI_normalized(subject,condition)
        % conditions = 'hum_body','hum_face','hum_obj','monk_body','monk_face','monk_obj';
        for c = 1:6
            PLI_normalized(:,c) = PLI_allconditions(:,c)-PLI_allconditions(:,c+6);
        end



        % save data in wide format (for statistical analysis in SPSS)
        variable_names = {'hum_body','hum_face','hum_obj','monk_body','monk_face','monk_obj'};
        T = array2table(PLI_normalized, 'VariableNames', variable_names);

        analysis = strcat(roi_string,'_Between-PLI_',freq_string);
        folder_path = strcat('statistics/SPSS/',analysis); 
        file_name = strcat(analysis,'.xlsx');
        full_file_path = fullfile(folder_path, file_name);  % combine folder path and file name




        if exist(full_file_path, 'file') % check if the file name already exists
            delete(full_file_path); % if it exists, delete it
        end

        writetable(T,full_file_path); % save new file



        % transform data from wide format to long format (for statistical analysis in R)
        % data_long{s}(:,subject | IV1(Species: 1.Human 2.Monkey) | IV2(Category: 1.Body 2.Face 3.Object) | Measurement(normalized PLI))

        for s = 1:29

            data_long{s}(1,4) = PLI_normalized(s,1); % human body
            data_long{s}(2,4) = PLI_normalized(s,2); % human face
            data_long{s}(3,4) = PLI_normalized(s,3); % human object

            data_long{s}(4,4) = PLI_normalized(s,4); % monkey body
            data_long{s}(5,4) = PLI_normalized(s,5); % monkey face
            data_long{s}(6,4) = PLI_normalized(s,6); % monkey object

            data_long{s}(:,1) = s; % subject number

            data_long{s}(1:3,2) = 1; % human
            data_long{s}(4:6,2) = 2; % monkey

            data_long{s}(1,3) = 1; % body
            data_long{s}(2,3) = 2; % face
            data_long{s}(3,3) = 3; % object

            data_long{s}(4,3) = 1; % body
            data_long{s}(5,3) = 2; % face
            data_long{s}(6,3) = 3; % object

        end

        % Concatenate all subject-level data
        % data_allsubj(:,subject | IV1(Species: 1.Human 2.Monkey) | IV2(Category: 1.Body 2.Face 3.Object) | Measurement(normalized PLI))
        data_allsubj = cat(1, data_long{:});
        variable_names = {'Subject' 'Species','Category','Measurement'};

        T = array2table(data_allsubj, 'VariableNames', variable_names);

        % save
        analysis = strcat(roi_string,'_Between-PLI_',freq_string);
        folder_path = strcat('statistics/R/',analysis); 
        file_name = strcat(analysis,'.xlsx');
        full_file_path = fullfile(folder_path, file_name);  % combine folder path and file name


        if exist(full_file_path, 'file') % check if the file name already exists
            delete(full_file_path); % if it exists, delete it
        end

        writetable(T,full_file_path); % save new file


    end

end

toc; 

