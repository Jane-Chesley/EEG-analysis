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
%   - Derivative data are saved to '/EEG-analysis/derivatives'
%   - Output is saved to '/EEG-analysis/output'
%
% Author:
%   Jane Chesley
%
% Affiliation:
%   Faculty of Psychology & Neuroscience, Maastricht University
%   Maastricht, the Netherlands
%
% Date:
%   8 November, 2023
%
% Version:
%   1.0
%
% Outline of script content:
%   1. Script setup
%   2. Analysis
%   3. Save 

% ------------------------------------------------------------------------



%% ------------------------------------------------------------------------
%  Part 1 - Script setup
%  ------------------------------------------------------------------------

% Part 1.1 clean working environment --------------------------------------------------

clear, clc, close all;
restoredefaultpath;

% Part 1.2 set up relevant directories --------------------------------------------------
[dir_parent, ~, ~] = fileparts(pwd);
cd(dir_parent);

% Part 1.3 initialize functions --------------------------------------------------

% add custom functions to path
addpath('functions');

% specify the path to the Fieldtrip Toolbox on your system
% dir_FT = '';
dir_FT = '/Users/jane_chesley/Library/Application Support/MathWorks/MATLAB Add-Ons/Collections/FieldTrip';

% initialize the toolbox
addpath(dir_FT); ft_defaults;



%% ------------------------------------------------------------------------
%  Part 2 - Analysis
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
% frequency bands are segmented as follows: {freq1: delta 1-3 Hz} {freq2: theta 4-7 Hz} {freq3: alpha 8-12 Hz} {freq4: beta 13-30 Hz} 
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



%% ------------------------------------------------------------------------
%  Part 3 - Save  
%  ------------------------------------------------------------------------

variables = who('*PLI*');  % Get all variables containing 'PLI' in their name
for i = 1:length(variables)
    save(fullfile('output', [variables{i} '.mat']), variables{i}); % save all to /output
end


%% ------------------------------------------------------------------------
%  Part 4 - Extract average PLI within subregion and frequency band of interest  
%  ------------------------------------------------------------------------

% identify all input data to analyze
% data to analyze are preprocessed EEG data, stratified by condition
input = dir(fullfile('output','*PLI*'));

% load all 
for i = 1: length(input)
    load(strcat(input(i).folder, '/',input(i).name));
end


all_channels = data_clean_hum_body_norm{1}.label; % channel labels are constant across all measurements 
roi = {'C3', 'CP3', 'P3'};

roi_idx = find(ismember(all_channels, roi)); 
freq = 2; % theta 


PLI_allconditions(:,1) = PLI_extraction(PLI_hum_body_norm, roi_idx, freq);
PLI_allconditions(:,2) = PLI_extraction(PLI_hum_face_norm, roi_idx, freq);
PLI_allconditions(:,3) = PLI_extraction(PLI_hum_obj_norm, roi_idx, freq);

PLI_allconditions(:,4) = PLI_extraction(PLI_monk_body_norm, roi_idx, freq);
PLI_allconditions(:,5) = PLI_extraction(PLI_monk_face_norm, roi_idx, freq);
PLI_allconditions(:,6) = PLI_extraction(PLI_monk_obj_norm, roi_idx, freq);

PLI_allconditions(:,7) = PLI_extraction(PLI_hum_body_scr, roi_idx, freq);
PLI_allconditions(:,8) = PLI_extraction(PLI_hum_face_scr, roi_idx, freq);
PLI_allconditions(:,9) = PLI_extraction(PLI_hum_obj_scr, roi_idx, freq);

PLI_allconditions(:,10) = PLI_extraction(PLI_monk_body_scr, roi_idx, freq);
PLI_allconditions(:,11) = PLI_extraction(PLI_monk_face_scr, roi_idx, freq);
PLI_allconditions(:,12) = PLI_extraction(PLI_monk_obj_scr, roi_idx, freq);



% Convert the matrix into a table with variable names
variable_names = {'hum_body_norm','hum_face_norm','hum_obj_norm','monk_body_norm','monk_face_norm','monk_obj_norm','hum_body_scr','hum_face_scr','hum_obj_scr','monk_body_scr','monk_face_scr','monk_obj_scr'};

T = array2table(PLI_allconditions, 'VariableNames', variable_names);


% Write the table to an Excel file
writetable(T, 'Stats_PLI_allconditions_check.xlsx');

%% compute normalization (normal - scramble)

for c = 1:6
    PLI_normalized(:,c) = PLI_allconditions(:,c)-PLI_allconditions(:,c+6);
end 
variable_names2 = {'hum_body','hum_face','hum_obj','monk_body','monk_face','monk_obj'};

T2 = array2table(PLI_normalized, 'VariableNames', variable_names2);

% Write the table to an Excel file
writetable(T2,'Stats_PLI_normalized.xlsx');


% %%
% % plotting
% PLI_grouplvl = mean(PLI_avg_allsubjects,3);
% 
% % Plot the channel*channel PLI matrix
% figure;
% imagesc(PLI_grouplvl);
% 
% % Customize the plot
% colorbar;
% caxis([0 1]);
% xticks(1:33)
% yticks(1:33)
% 
% axis square;  % Ensure that the aspect ratio is square
% title('PLI Matrix');
% 
% % Add axis labels
% xlabel('Channels');
% ylabel('Channels');






