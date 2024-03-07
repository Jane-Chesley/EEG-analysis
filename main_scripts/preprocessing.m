%% ------------------------------------------------------------------------
%  EEG Signal Preprocessing
%  ------------------------------------------------------------------------

% Script Description:
% This MATLAB script preprocesses Electroencephalogram (EEG) data for 29 subjects.
% It includes segmentation, extraction, cleaning, quality control and stratification of EEG data,
% as well as event labeling, signal re-referencing, and sampling reduction.
% After executing this script, the resulting data is ready for additional
% analyses such as event-related potential (ERP) analysis and spectral analysis.


% Usage:
%   - This script requires the following MATLAB toolboxes:
%       1. Image Processing Toolbox
%       2. Statistics and Machine Learning Toolbox
%       3. Bioinformatics Toolbox
%       4. FieldTrip Toolbox
%   - Note: check which toolboxes are required by [1] running script [2] executing: [fList,pList] = matlab.codetools.requiredFilesAndProducts('run_analysis.m'); disp({pList.Name}');


% Inputs:
%   - All input data is stored in '/EEG-ERP/input' and includes:
%       1. The subdirectory '/EEG_data', which contains the unprocessed EEG data for 30 subjects (.eeg, .vhdr and .vmrk files)
%       2. The subdirectory '/event_codes', which contains the condition codes (1-12) for each trial and for each subject

% Outputs:
%   - Derivative data are saved to '/EEG-ERP/derivatives'
%   - Output is saved to '/EEG-ERP/output'

% Author:
%   Jane Chesley

% Affiliation:
%   Faculty of Psychology & Neuroscience, Maastricht University
%   Maastricht, the Netherlands

% Date:
%   8 November, 2023

% Version:
%   1.0

% Outline of script content:
%   1. Script setup
%       1.1 clean working environment
%       1.2 set up relevant directories
%       1.3 initialize functions
%   2. Prepare data
%   3. EOG artifact removal
%       3.1 subject-level analysis
%       3.2 group-level analysis
%   4. Remove outlier trials
%       4.1 subject-level analysis
%       4.2 group-level analysis
%   5. Stratification
%       5.1 12-conditions
%       5.2 2-conditions

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
%  Part 2 - prepare data
%  ------------------------------------------------------------------------

% preprocessing steps: trial segmentation and extraction, event labeling, signal re-referencing, sampling reduction



% identify all input data to analyze
input = dir(fullfile(strcat('input','/EEG_data'), '*.eeg'));

for i = 1:length(input)

    % get .eeg file i
    file_EEG = strcat(input(i).folder, '/',input(i).name);

    % load corresponding event codes
    [~, subjectid, ~] = fileparts(input(i).name); % access subject id from current .eeg file
    load(strcat('input','/event_codes/',subjectid,'.mat')); % loads var 'event_codes'

    % preprocess file i
    data_preprocessed       = prepare_data(file_EEG, event_codes);

    % save derivative data i 'data_preprocessed'
    deriv                   = '/deriv01_data_preprocessed/';
    full_file_name          = strcat('derivatives',deriv,subjectid,'.mat');
    disp(strcat('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> saving_', deriv, subjectid,' <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<'));
    save(full_file_name, 'data_preprocessed')

end

clearvars -except dir_parent



%% ------------------------------------------------------------------------
%  Part 3 - EOG artifact removal
%  ------------------------------------------------------------------------

% remove EOG artifacts with ICA



% Part 3.1 subject-level --------------------------------------------------

% identify all data to analyze
% data to analyze is subject-level preprocessed data, constructed above in Part 2
input = dir(fullfile(strcat('derivatives','/deriv01_data_preprocessed'), '*.mat'));

for i = 1:length(input)

    % load input file i
    file_preprocessed = strcat(input(i).folder, '/',input(i).name);
    load(file_preprocessed)

    % run ICA on input file i
    [data_ICA, rejected_artifacts] = artifact_removal(data_preprocessed);

    % extract subject id
    [~, subjectid, ~] = fileparts(input(i).name); % access subject id from current input file

    % save derivative data [1] 'data_preprocessed'
    deriv                   = '/deriv02_data_ICA/'; % update accordingly
    full_file_name          = strcat('derivatives',deriv,subjectid,'.mat');
    disp(strcat('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> saving_', deriv, subjectid,' <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<'));
    save(full_file_name, 'data_ICA')

    % save derivative data [2] 'rejected_artifacts'
    deriv                   = '/deriv03_rejected_artifacts/'; % update accordingly
    full_file_name          = strcat('derivatives',deriv,subjectid,'.mat');
    disp(strcat('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> saving_', deriv, subjectid,' <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<'));
    save(full_file_name, 'rejected_artifacts')

end

clearvars -except dir_parent



% Part 3.2 group-level ----------------------------------------------------

% identify all data to analyze
% data to analyze is subject-level rejected artifacts, constructed above in Part 3.1
input = dir(fullfile(strcat('derivatives','/deriv03_rejected_artifacts'), '*.mat'));

% pre-allocate variables for speed
rejected_artifacts_all_subj = cell(1,length(input)); 
total_artifacts             = zeros(1,29);

for i = 1:length(input)

    % load input file i 'rejected_artifacts'
    file_artifacts = strcat(input(i).folder, '/',input(i).name);
    load(file_artifacts)

    % consolidate all subject data
    rejected_artifacts_all_subj{i}  = rejected_artifacts;
    total_artifacts(i) = length(rejected_artifacts);

end

% consolidate relevant results into one structure
results_artifact_removal.rejected_artifacts     = rejected_artifacts_all_subj;
results_artifact_removal.total_artifacts        = total_artifacts;
results_artifact_removal.group_mean             = round(mean(total_artifacts),2);
results_artifact_removal.group_SD               = round(std(total_artifacts),2);
results_artifact_removal.min                    = round(min(total_artifacts),2);
results_artifact_removal.max                    = round(max(total_artifacts),2);
results_artifact_removal.N                      = length(total_artifacts); % N = total # of subjects included in group-level results


% save group-level results 'results_artifact_removal'
full_file_name          = strcat('output','/','results_artifact_removal','.mat');
save(full_file_name, 'results_artifact_removal')

clearvars -except dir_parent



%% ------------------------------------------------------------------------
%  Part 4 - remove outlier trials
%  ------------------------------------------------------------------------



%  Part 4.1 subject-level -------------------------------------------------

% identify all data to analyze
% data to analyze is subject-level preprocessed data, constructed above in Part 2
input = dir(fullfile(strcat('derivatives','/deriv02_data_ICA'), '*.mat'));

for i = 1:length(input)

    % load input file i
    file_ICA = strcat(input(i).folder, '/',input(i).name);
    load(file_ICA)

    % extract subject id
    [~, subjectid, ~] = fileparts(input(i).name); % access subject id from current input file

    % run bad trial rejection on subject i
    [data_clean, bad_trials] = bad_trial_rej(data_ICA, subjectid);

    % save derivative data [1] 'data_clean'
    deriv                   = '/deriv04_data_clean/'; % update accordingly
    full_file_name          = strcat('derivatives',deriv,subjectid,'.mat');
    disp(strcat('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> saving_', deriv, subjectid,' <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<'));
    save(full_file_name, 'data_clean')

    % save derivative data [2] 'bad_trials'
    deriv                   = '/deriv05_bad_trials/'; % update accordingly
    full_file_name          = strcat('derivatives',deriv,subjectid,'.mat');
    disp(strcat('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> saving_', deriv, subjectid,' <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<'));
    save(full_file_name, 'bad_trials')


end



% Part 4.2 group-level ----------------------------------------------------

% identify all data to analyze
% data to analyze is subject-level rejected trials, constructed above in Part 4.1
input = dir(fullfile(strcat('derivatives','/deriv05_bad_trials'), '*.mat'));

% pre-allocate variables for speed
total_preserved_group           = zeros(1,length(input)); 
total_trials_group              = zeros(1,length(input)); 
perc_preserved_group            = zeros(1,length(input)); 

for i = 1:length(input)

    % load input file i 'bad_trials'
    file_bad_trials = strcat(input(i).folder, '/',input(i).name);
    load(file_bad_trials)

    % consolidate all subject data
    total_preserved_group(i)  = bad_trials.total_preserved;
    total_trials_group(i)     = bad_trials.total_trials;
    perc_preserved_group(i)   = bad_trials.perc_preserved;

end


% consolidate relevant results into one structure
results_bad_trials.total_preserved_trials   = total_preserved_group;
results_bad_trials.total_trials             = total_trials_group;
results_bad_trials.total_preserved_perc     = perc_preserved_group;
results_bad_trials.group_mean_trials        = round(mean(total_preserved_group),2);
results_bad_trials.group_mean_percent       = round(mean(perc_preserved_group),2);
results_bad_trials.group_SD_trials          = round(std(total_preserved_group),2);
results_bad_trials.group_SD_percent         = round(std(perc_preserved_group),2);
results_bad_trials.N                        = length(total_preserved_group); % N = total # of subjects included in group-level results


% save group-level results 'results_bad_trials'
full_file_name          = strcat('output','/','results_bad_trials','.mat');
save(full_file_name, 'results_bad_trials')

clearvars -except dir_parent



%% ------------------------------------------------------------------------
%  Part 5 - stratification
%  ------------------------------------------------------------------------

% segment clean data according to all 12 conditions (data_clean_12_cond)
% as well as by pooling across all normal and all scramble conditions (data_clean_pooled_NS)



% identify all data to analyze
% data to analyze is subject-level cleaned data, constructed above in Part 4
input = dir(fullfile(strcat('derivatives','/deriv04_data_clean'), '*.mat'));


% pre-allocate variables for speed
data_clean_12_cond              = cell(length(input),12); 
data_clean_pooled_NS            = cell(length(input),2); 

for i = 1:length(input)

    % load input file i
    file_clean = strcat(input(i).folder, '/',input(i).name);
    load(file_clean)

    % extract subject id
    [~, subjectid, ~] = fileparts(input(i).name); % access subject id from current input file



    %  Part 5.1 12-condition stratification -------------------------------------------------
    % segment data by all 12 conditions (bodies/faces/objects * human/monkey * normal/scramble)
    % re-structure data: store all clean data in a subject x condition cell array

    for condition = 1:12
        cfg = [];
        cfg.trials = data_clean.trialinfo == condition;
        data_clean_12_cond{i, condition} = ft_selectdata(cfg,data_clean);
    end



    %  Part 5.2 2-condition stratification -------------------------------------------------
    % segment data by 2 conditions (Normal v Scramble; pooled)
    % re-structure data: store all clean data in a subject x condition cell array

    % normal condition (column 1)
    cfg = [];
    cfg.trials = data_clean.trialinfo == 1 | data_clean.trialinfo == 2 | data_clean.trialinfo == 3 | data_clean.trialinfo == 4 | data_clean.trialinfo == 5 | data_clean.trialinfo == 6;
    data_clean_pooled_NS{i, 1} = ft_selectdata(cfg,data_clean);

    % scramble condition (column 2)
    cfg = [];
    cfg.trials = data_clean.trialinfo == 7 | data_clean.trialinfo == 8 | data_clean.trialinfo == 9 | data_clean.trialinfo == 10 | data_clean.trialinfo == 11 | data_clean.trialinfo == 12;
    data_clean_pooled_NS{i, 2} = ft_selectdata(cfg,data_clean);

end



% create separate variables for each condition (this reduces file size to optimize saving and loading)
data_clean_hum_body_norm    = data_clean_12_cond(:, 1);
data_clean_hum_face_norm    = data_clean_12_cond(:, 2); 
data_clean_hum_obj_norm     = data_clean_12_cond(:, 3); 

data_clean_monk_body_norm   = data_clean_12_cond(:, 4); 
data_clean_monk_face_norm   = data_clean_12_cond(:, 5); 
data_clean_monk_obj_norm    = data_clean_12_cond(:, 6); 

data_clean_hum_body_scr     = data_clean_12_cond(:, 7);  
data_clean_hum_face_scr     = data_clean_12_cond(:, 8); 
data_clean_hum_obj_scr      = data_clean_12_cond(:, 9); 

data_clean_monk_body_scr    = data_clean_12_cond(:, 10); 
data_clean_monk_face_scr    = data_clean_12_cond(:, 11); 
data_clean_monk_obj_scr     = data_clean_12_cond(:, 12); 

data_clean_pooled_normal    = data_clean_pooled_NS(:,1); 
data_clean_pooled_scramble  = data_clean_pooled_NS(:,2); 



% save all condition-level 'data_clean...' variables 
clear data_clean_12_cond
clear data_clean_pooled_NS
variables = who('*data_clean*');  % Get variables containing 'data_clean' in their name
for i = 1:length(variables)
    save(fullfile('output', [variables{i} '.mat']), variables{i});
end
















% %% Notes
%
%
% % % iterate through event numbers for each subject, save them as separate
% % % files, each names as the subject number e.g. 'subject01.mat':'subject31.mat'
% for i = 1:length(subjectinfo.eventNums)
%     event_codes = subjectinfo.eventNums{i};
%     subjectid = sprintf('subject%02d',i)
%     save(subjectid, 'event_codes')
% end
%


%       2. 'electrode_layout.mat', which describes the topographical parameters of the electrodes used in the present EEG recordings;
%           for more details, see: easycapM3 10%-based Electrode Layout
%           https://www.easycap.de/wp-content/uploads/2018/02/Easycap-10-based-electrode-layouts.pdf
%           https://www.fieldtriptoolbox.org/template/layout/#easycapm3---extended-1020-system-with-30-channels



