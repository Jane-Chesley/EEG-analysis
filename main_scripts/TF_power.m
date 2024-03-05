%% ------------------------------------------------------------------------
%  Time-frequency Power Analysis
%  ------------------------------------------------------------------------

% Script Description:
% This MATLAB script performs analysis of time-freuqency power of preprocessed EEG data.


% Usage:
%   - This script requires the following MATLAB toolboxes:
%       1. Image Processing Toolbox
%       2. Statistics and Machine Learning Toolbox
%       3. Bioinformatics Toolbox
%       4. FieldTrip Toolbox
%       5. fdr_bh
%   - Note: check which toolboxes are required by [1] running script [2] executing: [fList,pList] = matlab.codetools.requiredFilesAndProducts('run_analysis.m'); disp({pList.Name}');


% Inputs:
%   - Input data is stored in '/EEG-ERP/input' and includes:
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
%   2. ...

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

% issue a warning to the user if the path is not specified
if isempty(dir_FT)
    warning('dir_FT is empty. Please provide the path to the Fieldtrip Toolbox on your system.');
    % Prompt the user to assign a value to the variable
    dir_FT = input('Enter the path to Fieldtrip Toolbox on your system: ','s'); % 's' specifies input as string
end



% initialize the toolbox
addpath(dir_FT); ft_defaults;




%% ------------------------------------------------------------------------
%  Part 2 - load clean data
%  ------------------------------------------------------------------------

% identify all input data to analyze
% data to analyze are preprocessed EEG data, stratified by condition
input = dir(fullfile('output', '*data_clean*'));

% load all
disp('loading data...')
for i = 1: length(input)
    disp(strcat('loading data...', num2str(i),' of ', num2str(length(input))))
    load(strcat(input(i).folder, '/',input(i).name));
end






%% ------------------------------------------------------------------------
%  Part 3 - TF-specific preprocessing steps
%  ------------------------------------------------------------------------

% cfg: config TF-specific preprocessing steps
cfg                     = [];
cfg.lpfilter            = 'yes';
cfg.lpfreq              = 30;
cfg.hpfilter            = 'yes';
cfg.hpfreq              = 1;
cfg.detrend             = 'yes'; % removes drifts; use for TF analysis but not ERP

% execute preprocessing 
for s = 1:length(data_clean_hum_body_norm)

    % hum_body_norm
    data_TF_hum_body_norm{s}  = ft_preprocessing(cfg, data_clean_hum_body_norm{s});

    % hum_face_norm
    data_TF_hum_face_norm{s}  = ft_preprocessing(cfg, data_clean_hum_face_norm{s});

    % hum_obj_norm
    data_TF_hum_obj_norm{s}  = ft_preprocessing(cfg, data_clean_hum_obj_norm{s});

    % monk_body_norm
    data_TF_monk_body_norm{s}  = ft_preprocessing(cfg, data_clean_monk_body_norm{s});

    % monk_face_norm
    data_TF_monk_face_norm{s}  = ft_preprocessing(cfg, data_clean_monk_face_norm{s});

    % monk_obj_norm
    data_TF_monk_obj_norm{s}  = ft_preprocessing(cfg, data_clean_monk_obj_norm{s});

    % hum_body_scr
    data_TF_hum_body_scr{s}  = ft_preprocessing(cfg, data_clean_hum_body_scr{s});

    % hum_face_scr
    data_TF_hum_face_scr{s}  = ft_preprocessing(cfg, data_clean_hum_face_scr{s});

    % hum_obj_scr
    data_TF_hum_obj_scr{s}  = ft_preprocessing(cfg, data_clean_hum_obj_scr{s});

    % monk_body_scr
    data_TF_monk_body_scr{s}  = ft_preprocessing(cfg, data_clean_monk_body_scr{s});

    % monk_face_scr
    data_TF_monk_face_scr{s}  = ft_preprocessing(cfg, data_clean_monk_face_scr{s});

    % monk_obj_scr
    data_TF_monk_obj_scr{s}  = ft_preprocessing(cfg, data_clean_monk_obj_scr{s});

    % pooled_normal
    data_TF_pooled_normal{s}  = ft_preprocessing(cfg, data_clean_pooled_normal{s});

    % pooled_scramble
    data_TF_pooled_scramble{s}  = ft_preprocessing(cfg, data_clean_pooled_scramble{s});

end




%% ------------------------------------------------------------------------
%  Part 4 - TF transformation
%  ------------------------------------------------------------------------

% cfg: config TF analysis
cfg             = [];
cfg.output      = 'pow'; % power spectra
cfg.channel     = 'all'; % all channels (separately)
cfg.trials      = 'all';
cfg.foi         = 1:1:30; % all freqs (Hz) for now
cfg.keeptrials  = 'no'; % default, avg across condition
cfg.method      = 'wavelet'; % morlet wavelet transformation
cfg.width       = 3; % width of 7 yields time points only from -0.2, insufficient for baseline correction?
cfg.toi         = -0.6:0.05:1.7; % with this we go from 500 timepoints to 41


% execute TF transformation 
for s = 1:length(data_TF_hum_body_norm)

    % hum_body_norm
    data_TF_hum_body_norm{s}  = ft_freqanalysis(cfg, data_TF_hum_body_norm{s});

    % hum_face_norm
    data_TF_hum_face_norm{s}  = ft_freqanalysis(cfg, data_TF_hum_face_norm{s});

    % hum_obj_norm
    data_TF_hum_obj_norm{s}  = ft_freqanalysis(cfg, data_TF_hum_obj_norm{s});

    % monk_body_norm
    data_TF_monk_body_norm{s}  = ft_freqanalysis(cfg, data_TF_monk_body_norm{s});

    % monk_face_norm
    data_TF_monk_face_norm{s}  = ft_freqanalysis(cfg, data_TF_monk_face_norm{s});

    % monk_obj_norm
    data_TF_monk_obj_norm{s}  = ft_freqanalysis(cfg, data_TF_monk_obj_norm{s});

    % hum_body_scr
    data_TF_hum_body_scr{s}  = ft_freqanalysis(cfg, data_TF_hum_body_scr{s});

    % hum_face_scr
    data_TF_hum_face_scr{s}  = ft_freqanalysis(cfg, data_TF_hum_face_scr{s});

    % hum_obj_scr
    data_TF_hum_obj_scr{s}  = ft_freqanalysis(cfg, data_TF_hum_obj_scr{s});

    % monk_body_scr
    data_TF_monk_body_scr{s}  = ft_freqanalysis(cfg, data_TF_monk_body_scr{s});

    % monk_face_scr
    data_TF_monk_face_scr{s}  = ft_freqanalysis(cfg, data_TF_monk_face_scr{s});

    % monk_obj_scr
    data_TF_monk_obj_scr{s}  = ft_freqanalysis(cfg, data_TF_monk_obj_scr{s});

    % pooled_normal
    data_TF_pooled_normal{s}  = ft_freqanalysis(cfg, data_TF_pooled_normal{s});

    % pooled_scramble
    data_TF_pooled_scramble{s}  = ft_freqanalysis(cfg, data_TF_pooled_scramble{s});

end





%% ------------------------------------------------------------------------
%  Part 5 - baseline correction
%  ------------------------------------------------------------------------

% cfg: config baseline correction
cfg                 = [];
cfg.baseline        = [-0.6 -0.1];
cfg.baselinetype    = 'db';
cfg.parameter       = 'powspctrm';


% execute baseline correction 
for s = 1:length(data_TF_hum_body_norm)

    % hum_body_norm
    data_TF_hum_body_norm{s}  = ft_freqbaseline(cfg, data_TF_hum_body_norm{s});

    % hum_face_norm
    data_TF_hum_face_norm{s}  = ft_freqbaseline(cfg, data_TF_hum_face_norm{s});

    % hum_obj_norm
    data_TF_hum_obj_norm{s}  = ft_freqbaseline(cfg, data_TF_hum_obj_norm{s});

    % monk_body_norm
    data_TF_monk_body_norm{s}  = ft_freqbaseline(cfg, data_TF_monk_body_norm{s});

    % monk_face_norm
    data_TF_monk_face_norm{s}  = ft_freqbaseline(cfg, data_TF_monk_face_norm{s});

    % monk_obj_norm
    data_TF_monk_obj_norm{s}  = ft_freqbaseline(cfg, data_TF_monk_obj_norm{s});

    % hum_body_scr
    data_TF_hum_body_scr{s}  = ft_freqbaseline(cfg, data_TF_hum_body_scr{s});

    % hum_face_scr
    data_TF_hum_face_scr{s}  = ft_freqbaseline(cfg, data_TF_hum_face_scr{s});

    % hum_obj_scr
    data_TF_hum_obj_scr{s}  = ft_freqbaseline(cfg, data_TF_hum_obj_scr{s});

    % monk_body_scr
    data_TF_monk_body_scr{s}  = ft_freqbaseline(cfg, data_TF_monk_body_scr{s});

    % monk_face_scr
    data_TF_monk_face_scr{s}  = ft_freqbaseline(cfg, data_TF_monk_face_scr{s});

    % monk_obj_scr
    data_TF_monk_obj_scr{s}  = ft_freqbaseline(cfg, data_TF_monk_obj_scr{s});

    % pooled_normal
    data_TF_pooled_normal{s}  = ft_freqbaseline(cfg, data_TF_pooled_normal{s});

    % pooled_scramble
    data_TF_pooled_scramble{s}  = ft_freqbaseline(cfg, data_TF_pooled_scramble{s});

end



%% ------------------------------------------------------------------------
%  Part 6 - save subject-level TF data 
%  ------------------------------------------------------------------------

variables = who('*data_TF*');  % Get variables containing 'data_clean' in their name
for i = 1:length(variables)
    save(fullfile('derivatives', [variables{i} '.mat']), variables{i});
end



%% ------------------------------------------------------------------------
%  Part 7 - frequency grand average
%  ------------------------------------------------------------------------

% average power spectra data within subjects

% clear cfg
cfg = [];

% 12-condition FGA
FGA_hum_body_norm   = ft_freqgrandaverage(cfg,data_TF_hum_body_norm); 
FGA_hum_face_norm   = ft_freqgrandaverage(cfg,data_TF_hum_face_norm); 
FGA_hum_obj_norm   = ft_freqgrandaverage(cfg,data_TF_hum_obj_norm); 

FGA_monk_body_norm   = ft_freqgrandaverage(cfg,data_TF_monk_body_norm); 
FGA_monk_face_norm   = ft_freqgrandaverage(cfg,data_TF_monk_face_norm); 
FGA_monk_obj_norm   = ft_freqgrandaverage(cfg,data_TF_monk_obj_norm); 

FGA_hum_body_scr   = ft_freqgrandaverage(cfg,data_TF_hum_body_scr); 
FGA_hum_face_scr   = ft_freqgrandaverage(cfg,data_TF_hum_face_scr); 
FGA_hum_obj_scr   = ft_freqgrandaverage(cfg,data_TF_hum_obj_scr); 

FGA_monk_body_scr   = ft_freqgrandaverage(cfg,data_TF_monk_body_scr);
FGA_monk_face_scr   = ft_freqgrandaverage(cfg,data_TF_monk_face_scr); 
FGA_monk_obj_scr   = ft_freqgrandaverage(cfg,data_TF_monk_obj_scr); 

% 2-condition (Normal v. Scramble; pooled) FGA
FGA_pooled_normal   = ft_freqgrandaverage(cfg,data_TF_pooled_normal); 
FGA_pooled_scramble   = ft_freqgrandaverage(cfg,data_TF_pooled_scramble); 



%% ------------------------------------------------------------------------
%  Part 8 - save group-level TF data 
%  ------------------------------------------------------------------------

variables = who('*FGA*');  % Get variables containing 'data_clean' in their name
for i = 1:length(variables)
    save(fullfile('output', [variables{i} '.mat']), variables{i});
end












%% Notes
% % % % % REMEMBER

% 2. detrend for TF analyses
% cfg.detrend = 'yes'; % removes drifts, use for TFR but not ERPs

% 3. downsample AFTER lowpass filtering (or possibly not at all)
% % downsample
% cfg                     = [];
% cfg.resamplefs          = 250; % must be a multiple of sampling freq
% cfg.method              = 'downsample';
% data                    = ft_resampledata(cfg, data);






