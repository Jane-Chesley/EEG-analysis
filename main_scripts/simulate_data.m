%% Connectivity analysis on simulated data
%  ------------------------------------------------------------------------
%
% Script Description:
% This MATLAB script constructs simulated EEG data and performs connectivity analysis on it.
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
%   - Real, pre-processed EEG data
%
% Author:
%   Jane Chesley
%
% Affiliation:
%   Faculty of Psychology & Neuroscience, Maastricht University
%   Maastricht, the Netherlands
%
% Date:
%   15 April, 2024
%
% Outline of script content:
%   1. Script setup
%   2. Construct simulated data 
%   3. PLI computation on simulated data 
% ------------------------------------------------------------------------



% tic;

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
%  Part 2 - Construct simulated data
%  ------------------------------------------------------------------------

% identify input data: preprocessed EEG data
% data_clean_hum_body_normal contains all subject data (29) for one example condition (human body normal)
load('output/data_clean_hum_body_norm.mat');


% construct simulated data for all subjects, such that:
% channel 1 = random noise A (values between 0 and 1)
% channel 2 = random noise B (values between 0 and 1)
% channel 3 = sinusoid signal A (theta)
% channel 4 = sinusoid signal A (theta)
% channel 5 = sinusoid signal B (alpha)
% channel 6 = sinusoid signal B (alpha)
% all remaining channels = real data

% pre-allocate variable
% sim_data.trial{trial}(channel,timepts)
sim_data = data_clean_hum_body_norm;

for s = 1:length(sim_data) % loop for each subject

    for tr = 1:length(sim_data{s}.trial) % loop for each trial


        dim = size(sim_data{s}.trial{tr}(1,:)); % get dimensions of all time points for a given trial (tr) and channel 1

        % simulate random noise A (random values between 0 and 1)
        rand_A = randn(dim); % create a new vector of random noise, with the same length as the original time points
        sim_data{s}.trial{tr}(1,:) = rand_A; % assign it to channel 1

        % simulate random noise B
        rand_B = rand(dim); % create a new vector of random noise; this random noise is different from channel 1
        sim_data{s}.trial{tr}(2,:) = rand_B; % assign it to channel 2

        % simulate sinusoid signal A
        amplitude = 1; % we are not manipulating amplitude, but it is included to understand the equation
        frequency = 5; % theta
        timepts = sim_data{s}.time{tr}; % get time points (in seconds) for each trial (575 time points per trial)
        sin_A = amplitude*sin(2*pi*frequency*timepts); % construct sinusoid signal
        sim_data{s}.trial{tr}(3,:) = sin_A; % assign it to channel 3

        % simulate sinusoid signal B; similar phase to sin A (but not identical because PLI attenuates this)
        phase_difference = pi/128;  % Phase difference between the two signals
        sin_B = amplitude*sin(2*pi*frequency*timepts + phase_difference); % construct sinusoid signal with a phase difference
        sim_data{s}.trial{tr}(4,:) = sin_B; % assign it to channel 4

        % simulate sinusoid signal A
        amplitude = 1; % we are not manipulating amplitude, but it is included to understand the equation
        frequency = 10; % alpha
        timepts = sim_data{s}.time{tr}; % get time points (in seconds) for each trial (575 time points per trial)
        sin_C = amplitude*sin(2*pi*frequency*timepts); % construct sinusoid signal
        sim_data{s}.trial{tr}(5,:) = sin_C; % assign it to channel 5

        % simulate sinusoid signal B; similar phase to sin A (but not identical because PLI attenuates this)
        phase_difference = pi/128;  % Phase difference between the two signals
        sin_D = amplitude*sin(2*pi*frequency*timepts + phase_difference); % construct sinusoid signal with a phase difference
        sim_data{s}.trial{tr}(6,:) = sin_D; % assign it to channel 6

    end

end

%% ------------------------------------------------------------------------
%  Part 3 - PLI computation on simulated data
%  ------------------------------------------------------------------------

% expectations of simulated results:
% channel pairs 1-2 (noise A - noise B) not coherent (PLI = 0)
% channel pairs 1-3 (noise A - sin A) not coherent (PLI = 0)
% channel pairs 3-4 (sin A - sin B) close to perfectly coherent in theta band (PLI = 1)
% channel pairs 5-6 (sin C - sin D) close to perfectly coherent in alpha band (PLI = 1)

% define toi range in SECONDS
min_T1 = 0;
max_T1 = 1;

% define all frequency bands
all_freq = {'delta';'theta';'alpha';'beta';'gamma_A';'gamma_B'};

% get all channel labels for roi extraction
% channel labels are constant across all measurements
all_channels = data_clean_hum_body_norm{1}.label;

% run PLI computation
sim_PLI = PLI_all_subjects(sim_data, min_T1, max_T1);

%% expectations of simulated results:
% channel pairs 1-2 (noise A - noise B) not coherent (PLI = 0)
roi_idx = [1;2];

% compute within-region PLI separately for each frequency band
for freq = 1:length(all_freq)
    freq_string = all_freq{freq};
    sim_PLI_allfreq(freq)= mean(PLI_within(sim_PLI, roi_idx, freq));
end

disp(' ')
disp(' ')
disp('PLI of channel pairs 1-2 (noise A - noise B)')
disp('Expectation = not coherent (PLI = 0).')
disp(' ')
disp(strcat('simulated PLI delta:', num2str(sim_PLI_allfreq(1))))
disp(strcat('simulated PLI theta:', num2str(sim_PLI_allfreq(2))))
disp(strcat('simulated PLI alpha:', num2str(sim_PLI_allfreq(3))))
disp(strcat('simulated PLI beta:', num2str(sim_PLI_allfreq(4))))
disp(strcat('simulated PLI gamma A:', num2str(sim_PLI_allfreq(5))))
disp(strcat('simulated PLI gamma B:', num2str(sim_PLI_allfreq(6))))



%% expectations of simulated results:
% channel pairs 1-3 (noise A - sin A) not coherent (PLI = 0)
roi_idx = [1;3];
clear sim_PLI_allfreq 

% compute within-region PLI separately for each frequency band
for freq = 1:length(all_freq)
    freq_string = all_freq{freq};
    sim_PLI_allfreq(freq)= mean(PLI_within(sim_PLI, roi_idx, freq));
end

disp(' ')
disp(' ')
disp('PLI of channel pairs 1-3 (noise A - sin A)')
disp('Expectation = not coherent (PLI = 0).')
disp(' ')
disp(strcat('simulated PLI delta:', num2str(sim_PLI_allfreq(1))))
disp(strcat('simulated PLI theta:', num2str(sim_PLI_allfreq(2))))
disp(strcat('simulated PLI alpha:', num2str(sim_PLI_allfreq(3))))
disp(strcat('simulated PLI beta:', num2str(sim_PLI_allfreq(4))))
disp(strcat('simulated PLI gamma A:', num2str(sim_PLI_allfreq(5))))
disp(strcat('simulated PLI gamma B:', num2str(sim_PLI_allfreq(6))))



%% expectations of simulated results:
% channel pairs 3-4 (sin A - sin B) close to perfectly coherent in theta band (PLI = 1)
roi_idx = [3;4];
clear sim_PLI_allfreq 

% compute within-region PLI separately for each frequency band
for freq = 1:length(all_freq)
    freq_string = all_freq{freq};
    sim_PLI_allfreq(freq)= mean(PLI_within(sim_PLI, roi_idx, freq));
end

disp(' ')
disp(' ')
disp('PLI of channel pairs 3-4 (sin A - sin B)')
disp('Expectation = close to perfectly coherent in theta band (PLI = 1).')
disp(' ')
disp(strcat('simulated PLI delta:', num2str(sim_PLI_allfreq(1))))
disp(strcat('simulated PLI theta:', num2str(sim_PLI_allfreq(2))))
disp(strcat('simulated PLI alpha:', num2str(sim_PLI_allfreq(3))))
disp(strcat('simulated PLI beta:', num2str(sim_PLI_allfreq(4))))
disp(strcat('simulated PLI gamma A:', num2str(sim_PLI_allfreq(5))))
disp(strcat('simulated PLI gamma B:', num2str(sim_PLI_allfreq(6))))



%% expectations of simulated results:
% channel pairs 5-6 (sin C - sin D) close to perfectly coherent in theta band (PLI = 1)
roi_idx = [5;6];
clear sim_PLI_allfreq 

% compute within-region PLI separately for each frequency band
for freq = 1:length(all_freq)
    freq_string = all_freq{freq};
    sim_PLI_allfreq(freq)= mean(PLI_within(sim_PLI, roi_idx, freq));
end

disp(' ')
disp(' ')
disp('PLI of channel pairs 5-6 (sin C - sin D)')
disp('Expectation = close to perfectly coherent in alpha band (PLI = 1).')
disp(' ')
disp(strcat('simulated PLI delta:', num2str(sim_PLI_allfreq(1))))
disp(strcat('simulated PLI theta:', num2str(sim_PLI_allfreq(2))))
disp(strcat('simulated PLI alpha:', num2str(sim_PLI_allfreq(3))))
disp(strcat('simulated PLI beta:', num2str(sim_PLI_allfreq(4))))
disp(strcat('simulated PLI gamma A:', num2str(sim_PLI_allfreq(5))))
disp(strcat('simulated PLI gamma B:', num2str(sim_PLI_allfreq(6))))





