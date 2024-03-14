%% ------------------------------------------------------------------------
%  Time-frequency Power Analysis
%  ------------------------------------------------------------------------

% Script Description:
% This MATLAB script performs analysis of time-frequency power of preprocessed EEG data.


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
    disp(strcat('loading data...', num2str(i),' of ', num2str(length(input))));
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

% pre-allocate variables for speed
N = length(data_clean_hum_body_norm);
data_TF_hum_body_norm       = cell(1,N); 
data_TF_hum_face_norm       = cell(1,N); 
data_TF_hum_obj_norm        = cell(1,N); 
data_TF_monk_body_norm      = cell(1,N); 
data_TF_monk_face_norm      = cell(1,N); 
data_TF_monk_obj_norm       = cell(1,N); 
data_TF_hum_body_scr        = cell(1,N); 
data_TF_hum_face_scr        = cell(1,N); 
data_TF_hum_obj_scr         = cell(1,N); 
data_TF_monk_body_scr       = cell(1,N); 
data_TF_monk_face_scr       = cell(1,N); 
data_TF_monk_obj_scr        = cell(1,N); 
data_TF_pooled_normal       = cell(1,N); 
data_TF_pooled_scramble     = cell(1,N); 

% execute preprocessing 
for s = 1:N

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

    msgbox(strcat('analyzing data...', num2str(s)," of ", num2str(length(data_TF_hum_body_norm))));

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
    
    msgbox(strcat('analyzing data...', num2str(s)," of ", num2str(length(data_TF_hum_body_norm))));

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
    save(fullfile('output', [variables{i} '.mat']), variables{i});
end




%% ------------------------------------------------------------------------
%  Part 7 - cluster-based permutation analysis  
%  ------------------------------------------------------------------------

% 7.1
% input = TF data pooled across normal and scramble 
load(fullfile('output','data_TF_pooled_normal.mat'))
load(fullfile('output','data_TF_pooled_scramble.mat'))
load(fullfile('input','easycapM3.mat'))

% 7.2
% FGA 
cfg = [];
cfg.keepindividual = 'yes';

GA_TFR_NOR          = ft_freqgrandaverage(cfg,data_TF_pooled_normal{:});
GA_TFR_SCR          = ft_freqgrandaverage(cfg,data_TF_pooled_scramble{:});



% 7.3 
% define neighbours
cfg = [];
cfg.channel             = 'all';
cfg.method              = 'distance';
cfg.layout              = lay;
neighbours              = ft_prepare_neighbours(cfg, GA_TFR_NOR);

% cfg.method specifies method for constructing neighboring channels;
% 'triangulation' = nearest direct neighbors (more conservative); algorithm builds triangles between nearby nodes
% 'distance' = electrodes within 3-D Euclidean distance; draws a circle around each electrode
% 'template' = defined by ft_neighbourplot --> https://www.fieldtriptoolbox.org/template/neighbours/


% 7.4
% construct design matrix %NxConditions
subj = 29;
design = zeros(2,2*subj);
for i = 1:subj
    design(1,i) = i;
end
for i = 1:subj
    design(1,subj+i) = i;
end
design(2,1:subj)        = 1;
design(2,subj+1:2*subj) = 2;



% 7.5 
% Cluster analysis 
% input FGA
% output stat --> posclusters 

cfg = [];
cfg.channel          = 'all';

% ATTN 
cfg.frequency        = [4 7]; %[4 7]; % theta
cfg.avgoverfreq      = 'yes'; % 'yes';
cfg.latency          = [0 1]; % [0.2 0.4]; % in seconds %%%%%%%%% NOT MS 
cfg.avgovertime      = 'yes';  %'yes';


cfg.method           = 'montecarlo'; % Monte-Carlo estimates of critical vals from permutation distribution
cfg.statistic        = 'ft_statfun_depsamplesT'; % dependent samples T-statistic for within-subjects design
cfg.correctm         = 'cluster'; % correction for multiple comparisons ('cluster', 'bonferroni', 'fdr', 'no')
% 'cluster' calculates a  cluster-based test statistic and its significance probability --> maximum of the sum of t-values within every cluster
cfg.clusteralpha     = 0.05; % critical value
cfg.clusterstatistic = 'maxsum'; % how to combine samples that belong to a cluster ('maxsum','maxsize','wcm' /weighted cluster mass)
cfg.clusterthreshold = 'parametric'; % method for single-sample threshold ('parametric', 'nonparametric_individual', 'nonparametric_common')
cfg.minnbchan        = 2; % min # of channels in a cluster
cfg.tail             = 1; % -1 (one-tailed neg) / 1 (one-tailed pos) / 0 (two-tailed)
cfg.clustertail      = 1;
cfg.numrandomization = 10000; % above 1000 is good 
cfg.neighbours       = neighbours; % constructed above  from FT_PREPARE_NEIGHBOURS % data input can be any structure with channel info
cfg.design           = design; % NxConditions; constructed above
cfg.uvar             = 1;
cfg.ivar             = 2;

[stat] = ft_freqstatistics(cfg, GA_TFR_NOR, GA_TFR_SCR);

results_clusterperm.clusteralpha = cfg.clusteralpha;
results_clusterperm.minclustersize = cfg.minnbchan;
results_clusterperm.numrandomization = cfg.numrandomization; 
results_clusterperm.tail = cfg.tail; 
results_clusterperm.timewindow_seconds = cfg.latency;
results_clusterperm.freq = cfg.frequency;


% 7.6 
% Plot the results

cfg = [];
% cfg.zlim   = [-4 4];
cfg.highlightseries         = {'on', 'on'} ; % on for p < [0.01 0.05]
cfg.highlightsymbolseries   = ['*','*']; % asterisks for p < [0.01 0.05]
cfg.highlightsizeseries     = [10,10];
cfg.subplotsize             = [1 1];
cfg.layout                  = 'easycapM3.mat';
cfg.alpha                   = 0.05;
cfg.saveaspng               = 'Clusterplot_theta_0_1000ms';
cfg.colorbar                = 'yes';
cfg.colormap                = '-RdGy';
fig = ft_clusterplot(cfg, stat); 
% title('clusterplot (theta/200-400ms)');



% 7.7 
% Output

% channels > 0.05; 0-1000ms; normal vs scramble
channel_idx = find(stat.posclusterslabelmat == 1);

% all channels 
results_clusterperm.sig_channels = stat.cfg.channel(channel_idx); 

% save results 
save(fullfile('output', 'results_clusterperm.mat'));




%% ------------------------------------------------------------------------
%  Part 9 - Extract TF power from ROI (channels), FOI (frequency band), TOI
%  (time-window) for statistical analysis 
%  ------------------------------------------------------------------------

% identify all input data to analyze
% data to analyze are preprocessed EEG data, stratified by condition
input = dir(fullfile('output','*data_TF*'));

% load all 
for i = 1: length(input)
    load(strcat(input(i).folder, '/',input(i).name));
end

load(fullfile('output','results_clusterperm.mat'));


%% define parameters of interest 

% ROI
% roi = {'C3', 'CP3', 'P3'}; %% 'human body selective cluster'
% roi = results_clusterperm.sig_channels; % '*dynamic* object-level processing cluster'
all_channels = data_TF_hum_body_norm{1}.label; % channel labels are constant across all measurements 
roi_idx = find(ismember(all_channels, roi)); 
roi_string = 'roiHumanBody';

% FOI
foi_idx = 4:7; % theta 
foi_string = 'theta'; 

% TOI (in SECONDS)
min_T = 0.2;
max_T = 0.55;
toi_string = '200-550ms'; 

all_timepts = data_TF_hum_body_norm{1}.time; % all time points, which are constant across subjects 
all_timepts = round(all_timepts,2); % round time points 
toi_idx = find(all_timepts == min_T, 1, 'first'):find(all_timepts == max_T, 1, 'last'); % find indices of toi range 

% Extract parameters of interest 
for s = 1:length(data_TF_hum_body_norm) % length (subject numbers) is constant across datasets 

    data_allsubjects(s,1) = mean(data_TF_hum_body_norm{s}.powspctrm(roi_idx, foi_idx, toi_idx), 'all'); % average across ROI/FOI/TOI for each subject 
    data_allsubjects(s,2) = mean(data_TF_hum_face_norm{s}.powspctrm(roi_idx, foi_idx, toi_idx), 'all'); % average across ROI/FOI/TOI for each subject 
    data_allsubjects(s,3) = mean(data_TF_hum_obj_norm{s}.powspctrm(roi_idx, foi_idx, toi_idx), 'all'); % average across ROI/FOI/TOI for each subject 
    
    data_allsubjects(s,4) = mean(data_TF_monk_body_norm{s}.powspctrm(roi_idx, foi_idx, toi_idx), 'all'); % average across ROI/FOI/TOI for each subject 
    data_allsubjects(s,5) = mean(data_TF_monk_face_norm{s}.powspctrm(roi_idx, foi_idx, toi_idx), 'all'); % average across ROI/FOI/TOI for each subject 
    data_allsubjects(s,6) = mean(data_TF_monk_obj_norm{s}.powspctrm(roi_idx, foi_idx, toi_idx), 'all'); % average across ROI/FOI/TOI for each subject 
    
    data_allsubjects(s,7) = mean(data_TF_hum_body_scr{s}.powspctrm(roi_idx, foi_idx, toi_idx), 'all'); % average across ROI/FOI/TOI for each subject 
    data_allsubjects(s,8) = mean(data_TF_hum_face_scr{s}.powspctrm(roi_idx, foi_idx, toi_idx), 'all'); % average across ROI/FOI/TOI for each subject 
    data_allsubjects(s,9) = mean(data_TF_hum_obj_scr{s}.powspctrm(roi_idx, foi_idx, toi_idx), 'all'); % average across ROI/FOI/TOI for each subject 
    
    data_allsubjects(s,10) = mean(data_TF_monk_body_scr{s}.powspctrm(roi_idx, foi_idx, toi_idx), 'all'); % average across ROI/FOI/TOI for each subject 
    data_allsubjects(s,11) = mean(data_TF_monk_face_scr{s}.powspctrm(roi_idx, foi_idx, toi_idx), 'all'); % average across ROI/FOI/TOI for each subject 
    data_allsubjects(s,12) = mean(data_TF_monk_obj_scr{s}.powspctrm(roi_idx, foi_idx, toi_idx), 'all'); % average across ROI/FOI/TOI for each subject 

end 

% Convert the matrix into a table with variable names
variable_names = {'hum_body_norm','hum_face_norm','hum_obj_norm','monk_body_norm','monk_face_norm','monk_obj_norm','hum_body_scr','hum_face_scr','hum_obj_scr','monk_body_scr','monk_face_scr','monk_obj_scr'};
T = array2table(data_allsubjects, 'VariableNames', variable_names);

% Write the table to an Excel file
writetable(T, strcat('Stats_TFpower_allconditions_',roi_string,'_',foi_string,'_',toi_string,'.xlsx'));



% compute normalization (normal - scramble)
TFpower_normalized = zeros(29,6);
for c = 1:6
    TFpower_normalized(:,c) = data_allsubjects(:,c)-data_allsubjects(:,c+6);
end 
variable_names2 = {'hum_body','hum_face','hum_obj','monk_body','monk_face','monk_obj'};

T2 = array2table(TFpower_normalized, 'VariableNames', variable_names2);

% Write the table to an Excel file
writetable(T2, strcat('Stats_TFpower_normalized_',roi_string,'_',foi_string,'_',toi_string,'.xlsx'));












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






