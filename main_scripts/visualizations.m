%% ------------------------------------------------------------------------
%  Figures
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



%% identify all input data to analyze
% data to analyze are preprocessed EEG data, stratified by condition
input = dir(fullfile('output', '*data_TF*'));

% load all
disp('loading data...')
for i = 1: length(input)
    disp(strcat('loading data...', num2str(i),' of ', num2str(length(input))));
    load(strcat(input(i).folder, '/',input(i).name));
end

% channel layout for topoplots 
load(fullfile('input','easycapM3.mat'));
load(fullfile('input','results_clusterperm.mat'));





%% ------------------------------------------------------------------------
%  Part 2 - Time-series plots 
%  ------------------------------------------------------------------------

%% Define parameters of interest 

% ROI 
roi = {'C3', 'CP3', 'P3'}; %% 'human body selective cluster'
% roi = results_clusterperm.sig_channels; % '*dynamic* object-level processing cluster'

% FOI 
foi = [4 7]; % theta 

% TOI 
% all time points !! 



%% Extract parameters of interest 
% average across subjects; ROI; FOI 

cfg = [];
cfg.keepindividual = 'no';
cfg.foilim         = foi; 
cfg.toilim         = [-0.1 1]; 
cfg.channel        = roi;

% 12-condition FGA
FGA_hum_body_norm   = ft_freqgrandaverage(cfg,data_TF_hum_body_norm{1,:}); 
FGA_hum_face_norm   = ft_freqgrandaverage(cfg,data_TF_hum_face_norm{1,:}); 
FGA_hum_obj_norm   = ft_freqgrandaverage(cfg,data_TF_hum_obj_norm{1,:}); 

FGA_monk_body_norm   = ft_freqgrandaverage(cfg,data_TF_monk_body_norm{1,:}); 
FGA_monk_face_norm   = ft_freqgrandaverage(cfg,data_TF_monk_face_norm{1,:}); 
FGA_monk_obj_norm   = ft_freqgrandaverage(cfg,data_TF_monk_obj_norm{1,:}); 

FGA_hum_body_scr   = ft_freqgrandaverage(cfg,data_TF_hum_body_scr{1,:}); 
FGA_hum_face_scr   = ft_freqgrandaverage(cfg,data_TF_hum_face_scr{1,:}); 
FGA_hum_obj_scr   = ft_freqgrandaverage(cfg,data_TF_hum_obj_scr{1,:}); 

FGA_monk_body_scr   = ft_freqgrandaverage(cfg,data_TF_monk_body_scr{1,:}); 
FGA_monk_face_scr   = ft_freqgrandaverage(cfg,data_TF_monk_face_scr{1,:}); 
FGA_monk_obj_scr   = ft_freqgrandaverage(cfg,data_TF_monk_obj_scr{1,:}); 

% 2-condition (Normal v. Scramble; pooled) FGA
FGA_pooled_normal   = ft_freqgrandaverage(cfg,data_TF_pooled_normal{1,:}); 
FGA_pooled_scramble   = ft_freqgrandaverage(cfg,data_TF_pooled_scramble{1,:}); 


% relevant output of ft_freqgrandaverage is powspctrm (chan x freq x time)
% average across chan and freq (dim 1:2) for each condition 
time_series_hum_body_norm = squeeze(mean(FGA_hum_body_norm.powspctrm,1:2));
time_series_hum_face_norm = squeeze(mean(FGA_hum_face_norm.powspctrm,1:2));
time_series_hum_obj_norm = squeeze(mean(FGA_hum_obj_norm.powspctrm,1:2));

time_series_monk_body_norm = squeeze(mean(FGA_monk_body_norm.powspctrm,1:2));
time_series_monk_face_norm = squeeze(mean(FGA_monk_face_norm.powspctrm,1:2));
time_series_monk_obj_norm = squeeze(mean(FGA_monk_obj_norm.powspctrm,1:2));

time_series_hum_body_scr = squeeze(mean(FGA_hum_body_scr.powspctrm,1:2));
time_series_hum_face_scr = squeeze(mean(FGA_hum_face_scr.powspctrm,1:2));
time_series_hum_obj_scr = squeeze(mean(FGA_hum_obj_scr.powspctrm,1:2));

time_series_monk_body_scr = squeeze(mean(FGA_monk_body_scr.powspctrm,1:2));
time_series_monk_face_scr = squeeze(mean(FGA_monk_face_scr.powspctrm,1:2));
time_series_monk_obj_scr = squeeze(mean(FGA_monk_obj_scr.powspctrm,1:2));

time_series_pooled_normal = squeeze(mean(FGA_pooled_normal.powspctrm,1:2));
time_series_pooled_scramble = squeeze(mean(FGA_pooled_scramble.powspctrm,1:2));


%% normalization 
% normalization (subtraction Normal - Scramble)
time_series_diff_hum_body = time_series_hum_body_norm-time_series_hum_body_scr;
time_series_diff_hum_face = time_series_hum_face_norm-time_series_hum_face_scr;
time_series_diff_hum_obj = time_series_hum_obj_norm-time_series_hum_obj_scr;

time_series_diff_monk_body = time_series_monk_body_norm-time_series_monk_body_scr;
time_series_diff_monk_face = time_series_monk_face_norm-time_series_monk_face_scr;
time_series_diff_monk_obj = time_series_monk_obj_norm-time_series_monk_obj_scr;


%% time points 
% define vector of time points in ms to use for plots
timepts = FGA_hum_body_norm.time; % get timepts, which are constant across all datasets
timepts = round(timepts,2); % round to 2 decimals
timepts = timepts*1000; % ms to seconds


%% figure 1 

% plot difference bodies 
close all 
fig = figure; 
plot(timepts,time_series_diff_hum_body,"Color",'k','LineWidth',5) %%%%%%%%%%%%%%%%%%% plot line 1 HUMANS 

hold on
set(gca,'xtick',-200:200:1000,'FontSize',18,'FontName','Calibri','FontWeight','bold');
ylim([-0.7 1.3]);
xlim([-150 1100]); 
xticklabels({-200:200:1000}); 
set(gca,'XAxisLocation', 'origin');
set(gca,'YAxisLocation', 'origin');
% xlabel('Time (s)','FontSize',10)
% ylabel('Theta power (dB)','FontSize',10)
set(gca,'linewidth',1.5)
set(gca,'box','off');
title('Differential bodies'); 

hold on 
plot(timepts,time_series_diff_monk_body,"Color",'k','LineWidth',4,'LineStyle',':') %%%%%%%%%%%%%%%%%%% plot line 2 MONKEYS 
leg = legend({'Human', 'Monkey'});

% % save 
% saveas(fig,'TimeSeries_difference_Bodies.png')



%% figure 2 

% plot difference faces 
difference = 3;
close all 
fig = figure; 
plot(timepts,time_series_diff_hum_face,"Color",'k','LineWidth',5) %%%%%%%%%%%%%%%%%%% plot line 1 HUMANS 

hold on
set(gca,'xtick',-200:200:1000,'FontSize',18,'FontName','Calibri','FontWeight','bold');
set(gca,'ytick',-1:0.5:1.5)
ylim([-0.7 1.3]);
xlim([-150 1100]); 
xticklabels({-200:200:1000}); 
set(gca,'XAxisLocation', 'origin');
set(gca,'YAxisLocation', 'origin');
% xlabel('Time (s)','FontSize',10)
% ylabel('Theta power (dB)','FontSize',10)
set(gca,'linewidth',1.5)
set(gca,'box','off');
title('Differential faces'); 

hold on 
plot(timepts,time_series_diff_monk_face,"Color",'k','LineWidth',4,'LineStyle',':') %%%%%%%%%%%%%%%%%%% plot line 2 MONKEYS 
leg = legend({'Human', 'Monkey'});

% % save 
% saveas(fig,'TimeSeries_difference_Faces.png')



%% figure 3

% plot difference obj 
difference = 3;
close all 
fig = figure; 
plot(timepts,time_series_diff_hum_obj,"Color",'k','LineWidth',5) %%%%%%%%%%%%%%%%%%% plot line 1 HUMANS 

hold on
set(gca,'xtick',-200:200:1000,'FontSize',18,'FontName','Calibri','FontWeight','bold');
set(gca,'ytick',-1:0.5:1.5)
ylim([-0.7 1.3]);
xlim([-150 1100]); 
xticklabels({-200:200:1000}); 
set(gca,'XAxisLocation', 'origin');
set(gca,'YAxisLocation', 'origin');
% xlabel('Time (s)','FontSize',10)
% ylabel('Theta power (dB)','FontSize',10)
set(gca,'linewidth',1.5)
set(gca,'box','off');
title('Differential obj'); 

hold on 
plot(timepts,time_series_diff_monk_obj,"Color",'k','LineWidth',4,'LineStyle',':') %%%%%%%%%%%%%%%%%%% plot line 2 MONKEYS 
leg = legend({'Human', 'Monkey'});

% % save 
% saveas(fig,'TimeSeries_difference_Objects.png')



%% figure 4 

GA_diff_species_bodies = time_series_diff_hum_body-time_series_diff_monk_body; 

% plot 
close all 
fig = figure; 
plot(timepts,time_series_diff_hum_body,"Color",'#7D7D7D','LineWidth',5) %%%%%%%%%%%%%%%%%%% plot line 1 HUMAN 

hold on
set(gca,'xtick',-200:200:1000,'FontSize',18,'FontName','Calibri','FontWeight','bold');
ylim([-0.7 1.3]);
xlim([-150 1100]); 
xticklabels({-200:200:1000}); 
set(gca,'XAxisLocation', 'origin');
set(gca,'YAxisLocation', 'origin');
% xlabel('Time (s)','FontSize',10)
% ylabel('Theta power (dB)','FontSize',10)
set(gca,'linewidth',1.5)
set(gca,'box','off');
title('Differential bodies'); 

hold on 
plot(timepts,time_series_diff_monk_body,"Color",'#7D7D7D','LineWidth',4,'LineStyle',':') %%%%%%%%%%%%%%%%%%% plot line 2 MONKEY 


hold on 
plot(timepts,GA_diff_species_bodies,"Color",'#d72013','LineWidth',5) %%%%%%%%%%%%%%%%%%% plot line 3 HUMAN-MONKEY 
leg = legend({'Human', 'Monkey','Human-Monkey'}, 'FontSize',18);

% % save 
saveas(fig,'TimeSeries_difference_species_bodies.png')
% saveas(fig,'TimeSeries_difference_species_bodies.fig')




%% ------------------------------------------------------------------------
%  Part 3 - Time-frequency plots
%  ------------------------------------------------------------------------



cfg = [];

% 12-condition FGA
FGA_hum_body_norm   = ft_freqgrandaverage(cfg,data_TF_hum_body_norm{1,:}); 
FGA_hum_face_norm   = ft_freqgrandaverage(cfg,data_TF_hum_face_norm{1,:}); 
FGA_hum_obj_norm   = ft_freqgrandaverage(cfg,data_TF_hum_obj_norm{1,:}); 

FGA_monk_body_norm   = ft_freqgrandaverage(cfg,data_TF_monk_body_norm{1,:}); 
FGA_monk_face_norm   = ft_freqgrandaverage(cfg,data_TF_monk_face_norm{1,:}); 
FGA_monk_obj_norm   = ft_freqgrandaverage(cfg,data_TF_monk_obj_norm{1,:}); 

FGA_hum_body_scr   = ft_freqgrandaverage(cfg,data_TF_hum_body_scr{1,:}); 
FGA_hum_face_scr   = ft_freqgrandaverage(cfg,data_TF_hum_face_scr{1,:}); 
FGA_hum_obj_scr   = ft_freqgrandaverage(cfg,data_TF_hum_obj_scr{1,:}); 

FGA_monk_body_scr   = ft_freqgrandaverage(cfg,data_TF_monk_body_scr{1,:}); 
FGA_monk_face_scr   = ft_freqgrandaverage(cfg,data_TF_monk_face_scr{1,:}); 
FGA_monk_obj_scr   = ft_freqgrandaverage(cfg,data_TF_monk_obj_scr{1,:}); 

% 2-condition (Normal v. Scramble; pooled) FGA
FGA_pooled_normal   = ft_freqgrandaverage(cfg,data_TF_pooled_normal{1,:}); 
FGA_pooled_scramble   = ft_freqgrandaverage(cfg,data_TF_pooled_scramble{1,:}); 



%% TF plot; 0-1000 ms; 1-30 Hz; ROI=DynObjLvl


% zlim = [-2 2];
% 
% % human body normal 
% cfg = [];
% cfg.xlim           = [-0.1 1] % before -0.1 s is not visible (Morlet wavelets)
% cfg.ylim           = [4 30] % below 4 Hz is not visible (Morlet wavelets)
% cfg.zlim           = zlim; 
% cfg.fontsize       = 14
% cfg.title          = 'Human Body Normal'
% cfg.colormap       = '-RdGy';
% cfg.channel        = roi;
% 
% fig1 = figure;
% ft_singleplotTFR(cfg, FGA_hum_body_norm);
% 
% % % save
% % saveas(fig1,'Fig2_allchan_normal.png');
% 
% 
% % human body scramble
% cfg = [];
% cfg.xlim           = [-0.1 1] % before -0.1 s is not visible (Morlet wavelets)
% cfg.ylim           = [4 30] % below 4 Hz is not visible (Morlet wavelets)
% cfg.zlim           = zlim; 
% cfg.fontsize       = 14
% cfg.title          = 'Human Body Scramble'
% cfg.colormap       = '-RdGy';
% cfg.channel        = 'all';
% 
% 
% fig1 = figure;
% ft_singleplotTFR(cfg, FGA_hum_body_scr);
% 
% % % save
% % saveas(fig1,'Fig2_allchan_scramble.png');



zlim = [-.5 .5];



% human body normal - human body scramble 
cfg = [];
cfg.parameter       = 'powspctrm';
cfg.operation       = '(x1-x2)';

tfr_difference = ft_math(cfg, FGA_hum_body_norm, FGA_hum_body_scr);

cfg = [];
cfg.xlim           = [-0.1 1] % before -0.1 s is not visible (Morlet wavelets)
cfg.ylim           = [4 30] % below 4 Hz is not visible (Morlet wavelets)
cfg.zlim           = zlim;
cfg.fontsize       = 14
cfg.title          = 'Human body normal - Human body scramble'
cfg.colormap       = '-RdGy';
cfg.channel        = 'all';

fig1 = figure;
ft_singleplotTFR(cfg, tfr_difference);

% % save
% saveas(fig1,'Fig2_allchan_NS_diff.png');



% monkey body normal - monkey body scramble 
cfg = [];
cfg.parameter       = 'powspctrm';
cfg.operation       = '(x1-x2)';

tfr_difference = ft_math(cfg, FGA_monk_body_norm, FGA_monk_body_scr);

cfg = [];
cfg.xlim           = [-0.1 1] % before -0.1 s is not visible (Morlet wavelets)
cfg.ylim           = [4 30] % below 4 Hz is not visible (Morlet wavelets)
cfg.zlim           = zlim;
cfg.fontsize       = 14
cfg.title          = 'Monkey body normal - Monkey body scramble'
cfg.colormap       = '-RdGy';
cfg.channel        = 'all';

fig1 = figure;
ft_singleplotTFR(cfg, tfr_difference);

% % save
% saveas(fig1,'Fig2_allchan_NS_diff.png');



% human face normal - human face scramble 
cfg = [];
cfg.parameter       = 'powspctrm';
cfg.operation       = '(x1-x2)';

tfr_difference = ft_math(cfg, FGA_hum_face_norm, FGA_hum_face_scr);

cfg = [];
cfg.xlim           = [-0.1 1] % before -0.1 s is not visible (Morlet wavelets)
cfg.ylim           = [4 30] % below 4 Hz is not visible (Morlet wavelets)
cfg.zlim           = zlim;
cfg.fontsize       = 14
cfg.title          = 'Human face normal - Human face scramble'
cfg.colormap       = '-RdGy';
cfg.channel        = 'all';

fig1 = figure;
ft_singleplotTFR(cfg, tfr_difference);

% % save
% saveas(fig1,'Fig2_allchan_NS_diff.png');



% monkey face normal - monkey face scramble 
cfg = [];
cfg.parameter       = 'powspctrm';
cfg.operation       = '(x1-x2)';

tfr_difference = ft_math(cfg, FGA_monk_face_norm, FGA_monk_face_scr);

cfg = [];
cfg.xlim           = [-0.1 1] % before -0.1 s is not visible (Morlet wavelets)
cfg.ylim           = [4 30] % below 4 Hz is not visible (Morlet wavelets)
cfg.zlim           = zlim;
cfg.fontsize       = 14
cfg.title          = 'Monkey face normal - Monkey face scramble'
cfg.colormap       = '-RdGy';
cfg.channel        = 'all';

fig1 = figure;
ft_singleplotTFR(cfg, tfr_difference);

% % save
% saveas(fig1,'Fig2_allchan_NS_diff.png');


% human object normal - human object scramble 
cfg = [];
cfg.parameter       = 'powspctrm';
cfg.operation       = '(x1-x2)';

tfr_difference = ft_math(cfg, FGA_hum_obj_norm, FGA_hum_obj_scr);

cfg = [];
cfg.xlim           = [-0.1 1] % before -0.1 s is not visible (Morlet wavelets)
cfg.ylim           = [4 30] % below 4 Hz is not visible (Morlet wavelets)
cfg.zlim           = zlim;
cfg.fontsize       = 14
cfg.title          = 'Human object normal - Human object scramble'
cfg.colormap       = '-RdGy';
cfg.channel        = 'all';

fig1 = figure;
ft_singleplotTFR(cfg, tfr_difference);

% % save
% saveas(fig1,'Fig2_allchan_NS_diff.png');



% monkey object normal - monkey object scramble 
cfg = [];
cfg.parameter       = 'powspctrm';
cfg.operation       = '(x1-x2)';

tfr_difference = ft_math(cfg, FGA_monk_obj_norm, FGA_monk_obj_scr);

cfg = [];
cfg.xlim           = [-0.1 1] % before -0.1 s is not visible (Morlet wavelets)
cfg.ylim           = [4 30] % below 4 Hz is not visible (Morlet wavelets)
cfg.zlim           = zlim;
cfg.fontsize       = 14
cfg.title          = 'Monkey object normal - Monkey object scramble'
cfg.colormap       = '-RdGy';
cfg.channel        = 'all';

fig1 = figure;
ft_singleplotTFR(cfg, tfr_difference);

% % save
% saveas(fig1,'Fig2_allchan_NS_diff.png');