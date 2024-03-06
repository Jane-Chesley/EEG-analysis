% this function preprocesses raw EEG data:
% 1. trial segmentation
% 2. assign event codes to each trial (condition #s 1-12)
% 3. remove practice trials (1-30)
% 4. re-reference signal to linked mastoids
% 5. downsample 

function [data_preprocessed] = prepare(eegfile, eventNums)

% trial segmentation
cfg                     = [];
cfg.dataset             = eegfile;
cfg.trialfun            = 'ft_trialfun_general'; % default
cfg.trialdef.eventtype  = 'Stimulus';
cfg.trialdef.eventvalue = ('S  1'); % 'XXXX'
cfg.trialdef.prestim    = 0.6; % beginning of trial relative to S1 onset in seconds %%% this is already negative !! %%%
cfg.trialdef.poststim   = 1.7; % end of trial relative to S1 onset in seconds
cfg = ft_definetrial(cfg);
data_preprocessed       = ft_preprocessing(cfg);

% event codes 
cfg.trl(:, 4)           = eventNums; % assign new condition codes to the cfg
cfg2                    = [];
cfg2.trl                = cfg.trl;
data_preprocessed       = ft_redefinetrial(cfg2,data_preprocessed); % re-define trials according to new condition code

% remove practice trials 1-30
cfg                     = [];
cfg.trials              = 31:510;
data_preprocessed       = ft_preprocessing(cfg, data_preprocessed);

% re-reference
% here I re-reference to the average of the signals at left (M1/LMON) and right (M2/RMOF) mastoid
% however, there are many options: E.g. re-reference to Cz, re-reference to avg of EEG channels (excl. EOG)
cfg = [];
cfg.channel             = 'all';
cfg.reref               = 'yes';
cfg.refmethod           = 'avg'; % default
cfg.implicitref         = 'M1'; % the implicit (REF plug) reference channel is added to the data representation
cfg.refchannel          = {'M1', 'RMOF'}; % the average of these two is used as the new reference
data_preprocessed       = ft_preprocessing(cfg, data_preprocessed);

% exclude REF channels
cfg= [];
cfg.channel                         = {'all', '-M1','-RMOF'};
data_preprocessed                   = ft_selectdata(cfg, data_preprocessed);

% % downsample to save space  
cfg = [];
cfg.resamplefs = 250; % must be multiple of SF
data_preprocessed = ft_resampledata(cfg, data_preprocessed);



