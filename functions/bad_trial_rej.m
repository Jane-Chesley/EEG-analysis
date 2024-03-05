% this function rejects "bad" trials exceeding 3 SD of the variance per channel 

function [data_clean, bad_trials] = bad_trial_rej(data_ICA, subjectid)

% for a given subject, calculate the variance separately for each trial and channel
% variance = SD^2; distance from mean
for t = 1:length(data_ICA.trial) % t = trial 
    for c = 1:length(data_ICA.label) % c = channel 
        data = data_ICA.trial{t}(c,:); % get 1 data point for each trial-channel pair --> 1 data point = 1x500 time-samples for a given trial-channel pair (250 samples/second * 2s per trial including pre- and post- stim periods)
        data_variance(c,t) = var(data); % calculate the variance for each data point --> 33 channels x 480 trials 
    end 
end

% calculate the average variance of all trials (dim 2) for each channel 
% data_variance_mean = 33 channels x 1 mean variance value
data_variance_mean = mean(data_variance,2);   

% calculate the standard deviation of the variance of all trials (dim 2) for each channel
% data_variance_SDs = 33 channels x 1 SD variance value
data_variance_SDs = std(data_variance, 0, 2); % 0 = weight of 0 --> normalize by N-1 VS. 1 = weight of 1 --> normalize by N  

% identify outlier trials 
for c = 1:length(data_ICA.label) % c = channel
    crit        = 3*data_variance_SDs(c); % for each channel, define a critical cutoff (crit) for outlier trials; here, 3 SDs
    thr         = data_variance_mean(c)+crit; % for each channel, calculate 3 SDs above the mean variance (thr) % all variance values are positive; so no need to calculate 3 SDs below
    outl(c,:)   = data_variance(c,:)>thr; % for each channel, find trials with a variance exceeding the threshold ( == logical array)
    idx{c}      = find(outl(c,:)==1); % get indices of those trials and store in cell array (each cell = 1 channel)
end

% identify non-outlier trials to keep in data  
outlier_trials = [idx{:}]; % concatenate all indices of outlier trials across all channels 
outlier_trials = unique(outlier_trials);  % get only unique values (if trial N is an outlier for one channel, then it's an outlier trial for all channels)
all_trials = 1:length(data_ICA.trial); 
keep_trials = setdiff(all_trials,outlier_trials); % setdiff returns values in A that are not in B

% redefine trials in input data (include only non-outlier trials)
cfg = []; 
cfg.trials = keep_trials; 
data_clean = ft_redefinetrial(cfg,data_ICA); 

% record rejection info in a structure 
bad_trials.subjectid = subjectid; 
bad_trials.outlier_trials = outlier_trials; 
bad_trials.total_preserved = length(keep_trials); 
bad_trials.total_trials = length(all_trials); 
bad_trials.perc_preserved = length(keep_trials)/length(all_trials)*100; 