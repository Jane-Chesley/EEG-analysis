% this function rejects "bad" trials exceeding 3 SD of the variance per channel 

function [data_clean, bad_trials] = bad_trial_rej(data_ICA, subjectid)

% for a given subject, calculate the variance separately for each trial and channel
% variance = SD^2; distance from mean
for t = 1:length(data_ICA.trial) % t = trial 
    for c = 1:length(data_ICA.label) % c = channel 
        data = data_ICA.trial{t}(c,:); % get 1 data point for each trial i and channel i 
        % 1 data point = 1x575 time-samples for a given trial-channel pair 
        % here, sampling rate = 250 Hz 
        % 250 samples/second * 2.3s (0.6s pre-stim + 1s stim + 0.7s post-stim) = 575 time-samples
        data_variance(c,t) = var(data); % calculate the variance for each data point 
        % results in data points for 33 channels x 480 trials 
    end 
end



% calculate the average variance of all trials (dim 2) for each channel 
% data_variance_mean = 33 channels x 1 mean variance value
data_variance_mean = mean(data_variance,2);   



% calculate the standard deviation of the variance of all trials (dim 2) for each channel
% data_variance_SDs = 33 channels x 1 SD variance value
data_variance_SDs = std(data_variance, [], 2); % [] = weight for normalizing SD; 2 = dim 2 (trials) 



% define the critical cutoff for outliers in # of SDs (3 is typical)
crit = 3; 
clear outl
% identify outlier trials 
for c = 1:length(data_ICA.label) % c = channel
    thr         = data_variance_mean(c)+crit*data_variance_SDs(c); % for each channel, calculate 3 SDs above the mean variance (thr) % all variance values are positive; so no need to calculate 3 SDs below
    outl(c,:)   = data_variance(c,:)>thr; % for each channel, find trials with a variance exceeding the threshold ( == logical array)
    idx{c}      = find(outl(c,:)==1); % get indices of those trials and store in cell array (each cell = 1 channel)
end



% identify non-outlier trials to keep in data  
outlier_trials = [idx{:}]; % concatenate indices of outlier trials across all channels 
outlier_trials = unique(outlier_trials);  % get only unique values (e.g. if trial N is an outlier for one channel, then it's an outlier trial for all channels)
all_trials = 1:length(data_ICA.trial); 
keep_trials = setdiff(all_trials,outlier_trials); % return trial indices that are NOT outliers ( setdiff(A,B) returns values of in A that are not in B )



% redefine trials in data_clean to include only non-outlier trials
cfg = []; 
cfg.trials = keep_trials; 
data_clean = ft_redefinetrial(cfg,data_ICA); 



% record rejection info in a structure 
bad_trials.subjectid = subjectid; 
bad_trials.outlier_trials = outlier_trials; 
bad_trials.total_preserved = length(keep_trials); 
bad_trials.total_trials = length(all_trials); 
bad_trials.perc_preserved = length(keep_trials)/length(all_trials)*100; 