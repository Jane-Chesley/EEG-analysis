%% Network analysis 

% Input data is ...

% 1.    Segment data by condition = 12 datasets 
%       40 trials per condition 
%       Start by working with only one condition/dataset 

%       Data is 5-D
%       subject x trial x channel x time x signal 

% 2.    Filter the signal according to your frequency of interest (for now: 0.1-30 Hz)

% 3.    Organize the data:
%       For each subject and each trial, channel x time samples matrix (values = signal)

% 4.    For one subject and one trial, 
%       Compute phase values for each channel and each time sample (Hilbert transformation) 

% 5.    Compute PLI (= one value) for each channel pair at each time point
%       (E.g. tp1 ch1-ch2; tp1 ch1-ch3; ... tp2 ch1-ch2 tp2 ch1-ch3 ...)

% ... 