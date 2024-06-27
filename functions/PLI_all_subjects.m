function [PLI_allfreq] = PLI_all_subjects(input, min_T1, max_T1)

% pre-allocate vars
channels = length(input{1}.label); % constant across subjects
subjects = length(input);
PLI_allsubjects = zeros(channels,channels,subjects); % channel*channel*subject
PLI_allfreq = cell(1,7); % six freq bands (see below)

% get indices of time-window of interest (toi)
timepts = input{1}.time{1}; % all time points, which are constant across measurements
timepts = round(timepts,2); % round time points
T1_idx = find(timepts >= min_T1 & timepts <= max_T1); % find indices of toi range

% create a log file to document script progress, as these calculations can take several hours 
logFileName = 'logfile.txt';
logfile = fopen(logFileName, 'a'); % open or create a log file for writing (append mode)
fprintf(logfile, 'Script started\n'); % log a message indicating the start of the script


for freq = 1:7



    if freq == 1 % delta
        hpfreq = 0.5;
        lpfreq = 3.9;
        bsfilter = 'no';
        bsfreq = [];

    elseif freq == 2 % theta
        hpfreq = 4;
        lpfreq = 7;
        bsfilter = 'no';
        bsfreq = [];

    elseif freq == 3 % alpha
        hpfreq = 8;
        lpfreq = 12;
        bsfilter = 'no';
        bsfreq = [];

    elseif freq == 4 % beta
        hpfreq = 13;
        lpfreq = 30;
        bsfilter = 'no';
        bsfreq = [];

    elseif freq == 5 % gamma_A
        hpfreq = 35;
        lpfreq = 45;
        bsfilter = 'no';
        bsfreq = [];

    elseif freq == 6 % gamma_B
        hpfreq = 25;
        lpfreq = 45;
        bsfilter = 'yes'; % bandstop at 30 Hz for 30 Hz white noise in stimuli 
        bsfreq = [29 31];

    elseif freq == 7 % broadband signal 1-45 Hz 
        hpfreq = 1;
        lpfreq = 45;
        bsfilter = 'no';
        bsfreq = [];

    end



    for s = 1:subjects

        % get data for a single subject
        data = input{s};

        % filter according to frequency band
        cfg                     = [];
        cfg.lpfilter            = 'yes';
        cfg.lpfreq              = lpfreq;
        cfg.hpfilter            = 'yes';
        cfg.hpfreq              = hpfreq;
        cfg.bsfilter            = bsfilter;
        cfg.bsfreq              = bsfreq;

        data                    = ft_preprocessing(cfg, data);

        % extract data corresponding to toi
        % input = data.trial{trial}(channel,timepoints)
        % output = data_extracted{trial}(channel,timepoints)
        total_trials = length(data.trial); 
        data_extracted = cell(1,total_trials); % pre-allocate var
        for trial = 1:total_trials
            data_extracted{trial} = data.trial{trial}(:,T1_idx); % T1_idx is defined above
        end

        % for each trial, assign the (channel,timepoints) matrix as input for PLI calculation
        % PLI calculation outputs a channel*channel matrix, with cells representing PLI values
        % output = PLI_output(channel,channel,trial)
        PLI_output = zeros(channels,channels,total_trials); % pre-allocate var
        for trial = 1:total_trials
            PLI_input = data_extracted{trial}; 
            PLI_output(:,:,trial) = PLI_single_trial(PLI_input);
        end

        % average across all trials (3rd dimension)
        % input = PLI_output(channel,channel,trial)
        % output = PLI_avg(channel,channel), which represents the mean PLI matrix for ONE subject and ONE condition
        PLI_avg = mean(PLI_output,3);

        % consolidate data for all subjects
        % output = PLI_allsubjects(channel,channel,subject)
        PLI_allsubjects(:,:,s) = PLI_avg;

    end


    % consolidate data for all frequencies
    PLI_allfreq{freq} = PLI_allsubjects;
    
    % update log file 
    fprintf(logfile, 'frequency completed');


end

% close the log file
fclose(logfile); 

