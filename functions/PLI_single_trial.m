
function PLI = PLI_single_trial(EEG_data)
% Input: EEG time series data for a single trial (EEG_data = channels*time_points matrix)
% Output: PLI (channels*channels matrix, AKA phase lag index matrix for a single trial)

% get total # of channels 
num_channels = size(EEG_data,1);

% Calculate PLI for all possible channel pairs
PLI = zeros(num_channels); % pre-allocate variable 

for i = 1:num_channels 
    for j = 1:num_channels

        signal1 = EEG_data(i,:); % identify one signal (the EEG time-course for a given channel)
        signal2 = EEG_data(j,:); % identify a second signal (the EEG time-course for another given channel)
        
        % Calculate phase angles of the two signals  
        % hilbert() performs hilbert transformation, such that a real signal is transformed into an analytic signal 
        % angle() extracts the phase angle from the imaginary part of the analytic signal 
        % Note: if the input to hilbert() is a matrix, then the transformation is computed column-wise
        % This is counterintuitive for the EEG data, so I prefer to evaluate the signals one by one 
        phase_angle1 = angle(hilbert(signal1));
        phase_angle2 = angle(hilbert(signal2));
        
        % Calculate the instantaneous phase differences of the two signals over time  
        % When subtracting phase angles, you cannot do simple arithmetic subtraction 
        % Instead, you can do circular subtraction or sin(difference), as below 
        phase_differences = sin(phase_angle2 - phase_angle1);

        % Calculate Phase Lag Index (PLI) of the two signals 
        % PLI quantifies the absolute value of the average SIGNs (+ or -) of the phase angle differences, 
        % which represents the average of the amount of phase leads (+) and phase lags (-) 
        % To determine the directionality of the PLI (lead or lag), simply omit abs()
        current_PLI = abs(mean(sign(phase_differences))); 
        
        % assign the PLI for the current channel-pair to the respective indices of a channel*channel matrix 
        PLI(i,j) = current_PLI; 

        clear signal1 signal2 phase_angle1 phase_angle2 phase_differences current_PLI 

    end
end
















