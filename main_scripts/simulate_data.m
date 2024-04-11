
% create simulated data for one subject 
% sim_data.trial{trial}(channel,timepts)
sim_data = data_clean_hum_body_norm{1};



for tr = 1:length(sim_data.trial) % loop for each trial

    dim = size(sim_data.trial{tr}(1,:)); % get dimensions of time points for one trial and one channel
    
    % simulate random noise A (random values between 0 and 1)
    rand_A = rand(dim); % create a new vector of random noise, with the same length as the original time points 
    sim_data.trial{tr}(1,:) = rand_A; % assign it to channel 1 
    
    % simulate random noise B 
    rand_B = rand(dim); % create a new vector of random noise; this random noise is different from channel 1
    sim_data.trial{tr}(2,:) = rand_B; % assign it to channel 2 
    
    % simulate sinusoid signal A 
    amplitude = 1; % we are not manipulating amplitude, but it is included to understand the equation
    frequency = 5; % theta 
    timepts = sim_data.time{tr}; % get time points (in seconds) for each trial (575 time points per trial)
    sin_A = amplitude*sin(2*pi*frequency*timepts); % construct sinusoid signal 
    sim_data.trial{tr}(3,:) = sin_A; % assign it to channel 3 

    % simulate sinusoid signal B; similar phase to sin A (but not identical because PLI attenuates this)
    phase_difference = pi/128;  % Phase difference between the two signals
    sin_B = amplitude*sin(2*pi*frequency*timepts + phase_difference); % construct sinusoid signal with a phase difference  
    sim_data.trial{tr}(4,:) = sin_B; % assign it to channel 4
    
    % simulate sinusoid signal A 
    amplitude = 1; % we are not manipulating amplitude, but it is included to understand the equation
    frequency = 10; % alpha 
    timepts = sim_data.time{tr}; % get time points (in seconds) for each trial (575 time points per trial)
    sin_C = amplitude*sin(2*pi*frequency*timepts); % construct sinusoid signal 
    sim_data.trial{tr}(5,:) = sin_C; % assign it to channel 5

    % simulate sinusoid signal B; similar phase to sin A (but not identical because PLI attenuates this)
    phase_difference = pi/128;  % Phase difference between the two signals
    sin_D = amplitude*sin(2*pi*frequency*timepts + phase_difference); % construct sinusoid signal with a phase difference  
    sim_data.trial{tr}(6,:) = sin_D; % assign it to channel 6

end

% channel 1 = random noise A (values between 0 and 1)
% channel 2 = random noise B (values between 0 and 1)
% channel 3 = sinusoid signal A (theta)
% channel 4 = sinusoid signal A (theta)
% channel 5 = sinusoid signal B (alpha)
% channel 6 = sinusoid signal B (alpha)

%%%%%%%%% compute PLI on simulated data 

% expectations of simulated results:
% channel pairs 1-2 (noise A - noise B) not coherent (PLI = 0)
% channel pairs 1-3 (noise - sin) not coherent (PLI = 0)
% channel pairs 3-4 (sin A - sin B) close to perfectly coherent in theta band (PLI = 1)
% channel pairs 5-6 (sin C - sin D) close to perfectly coherent in alpha band (PLI = 1) 

sim_PLI = PLI_allfreq; 

roi_idx = [1;2];
sim_PLI_pair1 =  0.1456; 

roi_idx = [1;3]; 
sim_PLI_pair2 = 0.1220; 

roi_idx = [3;4]; 
sim_PLI_pair3 = 1; 

roi_idx = [5;6]; 
sim_PLI_pair4 = 1; 



