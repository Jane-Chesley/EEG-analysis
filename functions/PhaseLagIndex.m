function PLI = PhaseLagIndex(EEG)
% Input: EEG (EEG time series matrix), channels * time points
% Output: PLI (channels * channels, phase lag index matrix)
% author: Min Wu, 26-12-2020

[numChannels, ~] = size(EEG); 
PLI = zeros(numChannels);
h = hilbert(EEG'); % time*channels
angle_h = angle(h); % -pi <= angle <= pi, time*channels

for i = 1:numChannels-1
    for j = i+1:numChannels
        % E.g. 
        % channel 1 phase value - channel 2 phase value = phase difference for channel pair 1, calculated for each time point 
        % then take mean of phase differences for channel pair 1 for all time points  
        % output = 1 PLI for each channel pair 
        PLI(i, j) = abs(mean(sign(sin(angle_h(:, i) - angle_h(:, j))))); 
    end
end

% Symmetry Correction (transform triangle to matrix) by summing PLI and transposed PLI' 
% PLI values are calculated only for one direction (from channel i to channel j), 
% but the matrix should be symmetric  
% (E.g. PLI between channel i and j is the same as PLI between channel j and i)
PLI = PLI' + PLI; 