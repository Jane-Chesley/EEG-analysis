function [PLI_allsubjects] = PLI_within(data,roi_idx,freq)

% input data are PLI values organized as follows:
% data{freq}(channel,channel,subject)

PLI_extracted_matrix = data{freq}(roi_idx, roi_idx,:);
% PLI_extracted_matrix(channel,channel,subject)

% pre-allocate var
PLI_allsubjects = zeros(size(PLI_extracted_matrix,3),1);
% PLI_allsubjects(subject,mean_PLI)

for s = 1:size(PLI_extracted_matrix,3) % loop for all subjects

    single_matrix = PLI_extracted_matrix(:,:,s); % get PLI matrix for one subject; matrix is symmetrical along diagonal, and 0s are along the diagonal
    upper_triangle_idx = ~(tril(single_matrix)); % get logical indices of the matrix, such that 0s are on the lower triangle, 1s are along the diagonal and on the upper triangle 
    lower_triangle = single_matrix(~upper_triangle_idx); % extract PLIs of one triangle (lower), excluding diagonal 
    mean_PLI = mean(lower_triangle); % compute mean PLI for one subject
    
    PLI_allsubjects(s,1) = mean_PLI; % store all subject data

end