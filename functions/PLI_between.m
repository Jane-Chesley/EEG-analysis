function [PLI_allsubjects] = PLI_between(data,roi1_idx,roi2_idx,freq)

PLI_extracted_matrix = data{freq}(roi1_idx,roi2_idx,:); % extract relevant PLIs; this gives a 3D matrix (channels x channels x subjects)

for s = 1:size(PLI_extracted_matrix,3) % loop for all subjects 

    single_matrix = PLI_extracted_matrix(:,:,s); % extract PLIs for relevant channel-pairs for one subject 
    mean_PLI = mean(single_matrix,'all'); % compute the average; this gives one mean PLI value for each subject 

    PLI_allsubjects(s,1) = mean_PLI; % store all subject data 

end