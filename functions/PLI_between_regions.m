function [PLI_allsubjects] = PLI_between_regions(data,roi1_idx,roi2_idx,freq)

PLI_extracted_matrix = data{freq}(roi1_idx,roi2_idx,:); 

for s = 1:size(PLI_extracted_matrix,3)

    single_matrix = PLI_extracted_matrix(:,:,s); 
    mean_PLI = mean(single_matrix,'all'); 

    PLI_allsubjects(s,1) = mean_PLI; % store all subject data 

end
