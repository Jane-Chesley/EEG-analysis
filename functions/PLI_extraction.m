function [PLI_allsubjects] = PLI_extraction(data,roi_idx,freq)

PLI_extracted_matrix = data{freq}(roi_idx, roi_idx,:); 

for s = 1:size(PLI_extracted_matrix,3)

    single_matrix = PLI_extracted_matrix(:,:,s);
    lower_triangle = tril(single_matrix, -1);
    lower_triangle_vector = lower_triangle(lower_triangle ~= 0); % exclude zeros 
    mean_PLI = mean(lower_triangle_vector); 

    PLI_allsubjects(s,1) = mean_PLI; % store all subject data 

end
