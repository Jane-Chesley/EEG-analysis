% FDR correction 
clear;clc; close all;

% input original p-values 
orig_p = [0.016; 0.072;0.060];

% output FDR corrected p-values 'adj_p'
% 'pdep' (default) = independent tests; 'dep' = dependent tests 
% 'yes' or 'no' (default) = provide summary report to command line 
% [h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(orig_p,0.005);
q = 0.01;
[h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(orig_p,q,'pdep','yes');

disp(orig_p)
disp(adj_p)
% % adj_p < 0.05 = sig. after FDR correction 

% % % save output 
% T = array2table(adj_p);
% fileName = strcat('FDR_output_chanxchan_1to33.xlsx');
% writetable(T,fileName)
