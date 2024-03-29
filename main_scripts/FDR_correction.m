% FDR correction 
clear;clc; close all;

% input original p-values 
% orig_p = [0.016; 0.072;0.060; 0.110; 0.371; 0.006];
orig_p = [0.001;0.270;0.009];

% output FDR corrected p-values 'adj_p'
% 'pdep' (default) = independent tests; 'dep' = independent or dependent tests 
% 'yes' or 'no' (default) = provide summary report to command line 
% [h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(orig_p,0.005);
q = 0.01;
[h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(orig_p,q,'dep','yes');

disp(orig_p)
disp(adj_p)
% % adj_p < 0.05 = sig. after FDR correction 

% % % save output 
% T = array2table(adj_p);
% fileName = strcat('FDR_output_chanxchan_1to33.xlsx');
% writetable(T,fileName)
