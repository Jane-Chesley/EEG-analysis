%% Connectivity visualization 

%% ------------------------------------------------------------------------
%  Part 1 - construct .node and .edge files for BrainNet Viewer visualization  
%  ------------------------------------------------------------------------

% https://www.fieldtriptoolbox.org/template/electrode/
% https://github.com/fieldtrip/fieldtrip/blob/master/template/electrode/standard_1020.elc
 

standard_montage = ft_read_sens('standard_1020.elc');

all_channels = standard_montage.label;
channels_preserved = {'C3', 'CP3', 'P3'};

% easycapM3: all 33 channels are as follows:
% {'Fp1','F7','Fp2','F8','F3','Fz','F4','AFz','FT7','FC3','FCz','FC4','FT8','T8','C4','Cz','C3','T7','TP9','TP7','CP3','CPz','CP4','TP8','TP10','P8','P4','Pz','P3','P7','O1','Oz','O2'}; 

% Find the indices of the specified channels
idx = find(ismember(all_channels, channels_preserved)); 

% get MNI coordinates of the specified channels 
xyz_coordinates = standard_montage.elecpos(idx,:); 

% .node file requirements:
% rows = channels 
% col 1-3 = x,y,z coordinates, respectively 
% col 4 = node colors (E.g., a number 1-9, different colors for different subregions)
% col 5 = node sizes (E.g., nodal degree, centrality measure, t-value, etc.)
% col 6 = node labels (channel labels, '-' represents no label')
node_colors = ones(length(node_labels),1)
node_sizes = ones(length(node_labels),1)

% Create a table
T = table(xyz_coordinates(:,1), xyz_coordinates(:,2), xyz_coordinates(:,3), node_colors, node_sizes,channels_preserved','VariableNames', {'X', 'Y', 'Z', 'Color', 'Size','Label'});


% Save table 
writetable(T, 'connectivity.txt', 'Delimiter', ' ', 'WriteVariableNames', false); % Use space as delimiter


% manually change connectivity.txt to connectivity.node 

% construct .edge file
edges = ones(length(channels_preserved));
% UPDATE with matrix of PLI values 

% save table 
writetable(T,'connectivity.txt','Delimiter', ' ','WriteVariableNames', false); 
% manually change connectivity.txt to connectivity.edge


% Run Brain Net Viewer 
cd(    '/Users/jane_chesley/Documents/BrainNetViewer_20191031')

% not working with the files I made ... 




