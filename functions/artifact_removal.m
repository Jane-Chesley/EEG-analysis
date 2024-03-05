% this function performs independent component analysis (ICA) on preprocessed EEG data
% requires manual rejection of eye and movement artifacts
% https://www.fieldtriptoolbox.org/example/ica_eog/#example-dataset

function [data_ICA, rejected_artifacts] = artifact_removal(data_preprocessed)

% visualize all channels
cfg=[];
cfg.continuous = 'yes';
cfg.blocksize = 8;
ft_databrowser(cfg,data_preprocessed);

% decompose data 
cfg = [];
cfg.trials = 'all';
cfg.continuous = 'yes'; % ICA must be done on continuous data, so force all trials to be combined
cfg.channel = 'all';
cfg.method = 'runica';
comp = ft_componentanalysis(cfg,data_preprocessed); % comp = decompose data into independent components 

% identify artifacts
% plot spatial topography of components (topoplots) and visually inspect for artifacts
cfg = [];
cfg.component = 1:length(comp.label);
cfg.layout = 'easycapM3.mat'; % specify the layout file that should be used for plotting
cfg.comment  = 'no';
figure;
ft_topoplotIC(cfg,comp)


% plot time-course  of components and visually inspect for artifacts
figure('units','normalized','outerposition',[0.2 0.2 0.8 0.8])
cfg = [];
cfg.layout = 'easycapM3.mat'; % specify the layout file that should be used for plotting
cfg.viewmode = 'component';
cfg.plotlabels = 'yes';
cfg.continuous = 'yes';
cfg.blocksize = 8;
ft_databrowser(cfg, comp);

% disp(subject);
disp('press enter to submit components to delete');
beep
pause % pause the script to be able to scroll through the trials. continue by clicking in the command window and then enter

% load a gui in which the components that should be deleted can be typed in
% it should be inserted in numeric format and separated by a comma
% identify 0,1 or 2 components that clearly correspond to eye movements based on topotplot and continuous data
deletecompstring = inputdlg({ 'components to delete'},'ICA');
deletecomp = str2double(split(deletecompstring,','))'; % NEW


% remove the components that have been entered in the gui
% backproject the data 
cfg = [];
cfg.component = deletecomp; % removed component(s) % array format [1 2 3] 
data_ICA = ft_rejectcomponent(cfg, comp, data_preprocessed);

close all


% record rejected components 
rejected_artifacts = data_ICA.cfg.component; 


% remove EOG channels from data 
cfg= [];
cfg.channel             = {'all', '-VEOG1', '-VEOG2', '-HEOG1', '-HEOG2'}; 
data_ICA                = ft_selectdata(cfg, data_ICA);


