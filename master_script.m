% Master Script for behavioral and Ca imaging data analysis.

% The variables in S (now S_all) are described below:
%
% S.Running - Animal speed (V) (multiply by 80 to get cm/s).
%
% 
% S.acq is the data acquisition rate
% S.wsALL(:,1) = lap start (V)
% S.wsALL(:,2) = running including stationary periods (V) (multiply by 80 to get to cm/s)
% S.wsALL(:,3) = valve opening (V)
% S.wsALL(:,4) = licking (V)
% S.wsALL(:,5) = frame clock (V)
% S.wsALL(:,6) = position (V) integrated signal from running
%
% S.wsL - Same information as wsALL, just divided up by lap.
%
% S.dataset - Fl values in a 1D column, one for each cell. Recorded for a given
% number of frames (72000 frames in the case of CG004_220307).
% S.datasetSm - smoothed fl values, same structure as S.dataset.
% S.dataperlap - Fl values seperated up by lap and cell. This array is a
% cell array, where the first indices indicate the cell number, and the
% second indices the lap number. Fl measurements are still recorded for
% every frame. Because laps take a variable amount of time to complete,
% laps have a varying number of frames within them.
%
% S.R - mean fluorescence per bin. Each lap has 50 bins.
% S.Rmax - max fluorescence per bin.
%
% S.Fzero - The baseline fluorescence (fl). Gaussian fitting is done on a
% histogram of fluorescence values over time, and since cells are most
% often at their baseline fl, there will be a large peak at the baseline fl
% value. So gaussian fitting to this finds Fzero for each cell. This is
% then needed for the dF/F plots, because in reality the equation is: 
% (F - Fzero)/Fzero. Thus dF = F - Fzero, and the F in dF is actually
% Fzero.
%
%

% Directory where both behavioral and fluorescence data is located.
data_dir = 'D:\Michael Satchell Rotation\CG004_220307'; % Must be a text scalar; use '' instead of "". 
% The directory where data from Matlab processing will be saved.
save_dir = [data_dir, filesep, 'Matlab Data'];


if not(isfolder(save_dir)) % If there is no folder at save_dir, run 
    % wavesurfer and suite2p processing functions.

    S_ws = struct; % Creates empty scalar structure for behavioral data.
    
    [S_ws, control] = BRim0_MS(S_ws, data_dir); % Add behavioral data to S_ws.
    
    [S_all, control] = BRim1_MS(S_ws, control); % Create new structure that holds both 
    % wavesurfer and suite2p data.
    
    
    mkdir(save_dir) % Create a folder to hold the processed data.
    
    
    save([save_dir, filesep, 'S_all.mat'], 'S_all', '-v7.3') % Save ws and fl data.
    % Necessary to use version 7.3 to save data larger than 2GB.
    save([save_dir, filesep, 'control.mat'], 'control', '-v7.3') % Save lap control data.
    
else 
    % If file already exists, load existing data.
    load([save_dir, filesep, 'S_all.mat'])
    load([save_dir, filesep, 'control.mat'])

end


%%
num_nrns = size(S_all.dataset,2); % Number of neuron ROIs
num_ROIs = size(S_all.datasetorig,2); % Number of total ROIs (neurons + neuropil)
num_neu = num_ROIs - num_nrns;
num_laps = size(S_all.Running,2); % Number of laps the mouse has run
num_bins = size(S_all.Running,1); % Number of bins each lap is tiled into.

S_all = noise_dist(S_all, 0); % Calculate standard deviations for the noise 
% of the dF/F traces of each neuron. 



%% Save data

save([save_dir, filesep, 'S_all.mat'], 'S_all', '-v7.3') % Save ws and fl data.
% Necessary to use version 7.3 to save data larger than 2GB.
save([save_dir, filesep, 'control.mat'], 'control', '-v7.3') % Save lap control data.
