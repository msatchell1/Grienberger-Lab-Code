function [S_all] = zero_Ca_decay(S_all)
% Adds a new array to S_all that is the same as S_all.datasetSm, but with
% transient decaying Ca activity points set to zero. Decaying dF/F data
% points are determined by analyzing dF/F slope and magnitude.

S_all.datasetSm_nodecay = zeros(size(S_all.datasetSm));

dFF_diffs = diff(S_all.datasetSm);

% Plot histogram showing the distribution of differences between adjacent
% points (see documentation on diff()). 
figure;
H = histogram(dFF_diffs,500);

% Makes logical array for points where the slope is 

end