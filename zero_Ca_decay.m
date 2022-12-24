function [S_all] = zero_Ca_decay(S_all)
% Adds a new array to S_all that is the same as S_all.datasetSm, but with
% transient decaying Ca activity points set to zero. Decaying dF/F data
% points are determined by analyzing dF/F slope and magnitude.

S_all.datasetSm_nodecay = zeros(size(S_all.datasetSm));

dFF_diffs = diff(S_all.datasetSm, 1, 1);

% Plot histogram showing the distribution of differences between adjacent
% points (see documentation on diff()).
figure;
H = histogram(dFF_diffs,500);

% Makes logical array for points where the slope is nonnegative and magnitude is above
% threshold. Specifically, for any pair of points x(1) and x(2), the
% difference between the two is taken as the slope which must be less
% than zero, and the point x(1) is then evaluated as to whether it is over
% three times the standard deviation of noise. If both are true, x(1) is
% marked as a logical true (1).
diffs_log = dFF_diffs < 0 & S_all.datasetSm(1:end-1) >= S_all.dFF_noise_std(1:end-1).*3;






end