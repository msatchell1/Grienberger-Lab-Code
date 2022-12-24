function [S_all] = zero_Ca_decay(S_all)
% Adds a new array to S_all that is the same as S_all.datasetSm, but with
% transient decaying Ca activity points set to zero. Decaying dF/F data
% points are determined by analyzing dF/F slope and magnitude.

dFF_data = S_all.datasetSm;

S_all.datasetSm_nodecay = zeros(size(dFF_data));

dFF_diffs = diff(dFF_data, 1, 1);

% Plot histogram showing the distribution of differences between adjacent
% points (see documentation on diff()).
figure;
H = histogram(dFF_diffs,500);
title('Histogram of Slopes Between Adjacent Points')
ylabel('Counts')
xlabel('Slope')

for i = 1:size(dFF_diffs,2) % Loops each neuron.

    % Makes logical array for points where the slope is nonnegative and magnitude is above
    % threshold. Specifically, for any pair of points x(1) and x(2), the
    % difference between the two is taken as the slope which must be less
    % than zero, and the point x(1) is then evaluated as to whether it is over
    % three times the standard deviation of noise. If both are true, x(1) is
    % marked as a logical true (1).
    diffs_log = dFF_diffs(:,i) < 0 & dFF_data(1:end-1,i) >= S_all.dFF_noise_std(i).*3;

    % Add a zero to the end of the array so it is the same size as
    % datasetSm.
    diffs_log(end+1) = 0;
    
    % Sets value of all these points to 0.
    dFF_zeroed = dFF_data(:,i);
    dFF_zeroed(find(diffs_log)) = 0;

    S_all.datasetSm_nodecay(:,i) = dFF_zeroed;
    
end



end