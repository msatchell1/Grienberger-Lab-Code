function [S_all] = zero_Ca_decay(S_all)
% Adds a new array to S_all that is the same as S_all.datasetSm, but with
% transient decaying Ca activity points set to zero. Decaying dF/F data
% points are determined by analyzing dF/F slope and magnitude. New array is
% called S_all.datasetSm_nodecay.

dFF_data = S_all.datasetSm;

S_all.datasetSm_nodecay = zeros(size(dFF_data));

dFF_diffs = diff(dFF_data, 1, 1);

% Plot histogram showing the distribution of differences between adjacent
% points (see documentation on diff()).
% figure;
% H = histogram(dFF_diffs,500);
% title('Histogram of Slope Between Adjacent Points For All Neurons')
% ylabel('Counts')
% xlabel('Slope')


for i = 1:size(dFF_diffs,2) % Loops each neuron.

    % Makes 1-D logical array for points where the slope is nonnegative and magnitude is above
    % threshold. Specifically, for any pair of points x(1) and x(2), the
    % difference between the two is taken as the slope which must be less
    % than zero, and the point x(1) is then evaluated as to whether it is over
    % three times the standard deviation of noise. If both are true, x(1) is
    % marked for zeroing with a logical true (1).
    diffs_log = dFF_diffs(:,i) < 0 & dFF_data(1:end-1,i) >= S_all.dFF_noise_std(i).*3;

    % Add a zero to the end of the array so it is the same size as
    % datasetSm.
    diffs_log(end+1) = 0;


    
    max_group_decay = 3; % Maximum group size for "smoothing out" Ca decay 
    % identifier array points.
    max_group_rise = 1; % Maximum group size for "smoothing out" Ca rise
    % identifier array points. Note that smoothing out rise data may not be
    % a good idea at all, and needs further investigation.

    for j = 1:(length(diffs_log) - (max_group_decay+1)) % Loops each frame
        % up until the edge case limit.
        

        % There are many instances in a Ca spike where the slope briefly
        % flattens out or goes positive in the decaying signal which we would
        % like to remove. To do this, groups of max_group points or fewer that are marked with logical 0 but
        % have 1s on the boundaries are turned to 1s.
        if diffs_log(j) == 1 % If on a decaying portion of Ca signal.
            
            % If the next point is marked as "not decaying singal" with 0 and the 
            % group size of 0s is no larger than max_group.
            if diffs_log(j+1) == 0 && any(diffs_log(j+1:j+max_group_decay+1) == 1)
                
                % Finds the first instance where the logical array again
                % has a 1. (note edge_ind only uses indices between 1 and
                % max_group). 
                edge_ind = find(diffs_log(j+1:j+max_group_decay+1) == 1, 1, "first");
                    
                diffs_log(j+1:j+edge_ind) = 1; % Turns all points inbetween to 1.

            end
        
        end
        
        
        % The same thing needs to be done for rising Ca signal. Below is
        % performed the inverse of above - groups of 1s size max_group and
        % below bordered by 0s are converted to 0s.
        if diffs_log(j) == 0 % If on a rising portion of Ca signal.
            
            % If the next point is marked as "decaying singal" with 1 and the 
            % group size of 1s is no larger than max_group.
            if diffs_log(j+1) == 1 && any(diffs_log(j+1:j+max_group_rise+1) == 0)
                
                % Finds the first instance where the logical array again
                % has a 0. (note edge_ind only uses indices between 1 and
                % max_group). 
                edge_ind = find(diffs_log(j+1:j+max_group_rise+1) == 0, 1, "first");
                    
                diffs_log(j+1:j+edge_ind) = 0; % Turns all points inbetween to 0.

            end
        
        end
    


    end
    
    % Sets value of all indices where diffs_log == 1 to 0.
    dFF_zeroed = dFF_data(:,i);
    dFF_zeroed(find(diffs_log)) = 0;

    S_all.datasetSm_nodecay(:,i) = dFF_zeroed;
    
end



end