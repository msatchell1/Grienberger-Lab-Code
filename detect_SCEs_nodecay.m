function [S_all] = detect_SCEs_nodecay(S_all)
% Performs the sametask as detect_SCEs.m, but uses S_all.datasetSm_nodecay
% for labelling neuronal Ca activity.
%
% Detects synchronous calcium events (SCEs). The method of detection is to
% monitor the number of cells with dF/F above
% a threshold act_thres determined by the standard deviation of that cell's noise. If
% at any time the number of active cells reaches the value SCE_thres, an SCE
% is said to be occuring until the number of cells drops back down. To
% prevent brief changes from splitting up a SCE, if an SCE is dropped out
% of and then re-entered within a time(frames) window of
% min_SCE_gap, the SCEs will be joined by including the brief non-SCE
% period into the SCE.
%
% This function adds which neurons have dF/F above threshold per frame
% (act_nrns_perframe), the number of these neurons per frame
% (num_act_nrns_perframe), and whether an SCE is occuring per frame
% (isSCE) to S_all.
% 
%
% ----- Michael Satchell 12/14/22 -----

% Use the smoothed dF/F data:
dFF_data = S_all.datasetSm_nodecay;
num_frames = size(dFF_data,1);
num_nrns = size(dFF_data,2);


noise_stds = S_all.dFF_noise_std; % Standard deviation of noise dF/F.
act_thres = noise_stds.*3; % Setting the activity thresholds.

S_all.act_nrns_perframe = cell(num_frames,1); % A cell array for each frame,
% where each cell holds the indices of neurons (as ordered in datasetSm) that
% have dF/F above threshold at that time.

num_act_nrns_perframe = zeros(num_frames,1);

for i = 1:num_frames
    
    % Find the active neurons at frame i by comparing the activity of all neurons
    % with their thresholds simultaneously. act_nrns is an array of 1s
    % and 0s: 1s at indices of active neurons and 0 at inactives.
    act_nrns = dFF_data(i,:) > act_thres;
    act_nrn_inds = find(act_nrns); % Indices of active neurons
    
    % Add the indices of these neurons to S_all.
    S_all.act_nrns_perframe{i} = act_nrn_inds;

    % Records the number of active cells at each frame.
    num_act_nrns_perframe(i) = length(act_nrn_inds);

end

% figure;
% hold on;
% plot(num_act_nrns_perframe, 'DisplayName', 'PMT-artifact included')
% title('Ensemble Activity Before and After Removing PMT Artifact')
% ylabel('Number of Active Neurons')
% xlabel('Imaging Frames')
% legend();

% The PMT-off periods drop the number of active cells to 0 about every 400
% frames (but not exactly every 400 frames, sometimes it's a frame or two off.
% To fix this, every location where the number of above-threshold nrns is 0 is
% located, and those values replaced by an average of the previous data. It
% is important to do this to avoid dropping in and out of SCEs
% unnecessarily.
for i = 1:num_frames

    if num_act_nrns_perframe(i) == 0
        
        % Whenever a value of zero is come upon, assign it the avg of the
        % last four values. This works because the loop runs
        % progressively through the data.
        num_act_nrns_perframe(i) = mean(num_act_nrns_perframe(i-4:i-1));

    end
end

% plot(num_act_nrns_perframe, 'DisplayName', 'PMT-artifact removed')

% Save the number of active neurons per frame to S_all.
S_all.num_act_nrns_perframe = num_act_nrns_perframe;

%% Detecting SCEs
% Now a threshold is defined for the number of neurons required to be
% active simultaneously in order to initiate an SCE. 

SCE_thres = 30; % SCE detection threshold (number of active neurons).
min_SCE_gap = 50; % The minimum gap allowed between SCEs (imaging frames).

isSCE = zeros(num_frames,1); % Array to hold whether or not an SCE is occuring 
% at each time frame.


frames_since_SCE = 0; % Counter to record the number of frames between SCEs.
for i = 1:num_frames
    
    if num_act_nrns_perframe(i) >= SCE_thres % SCE threshold reached.

        if frames_since_SCE < min_SCE_gap % If an SCE occured recently, designates to
            % previous SCE.
            isSCE(i-frames_since_SCE:i) = 1; % Sets all values since last SCE to 1.
        else % Designates to new SCE.
            isSCE(i) = 1;
        end
        
        frames_since_SCE = 0; % Reset frame counter.

    end
    
    frames_since_SCE = frames_since_SCE + 1;
end

S_all.isSCE = isSCE;


figure;
hold on;
plot(num_act_nrns_perframe, 'DisplayName', '# nrns')
plot(zeros(size(num_act_nrns_perframe))+SCE_thres, 'DisplayName', 'SCE threshold')
plot(isSCE.*35, 'DisplayName', 'isSCE')
title('SCE Event Detection')
ylabel('Number of Active Neurons')
xlabel('Imaging Frames')
legend();

end