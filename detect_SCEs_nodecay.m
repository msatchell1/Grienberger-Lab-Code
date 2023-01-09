function [S_all] = detect_SCEs_nodecay(S_all)
% Performs the sametask as detect_SCEs.m, but uses S_all.datasetSm_nodecay
% for labelling neuronal Ca activity. 
%
% Detects synchronous calcium events (SCEs). The method of detection is to
% monitor the number of cells with dF/F above
% a threshold act_thres determined by the standard deviation of that cell's noise. If
% at any time the number of active neurons (NAN) reaches the value SCE_thres, an SCE
% is said to be occuring until the number of cells drops back down. To
% prevent brief changes from splitting up a SCE, if an SCE is dropped out
% of and then re-entered within a time(frames) window of
% min_SCE_gap, the SCEs will be joined by including the brief non-SCE
% period into the SCE.
%
% This function adds which neurons have dF/F above threshold per frame
% (act_nrns_perframe_nd), the number of these neurons per frame
% (NAN_perframe_nd), and whether an SCE is occuring per frame
% (isSCE_nd) to S_all.
% 
%
% ----- Michael Satchell 12/14/22 -----

% Use the smoothed dF/F data with no Ca decay period:
dFF_data = S_all.datasetSm_nodecay;
num_frames = size(dFF_data,1);
num_nrns = size(dFF_data,2);


noise_stds = S_all.dFF_noise_std; % Standard deviation of noise dF/F.
act_thres = noise_stds.*3; % Setting the activity thresholds.

S_all.act_nrns_perframe_nd = cell(num_frames,1); % A cell array for each frame,
% where each cell holds the indices of neurons (as ordered in datasetSm) that
% have dF/F above threshold at that time.

NAN_perframe_nd = zeros(num_frames,1);

for i = 1:num_frames
    
    % Find the active neurons at frame i by comparing the activity of all neurons
    % with their thresholds simultaneously. act_nrns is an array of 1s
    % and 0s: 1s at indices of active neurons and 0 at inactives.
    act_nrns = dFF_data(i,:) > act_thres;
    act_nrn_inds = find(act_nrns); % Indices of active neurons
    
    % Add the indices of these neurons to S_all.
    S_all.act_nrns_perframe_nd{i} = act_nrn_inds;

    % Records the number of active cells at each frame.
    NAN_perframe_nd(i) = length(act_nrn_inds);

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

    if NAN_perframe_nd(i) == 0
        
        % Whenever a value of zero is come upon, assign it the avg of the
        % last four values. This works because the loop runs
        % progressively through the data.
        NAN_perframe_nd(i) = mean(NAN_perframe_nd(i-4:i-1));

    end
end

% plot(num_act_nrns_perframe, 'DisplayName', 'PMT-artifact removed')

% Save the number of active neurons per frame to S_all.
S_all.num_act_nrns_perframe_nd = NAN_perframe_nd;

%% Detecting SCEs
% Now a threshold is defined for the number of neurons required to be
% active simultaneously in order to initiate an SCE. 

% To define the threshold, the NAN trace is shuffled 20 times and the 95th
% percentile of the shuffle is taken as the threshold value. This follows
% the methods of Liu, Xin et al. (2022). “E-Cannula reveals anatomical 
% diversity in sharp-wave ripples and a driver for the recruitment of
% distinct hippocampal assemblies".

num_shuffles = 20; %The number of times to shuffle the data

NAN_shuffles = zeros([length(NAN_perframe_nd), num_shuffles]);

for i = 1:num_shuffles
    NAN_shuffles(:,i) = NAN_perframe_nd(randperm(length(NAN_perframe_nd)));
end

avg_NAN_shuffle = mean(NAN_shuffles,2);

SCE_thres = prctile(avg_NAN_shuffle, 95) % SCE detection threshold as 95th
% percentile of the avg shuffled data.

% figure;
% hold on;
% plot(S_all.frametimingOriginal./S_all.acq, avg_NAN_shuffle)
% title("Avg Shuffled NAN")
% ylabel("Active Neuron Count")
% xlabel("Time (s)")

% Unfortunately there is a problem with this method, as it defines
% thresholds that are too low.



% Another way to define SCE_thres is to use histograms, similar to what is
% done to find the noise std dF/F for each neuron. Here, only NAN during non-running
% is considered - the noise in NAN is found during standing times by
% binning all standing time NAN data into a histogram. Then the SCE_thres
% can be defined at 2 or 3 stds above the noise. Note this threshold will
% only work when applied to standing periods, so I can't use this threshold
% everywhere. Also I will need to create the threshold using data that has
% enough standing in it. 

vel = S_all.wsALL(:,2).*80; % Mouse velocity in cm/s.

stand_thres = 1; % Everything below this value qualifies as standing (cm/s).
stand_log = vel < stand_thres; % Logical array of frames where v < stand_thres cm/s.
stand_inds = find(stand_log); % Indices where v < stand_thres.

% The percentage of frames spent standing
prct_stand = 100*(length(stand_inds)/length(vel))

% Each imaging frame is mapped onto a wavesurfer frame with
% S_all.frametimingOriginal. Below takes the knowledge of which ws frames
% occur when vel < stand_thres, and uses frametimingOriginal to convert
% that to which imaging frames occur when vel < stand_thres.

 % Imaging frames that occur during standing in [ws frame timing, imaging
 % frame timing].
[stand_im_frames_ws, stand_im_frames_im, ~] = intersect(S_all.frametimingOriginal, stand_inds, 'stable');

stand_NAN = NAN_perframe_nd(stand_im_frames_im); % Selects the NAN values that 
% occur during standing imaging frames.

ws_time_sec = (1:size(S_all.wsALL,1))./S_all.acq; % Gets ws frame times in seconds.

figure;
hold on;
yyaxis left
plot(ws_time_sec, S_all.wsALL(:,2).*80,'DisplayName', 'velocity')
ylabel('Velocity (cm/s)')

yyaxis right
plot(stand_im_frames_ws./S_all.acq, stand_NAN, 'DisplayName', 'num active nrns')
ylabel("Number of Active Neurons")

% legend();
title("NAN During Standing Only and Velocity, No Ca Decay")
xlabel("Time (s)")


% Now the NAN values during standing are binned onto a histogram and a
% gauss fit to the noise part of the signal is used to define a std.






min_SCE_gap = 50; % The minimum gap allowed between SCEs (imaging frames).

isSCE_nd = zeros(num_frames,1); % Array to hold whether or not an SCE is occuring 
% at each time frame.

frames_since_SCE = 0; % Counter to record the number of frames between SCEs.
for i = 1:num_frames
    
    if NAN_perframe_nd(i) >= SCE_thres % SCE threshold reached.

        if frames_since_SCE < min_SCE_gap % If an SCE occured recently, designates to
            % previous SCE.
            isSCE_nd(i-frames_since_SCE:i) = 1; % Sets all values since last SCE to 1.
        else % Designates to new SCE.
            isSCE_nd(i) = 1;
        end
        
        frames_since_SCE = 0; % Reset frame counter.

    end
    
    frames_since_SCE = frames_since_SCE + 1;
end

S_all.isSCE_nd = isSCE_nd;


figure;
hold on;
plot(S_all.frametimingOriginal./S_all.acq, NAN_perframe_nd, 'DisplayName', '# nrns')
plot(S_all.frametimingOriginal./S_all.acq, zeros(size(NAN_perframe_nd))+SCE_thres, 'DisplayName', 'SCE threshold')
plot(S_all.frametimingOriginal./S_all.acq, isSCE_nd.*35, 'DisplayName', 'isSCE')
title('SCE Detection')
ylabel('Number of Active Neurons')
xlabel('Time (s)')
legend();

end