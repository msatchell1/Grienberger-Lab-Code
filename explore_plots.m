% Script to explore fluorescence and ws data with plots
%
% ------ Michael Satchell 12/01/22 ----------

num_nrns = size(S_all.dataset,2); % Number of neuron ROIs
num_ROIs = size(S_all.datasetorig,2); % Number of total ROIs (neurons + neuropil)
num_neu = num_ROIs - num_nrns;
num_laps = size(S_all.Running,2); % Number of laps the mouse has run
num_bins = size(S_all.Running,1); % Number of bins each lap is tiled into.

%% Show difference after BRim1_MS processing

for i = 1:1
    figure();
    hold on
    plot(S_all.frametimingOriginal./S_all.acq, S_all.datasetorig(:,i))
    title(strcat("Unaltered Fluorescence for Cell ", num2str(i)))
    ylabel("Fluorescence (V)")
    xlabel("Time (s)")
end

for i = 1:1
    figure();
    hold on
    plot(S_all.frametimingOriginal./S_all.acq, S_all.datasetSm(:,i))
    title(strcat("Bleaching/PMT Artifact Corrected Smooth dF/F for Cell ", num2str(i)))
    ylabel("dF/F")
    xlabel("Time (s)")
end

%%
% figure; 
% imagesc(transpose(S_all.Running)*80); 
% ax=gca;
% set(ax, 'XTick', [0.5 13 25.5 38 50.5]); 
% set(ax, 'XTickLabel',({'0', '45', '90','135','180'}));
% ylabel('Lap #'); 
% xlabel("Position (cm)")
% c = colorbar;
% c.Label.String = "Velocity (cm/s)";
% caxis([0 inf])
% title ('Velocity across laps')
% 
% figure();
% plot(S_all.MeanRunning*80, 'b', 'Linewidth',2); 
% ax=gca;
% set(ax, 'XTick', [1 6.25 12.5 18.75 25 31.25 37.5 43.75 50]); 
% set(ax, 'XTickLabel',({'0', '', '45', '', '90','','135','','180'}));
% xlabel('Position (cm)'); ylabel('Velocity (cm/s)'); 
% title("Mean Velocity Across All Laps")
% ylim([0 inf]);

%% Plot individual laps with only running
% figure();
% for i = 2:2
%     hold on;
%     plot(S_all.Running(:,i)*80)
%     ax=gca;
%     set(ax, 'XTick', [1 6.25 12.5 18.75 25 31.25 37.5 43.75 50]); 
%     set(ax, 'XTickLabel',({'0', '', '45', '', '90','','135','','180'}));
%     xlabel('Position (cm)'); ylabel('Velocity (cm/s)'); 
%     title("Velocity of Individual Laps, Running Only ");
%     ylim([0 inf]);
% end
% 
% % Plot individual laps with running and standing
% figure;
% for i = 2:2
%     hold on;
%     plot(S_all.wsL{1,2}{1,i}) % wsL{1,2} holds velocity data for each lap.
% %     ax=gca;
% %     set(ax, 'XTick', [0, max()]);
% %     set(ax, 'XTickLabel',({'0', '', '45', '', '90','','135','','180'}));
%     xlabel('Frames'); ylabel('Velocity (cm/s)'); 
%     title("Velocity of Individual Laps, Standing and Running ")
% end

%% Plot the mouse velocity and single cell fluorescence activity for entire recording.

ws_time_sec = (1:size(S_all.wsALL,1))./S_all.acq; % Gets ws frame times in seconds.
cellnum = 200; % Cell to plot.

vel_data = S_all.wsALL(:,2); % velocity of entire recording.
fl_data = S_all.dataset(:,cellnum); % dF/F data for that cell.
fl_data_sm = S_all.datasetSm(:,cellnum); % smoothed dF/F data for that cell.

% Include the threshold value for identifying cell activity.
act_thres = S_all.dFF_noise_std(1, cellnum)*3;

figure;
hold on;

yyaxis left
plot(ws_time_sec, vel_data.*80, 'DisplayName', 'velocity' )
ylabel("Velocity (cm/s)")
xlabel('Time (s)')

yyaxis right
plot(S_all.frametimingOriginal./S_all.acq, fl_data_sm, 'DisplayName', 'dF/F smoothed' )
% plot(S_all.frametimingOriginal./S_all.acq, fl_data, 'DisplayName', 'dF/F' )
% plot([1, ws_time_sec(end)], [act_thres, act_thres], 'DisplayName', 'activity thres')
title(strcat("Fluorescence for Cell ", num2str(cellnum), " and Mouse Velocity"))
ylabel("dF/F")
% legend()

%% Original Fluorescence plots

% for i = 3:3
%     figure();
%     hold on
%     plot(S_all.datasetorig(:,i))
%     title(strcat("Unaltered Data for Cell ", num2str(i)))
%     ylabel("Fluorescence (V)")
%     xlabel("Frames")
% %     plot(datasetNeu(:,i))
% %     plot(datasetSm(:,i)-datasetNeu(:,i))
% %     legend(["Raw F","Neuropil F","F-Fneu"])
% end


%% Plot velocity and average fluorescence trace with shaded error bars.
% I don't want to do this independently for each lap. To look for SCEs, I
% can just look at the whole fluorescence trace over the entire training
% period and compare that with the velocity.

% --------------------
% NOTE: What I did here below was wrong. I cannot simply compress the
% wavesurfer data to be the size of the site2p imaging data to align them
% to the same time scale. This is because the suite2p data was not taken
% continuously, nor did it start or end recording at the same as the ws
% data. Instead, I need to use S_all.frametimingOriginal to align the
% suite2p data to the ws data.
% -----------------------
%
% vel_exp = S_all.wsALL(:,2); % Velocity across entire experiment.
% 
% % I need to divide the velocity into chunks and take the mean of those chunks
% % so that the result is an array of equal length to the fl data.
% n_chunks = floor(length(vel_exp)/length(fl_exp));
% 
% % This line shrinks the velocity data by reshaping to a 2D matrix that is 
% % n_chunks x length(fl_data)+1 where the last column is partially filled
% % with NaNs to account for the uneven splitting into chunks. The mean of
% % this matrix is then computed along the first axis, giving a reduced size
% % mean velocity array.
% vel_exp_shk = nanmean(reshape( [vel_exp(:);nan(mod(-length(vel_exp),n_chunks),1)], n_chunks,[]));
% % vel_exp_shk(length(vel_exp_shk)-length(fl_exp):end) = []; % Delete last 
% % elements so that array has same size as fl_exp.
% 
% neg_inds = vel_exp_shk < 0; % Indices in velocity data where velocity is less than 0.
% vel_exp_shk(neg_inds) = 0; % Changes negative values to zero. Note I can no longer
% % integrate this data to accurately get position. 
% 
% vel_exp_mf = medfilt1(vel_exp_shk); % 1D median filter to smooth data.
% 
% % figure;
% % hold on;
% % plot(vel_exp_shk, 'DisplayName', 'shrunk data')
% % plot(vel_exp_mf, 'DisplayName', 'median filtered')
% % ylabel('Velocity (V)')
% % legend()


% Correct way to combine ws and suite2p data:

% Find average of fluorescence across all cells.
fl_exp = S_all.datasetSm; % Fluorescence across enitre experiment for all nrns.
fl_exp_avg = mean(fl_exp,2);
fl_exp_std = std(fl_exp,0,2);

ws_time_sec = (1:size(S_all.wsALL,1))./S_all.acq; % Gets ws frame times in seconds.

figure;
hold on;
colororder({'[0 0.4470 0.7410]','k'})

yyaxis left
plot(ws_time_sec, S_all.wsALL(:,2).*80, 'DisplayName', 'velocity');
ylabel("Velocity (cm/s)")

yyaxis right
eb_ax = shadedErrorBar(S_all.frametimingOriginal./S_all.acq, fl_exp_avg, fl_exp_std);
eb_ax.mainLine.DisplayName = 'avg dF/F';
ylabel("Average dF/F")

xlabel('Time (s)')
title('Velocity and Average dF/F')
% legend();

%% Fl of individual cells vs velocity whole experiment

% figure;
% hold on;
% plot(vel_exp_mf, 'DisplayName', 'velocity (V)');
% for i = 40:40
%     plot(S_all.datasetSm(:,i), 'DisplayName', strcat('Cell', num2str(i)));
% end
% title('Fl of Cells and Mouse Velocity - INCORRECT FRAME CONVERSION')
% xlabel('Frames')
% legend();

% % Compare dF/F with raw fl trace
% rawfl_nrns = S_all.datasetorig(:,S_all.listofneurons);
% plot(rawfl_nrns(:,i)./5000, "DisplayName", 'Cell Raw Fl')


%% SCE Events and velocity

% figure;
% hold on;
% plot(vel_exp_mf, 'DisplayName', 'velocity (V)');
% plot(S_all.isSCE.*S_all.num_act_nrns_perframe./40, 'DisplayName', 'isSCE');
% title('SCEs and Velocity (INCORRECT FRAME CONVERSION)');
% xlabel('Frame')
% legend();


%% Cell stds

r2_thres = 0.23; % Threshold for seperating the stds of cells.

n_above = length(find(S_all.dFF_noise_std >= r2_thres));
n_below = length(find(S_all.dFF_noise_std < r2_thres));
p_above = 100*n_above/size(S_all.dFF_noise_std,2);
p_below = 100*n_below/size(S_all.dFF_noise_std,2);


figure;
hold on;
plot(S_all.dFF_noise_std, 'o', 'DisplayName', 'Cell Stds')
plot(zeros(size(S_all.dFF_noise_std))+r2_thres, 'DisplayName', 'Std Threshold')
title('Cell Standard Deviations (NOT A GOOD INDICATOR OF NRN TYPE)');
ylabel('Std (dF/F)')
xlabel('Cell Number in Suite2p Order');
legend()
txt = ['above thr = ', num2str(round(p_above,0)), '%, below thr = ', num2str(round(p_below,0)), '%'];
text(0.25*length(S_all.dFF_noise_std), 0.75*max(S_all.dFF_noise_std), txt)


%% Cell r-squared values

r2_thres = 0.94; % Threshold for seperating the stds of cells.

n_above = length(find(S_all.dFF_noise_r2 >= r2_thres));
n_below = length(find(S_all.dFF_noise_r2 < r2_thres));
p_above = 100*n_above/size(S_all.dFF_noise_r2,2);
p_below = 100*n_below/size(S_all.dFF_noise_r2,2);


figure;
hold on;
plot(S_all.dFF_noise_r2, 'o', 'DisplayName', 'Cell R-sqrd')
plot(zeros(size(S_all.dFF_noise_r2))+r2_thres, 'DisplayName', 'R-sqrd Threshold')
title('Goodness of Gaussian Fits to dF/F Histograms');
ylabel('R-sqrd')
xlabel('Cell Number in Suite2p Order');
legend()
txt = ['above thr = ', num2str(round(p_above,0)), '%, below thr = ', num2str(round(p_below,0)), '%'];
text(0.25*length(S_all.dFF_noise_r2), 0.9*r2_thres, txt)



%% PLotting frameclock and velocity together

figure;
hold on;
plot(ws_time_sec, S_all.wsALL(:,5),'DisplayName', 'frameclock')
plot(ws_time_sec, S_all.wsALL(:,2), 'DisplayName', 'velocity')
title("Original ws Data for Frameclock and Velocity")
xlabel("Time (s)")
ylabel('Wavesurfer Readout (V)')
legend();

%% Plotting num active neurons with velocity

% Create array of nans the size of ws frame data.
num_act_nrns_wsf = zeros(size(S_all.wsALL(:,1)));
num_act_nrns_wsf(:) = nan; 

% Fill indices where imaging occurs with data collected during imaging. In
% this case, the number of active neurons.
num_act_nrns_wsf(S_all.frametimingOriginal) = S_all.num_act_nrns_perframe;

% figure;
% hold on;
% plot(ws_time_sec, S_all.wsALL(:,5),'DisplayName', 'frameclock')
% % plot(S_all.frametimingOriginal, num_act_nrns_wsf(~isnan(num_act_nrns_wsf)), 'DisplayName', 'num active nrns')
% plot(S_all.frametimingOriginal./S_all.acq, S_all.num_act_nrns_perframe, 'DisplayName', 'num active nrns')
% legend();
% title("Number of Active Neurons Fit Over Frameclock")
% xlabel("Time (s)")
% ylabel("Number of Active Neurons")


% And plotting the same thing but with velocity:
figure;
hold on;
yyaxis left
plot(ws_time_sec, S_all.wsALL(:,2).*80,'DisplayName', 'velocity')
ylabel('Velocity (cm/s)')

yyaxis right
plot(S_all.frametimingOriginal./S_all.acq, S_all.num_act_nrns_perframe, 'DisplayName', 'num active nrns')
ylabel("Number of Active Neurons")

% legend();
title("Number of Active Neurons with Velocity")
xlabel("Time (s)")


% 
% % And finally plotting the SCE events with velocity:
% figure;
% hold on;
% plot(ws_time_sec, S_all.wsALL(:,2),'DisplayName', 'velocity')
% plot(S_all.frametimingOriginal./S_all.acq, S_all.isSCE.*S_all.num_act_nrns_perframe./40, 'DisplayName', 'SCE events')
% legend();
% title("SCE Events with Velocity")
% xlabel("Time (s)")
% ylabel("(Number of Active Neurons)/40")

%% Number of active neurons and velocity no Ca decay


ws_time_sec = (1:size(S_all.wsALL,1))./S_all.acq; % Gets ws frame times in seconds.

figure;
hold on;
yyaxis left
plot(ws_time_sec, S_all.wsALL(:,2).*80,'DisplayName', 'velocity')
ylabel('Velocity (cm/s)')

yyaxis right
plot(S_all.frametimingOriginal./S_all.acq, S_all.num_act_nrns_perframe_nd, 'DisplayName', 'num active nrns')
ylabel("Number of Active Neurons")

% legend();
title("Number of Active Neurons and Velocity, No Ca Decay")
xlabel("Time (s)")
%% Compare normal smoothed dF/F to no decay dF/F.

cellnum = 100; % Cell to plot.

ws_time_sec = (1:size(S_all.wsALL,1))./S_all.acq; % Gets ws frame times in seconds.
fl_data_sm = S_all.datasetSm(:,cellnum); % smoothed dF/F data for that cell.

% Include the threshold value for identifying cell activity.
act_thres = S_all.dFF_noise_std(1, cellnum)*3;

figure;
hold on;
plot(S_all.frametimingOriginal./S_all.acq, fl_data_sm, 'DisplayName', 'dF/F' )
plot(S_all.frametimingOriginal./S_all.acq, S_all.datasetSm_nodecay(:,cellnum), 'DisplayName', 'dF/F no decay' )
plot([1, ws_time_sec(end)], [act_thres, act_thres], 'DisplayName', 'activity thres')
title(strcat("No-Decay Algorithm Analysis with Fluorescence from Cell ", num2str(cellnum)))
ylabel("dF/F")
xlabel('Time (s)')
legend()


%% Plot velocity and average dF/F no decay.

% Find average of fluorescence across all cells.
fl_exp = S_all.datasetSm_nodecay; % Fluorescence across enitre experiment for all nrns.
fl_exp_avg = mean(fl_exp,2);
fl_exp_std = std(fl_exp,0,2);

ws_time_sec = (1:size(S_all.wsALL,1))./S_all.acq; % Gets ws frame times in seconds.

figure;
hold on;
colororder({'[0 0.4470 0.7410]','k'})

yyaxis left
plot(ws_time_sec, S_all.wsALL(:,2).*80, 'DisplayName', 'velocity');
ylabel("Velocity (cm/s)")

yyaxis right
eb_ax = shadedErrorBar(S_all.frametimingOriginal./S_all.acq, fl_exp_avg, fl_exp_std);
eb_ax.mainLine.DisplayName = 'avg dF/F';
ylabel("Average dF/F")

xlabel('Time (s)')
title('Velocity and Avg dF/F No Ca Decay')
% legend();
