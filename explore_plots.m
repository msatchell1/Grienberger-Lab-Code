% Script to explore fluorescence and ws data with plots
%
% ------ Michael Satchell 12/01/22 ----------

num_nrns = size(S_all.dataset,2); % Number of neuron ROIs
num_ROIs = size(S_all.datasetorig,2); % Number of total ROIs (neurons + neuropil)
num_neu = num_ROIs - num_nrns;
num_laps = size(S_all.Running,2); % Number of laps the mouse has run
num_bins = size(S_all.Running,1); % Number of bins each lap is tiled into.

%%
figure; 
imagesc(transpose(S_all.Running)*80); 
ax=gca;
set(ax, 'XTick', [0.5 13 25.5 38 50.5]); 
set(ax, 'XTickLabel',({'0', '45', '90','135','180'}));
ylabel('Lap #'); 
xlabel("Position (cm)")
c = colorbar;
c.Label.String = "Velocity (cm/s)";
caxis([0 inf])
title ('Velocity across laps')

figure();
plot(S_all.MeanRunning*80, 'b', 'Linewidth',2); 
ax=gca;
set(ax, 'XTick', [1 6.25 12.5 18.75 25 31.25 37.5 43.75 50]); 
set(ax, 'XTickLabel',({'0', '', '45', '', '90','','135','','180'}));
xlabel('Position (cm)'); ylabel('Velocity (cm/s)'); 
title("Mean Velocity Across All Laps")
ylim([0 inf]);

%% Plot individual laps with only running
figure();
for i = 2:2
    hold on;
    plot(S_all.Running(:,i)*80)
    ax=gca;
    set(ax, 'XTick', [1 6.25 12.5 18.75 25 31.25 37.5 43.75 50]); 
    set(ax, 'XTickLabel',({'0', '', '45', '', '90','','135','','180'}));
    xlabel('Position (cm)'); ylabel('Velocity (cm/s)'); 
    title("Velocity of Individual Laps, Running Only ");
    ylim([0 inf]);
end

% Plot individual laps with running and standing
figure;
for i = 2:2
    hold on;
    plot(S_all.wsL{1,2}{1,i}) % wsL{1,2} holds velocity data for each lap.
%     ax=gca;
%     set(ax, 'XTick', [0, max()]);
%     set(ax, 'XTickLabel',({'0', '', '45', '', '90','','135','','180'}));
    xlabel('Frames'); ylabel('Velocity (cm/s)'); 
    title("Velocity of Individual Laps, Standing and Running ")
end

%% Plot the mouse velocity and fluorescence activity
cellnum = 46;
lapnum = 100;

% First, because the ws data is taken at a much higher sampling frequency
% that the fl image data, I need to reduce the size of the ws data to the
% size of the fl data.
vel_data = S_all.wsL{1,2}{1,lapnum}; % velocity data by lap, high
% frequency recording.
fl_data = S_all.dataperlap{1,cellnum}{1,lapnum}; % dF/F data, lower freq recording.

% I need to divide vel_data into chunks and take the mean of those chunks
% so that the result is an array of equal length to the fl data. Thus the
% number of chunks should be floor(length(vel_data)/length(fl_data)).
n_chunks = floor(length(vel_data)/length(fl_data));

% This line shrinks the velocity data by reshaping to a 2D matrix that is 
% n_chunks x length(fl_data)+1 where the last column is partially filled
% with NaNs to account for the uneven splitting into chunks. The mean of
% this matrix is then computed along the first axis, giving a reduced size
% mean velocity array.
vel_shrunk = nanmean(reshape( [vel_data(:);nan(mod(-length(vel_data),n_chunks),1)],n_chunks,[]));

% For this plot, I want to drop the last element of vel_shrunk so that the
% data have equal lengths.
vel_shrunk(end) = [];

figure;

hold on;
plot(vel_shrunk, 'DisplayName', 'velocity' )
plot(fl_data, 'DisplayName', 'dF/F' )
title(strcat("Velocity and Fluorescence for Cell ", num2str(cellnum), " on Lap ", num2str(lapnum)))
legend()

%% Original Fluorescence plots

for i = 3:3
    figure();
    hold on
    plot(S_all.datasetorig(:,i))
    title(strcat("Unaltered Data for Cell ", num2str(i)))
    ylabel("Fluorescence (V)")
    xlabel("Frames")
%     plot(datasetNeu(:,i))
%     plot(datasetSm(:,i)-datasetNeu(:,i))
%     legend(["Raw F","Neuropil F","F-Fneu"])
end


%% Plot velocity and average fluorescence trace with shaded error bars.
% I don't want to do this independently for each lap. To look for SCEs, I
% can just look at the whole fluorescence trace over the entire training
% period and compare that with the velocity.

vel_exp = S_all.wsALL(:,2); % Velocity across entire experiment.
fl_exp = S_all.datasetSm; % Fluorescence across enitre experiment for all nrns.

% figure;
% plot(S_all.wsALL(:,2), 'DisplayName','original')
% hold on;
% plot(vel_mf, 'DisplayName','median filter')
% legend()

% I need to divide the velocity into chunks and take the mean of those chunks
% so that the result is an array of equal length to the fl data.
n_chunks = floor(length(vel_exp)/length(fl_exp));

% This line shrinks the velocity data by reshaping to a 2D matrix that is 
% n_chunks x length(fl_data)+1 where the last column is partially filled
% with NaNs to account for the uneven splitting into chunks. The mean of
% this matrix is then computed along the first axis, giving a reduced size
% mean velocity array.
vel_exp_shk = nanmean(reshape( [vel_exp(:);nan(mod(-length(vel_exp),n_chunks),1)], n_chunks,[]));
% vel_exp_shk(length(vel_exp_shk)-length(fl_exp):end) = []; % Delete last 
% elements so that array has same size as fl_exp.

neg_inds = vel_exp_shk < 0; % Indices in velocity data where velocity is less than 0.
vel_exp_shk(neg_inds) = 0; % Changes negative values to zero. Note I can no longer
% integrate this data to accurately get position. 

vel_exp_mf = medfilt1(vel_exp_shk); % 1D median filter to smooth data.

% figure;
% hold on;
% plot(vel_exp_shk, 'DisplayName', 'shrunk data')
% plot(vel_exp_mf, 'DisplayName', 'median filtered')
% ylabel('Velocity (V)')
% legend()

% Find average of fluorescence across all cells.
fl_exp_avg = mean(fl_exp,2);
fl_exp_std = std(fl_exp,0,2);

figure;
hold on;
plot(vel_exp_mf, 'DisplayName', 'velocity (V)');
eb_ax = shadedErrorBar([],fl_exp_avg, fl_exp_std);
xlabel('frames')
title('Velocity and Avg Fluorescence for Entire Experiment')
legend();

%% Fl of individual cells vs velocity whole experiment

figure;
hold on;
plot(vel_exp_mf, 'DisplayName', 'velocity (V)');
for i = 40:40
    plot(S_all.datasetSm(:,i), 'DisplayName', strcat('Cell', num2str(i)));
end
title('Fl of Cells and Mouse Velocity')
xlabel('Frames')
legend();

% % Compare dF/F with raw fl trace
% rawfl_nrns = S_all.datasetorig(:,S_all.listofneurons);
% plot(rawfl_nrns(:,i)./5000, "DisplayName", 'Cell Raw Fl')


%% SCE Events and velocity

figure;
hold on;
plot(vel_exp_mf, 'DisplayName', 'velocity (V)');
plot(S_all.isSCE.*S_all.num_act_nrns_perframe./40, 'DisplayName', 'isSCE');
title('SCEs and Velocity');
xlabel('Frame')
legend();
