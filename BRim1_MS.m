function [S, control] = BRim1_MS(S, control)
% To analyze suite2p flourescence data in context of behavioral data from BRim0_MS.
% 
% Inputs:
% S - structure holding all behavioral data.
% control -
%
% Outputs:
% S - structure holding all behavioral and fluorescence data.
% control - 
%
% Note: BRim0_MS must be run first.
%
% Note, suite2p can save data as a .mat file, which is loaded by this
% function. The imaging data passed to suite2p, and thus the data suite2p
% saves in the .mat file, is recorded at interals of 30 ms (thus at 33.33
% Hz). 


%% Load data
% The .mat file Fall contains all the suite2p data. Variables loaded are:
% F, Fneu, iscell, ops, redcell, spks, and stat. 
s2p_vars = load([S.data_dir, '\suite2p\plane0\Fall.mat']);


% Adding suite2p data to S.
S.medfiltC=500; %normally 500 
S.polynom=1; 
S.medfiltSm=5;
%S.baselineCorrection=[28425, 47855, 54939, 60500];
S.listofneurons=vertcat(1, find(s2p_vars.iscell(:,1)==1)+1); % List of only neurons
S.F=s2p_vars.F; % Cell flourenscence
S.Fneu=s2p_vars.Fneu; % Neuropil flourescence
S.iscell=s2p_vars.iscell; % Binary is cell or neuropil
S.ops=s2p_vars.ops; % Processing options
S.spks=s2p_vars.spks; % Denconvolved spikes (don't use until reading Marius/Carsen papers)
S.stat=s2p_vars.stat; % ROI information

%dataset=[transpose(S.F(1,:)), transpose(S.F)]; 
dataset = transpose(S.F); % Commenting out above line and replacing with this 
% to not have the first neuron included twice in dataset.
dataset=double(dataset); 
S.datasetorig=dataset; 
dataset=dataset(:,transpose(S.listofneurons));
datasetSm=dataset; % Unedited fluorescence of neurons. Each column is a neuron,
% and each row is a frame.

datasetNeu=[transpose(S.Fneu(1,:)), transpose(S.Fneu)]; 
datasetNeu=double(datasetNeu); 
S.datasetNeuorig=datasetNeu; 

% --------------------
% Isn't this again just selecting the neuron fluorescence traces? It seems
% like it was intended to get the neuropil traces. Or is this the
% "neuropil" background that goes along with each neuron? Also, what
% exactly is then the F fluorescence? Is it raw fluorescence? F - 0.7*Fneu?
% ** It seems like it is raw F, because Fneu is subtracted from F later on.
% 
% This is the background neuropil signal around ROIs. 
% --------------------
datasetNeu=datasetNeu(:,transpose(S.listofneurons)); 


S.Fzero=[]; S.FzeroNeu=[];


%% Michael: Exploring the data
% -------------------------------
% For some reason the first neuron trace is not the first ROI in suite2p.
% Instead the second trace in the file datasetSm is the first ROI. After
% that everything lines up, so it seems there is just an extra neuron in
% datasetSm. I think this first trace is from the large central interneuon
% found in CG004_220307.
% -------------------------------

for i = 1:1
    figure();
    hold on
    plot(dataset(:,i))
    title(strcat("Unaltered Data for Cell ", num2str(i)))
    ylabel("Fluorescence (V)")
    xlabel("Frames")
%     plot(datasetNeu(:,i))
%     plot(datasetSm(:,i)-datasetNeu(:,i))
%     legend(["Raw F","Neuropil F","F-Fneu"])
end
% This datasetSm data has artifacts in it of the laser shutter opening and
% closing every 400 frames.

%% correct dataset variable for bleaching and remove artifact

% When I load old data that has the PMT-off, I just need to look and find a
% good threshold number. S.bar is just a threshold to detect and remove
% the PMT-off period artifacts. I could also find a better way to code this.
S.bar=zeros(size(s2p_vars.F,2),1); S.bar(:,:)=100; hold on; plot(S.bar);

% I need to figure out a better way to do the threshold selection, because
% I'm not sure if the PMT-off signal will be the same between different
% experiments. Will need to check. 
thres=dataset(:,1)<S.bar(1); % Using the first cell to detect when the fl trace
% crosses the threshold set by S.bar. The intent is to only capture frames
% when the PMT-off signal occurs.

for k=1:size(dataset,2) % Loops over each neuron. This loop transforms 
    % dataset so that the bleaching is corrected for, the PMT shutter off
    % periods are fixed, and the neuropil signal is subtracted from the raw
    % fluorescence.
    
    if isfield(S,'baselineCorrection')
        for i=1:size(S.baselineCorrection,2)
            dataset(1:S.baselineCorrection(1,i),k)=dataset(1:S.baselineCorrection(1,i),k)+((median(dataset(S.baselineCorrection(1,i)+1:S.baselineCorrection(1,i)+1000,k),1,'omitnan'))-median(dataset(S.baselineCorrection(1,i)-999:S.baselineCorrection(1,i),k),1,'omitnan'));
            datasetNeu(1:S.baselineCorrection(1,i),k)=datasetNeu(1:S.baselineCorrection(1,i),k)+((median(datasetNeu(S.baselineCorrection(1,i)+1:S.baselineCorrection(1,i)+1000,k),1,'omitnan'))-median(datasetNeu(S.baselineCorrection(1,i)-999:S.baselineCorrection(1,i),k),1,'omitnan'));
        end
    end
    
    x=1:size(dataset,1); x=transpose(x); 

    
    y=dataset(:,k); y(thres)=NaN; 
    % dataset(:,k)=dataset including initial frames and PMT-off periods/used for
    % smoothing
    % y=dataset with initial frames and PMT-off periods set to NaN
    
    test2=y(~thres); 
    % test2=dataset with initial frames and PMT-off periods removed; number
    % of points different; used only for Fzero calculation
    
    % ------------------------
    % What is Fzero and why do we fit a gaussian to a histogram of the
    % PMT-off removed data?
    % dF = F - Fzero. Fzero is just the baseline fluorescence. THe gaussian
    % fitting is then done on the histogram of fluorescence to find the
    % mean location of the peak representing Fzero. 
    % ------------------------
    edges=-1000:1:30000;
    histoX=-999.95:1:29999.95;
    histoY=histcounts(test2,edges); hold on; 
    f = fit(transpose(histoX),transpose(histoY),'gauss1');
    S.Fzero(k,1)=f.b1; 
    
    dataset(:,k)=(dataset(:,k)-S.Fzero(k,1))/(S.Fzero(k,1)); %dataset(:,k)=with artifact
    
    y=(y-S.Fzero(k,1))/(S.Fzero(k,1)); %y=without artifact
    
    %filtdata=medfilt1(y,S.medfiltC, 'omitnan'); 
    %c = polyfit(x,filtdata,S.polynom);
    dataset(:,k)=dataset(:,k)-polyval(polyfit(x,medfilt1(y,S.medfiltC, 'omitnan'),S.polynom),x); 
    % --------------
    % Point of using polyfit?
    % Bleaching occurs over time for cells so this corrects that. 
    % OR is this correcting the PMT-off periods? It doesn't seem to remove
    % the PMT-off periods...
    %----------------

    % neuropil correction
    
    y=datasetNeu(:,k); y(thres)=NaN; 
    % dataset(:,k)=dataset including initial frames and PMT-off periods/used for
    % smoothing
    % y=dataset with initial frames and PMT-off periods set to NaN
    
    test2=y(~thres); 
    % test2=dataset with initial frames and PMT-off periods removed; number
    % of points different; used only for Fzero calculation
    
    histoY=histcounts(test2,edges); hold on; 
    %histoY=ha.Values;
    f = fit(transpose(histoX),transpose(histoY),'gauss1');
    S.FzeroNeu(k,1)=f.b1; 
    
    datasetNeu(:,k)=(datasetNeu(:,k)-S.FzeroNeu(k,1))/(S.FzeroNeu(k,1)); %dataset(:,k)=with artifact
    
    y=(y-S.FzeroNeu(k,1))/(S.FzeroNeu(k,1)); %y=without artifact
    
    %filtdata=medfilt1(y,S.medfiltC, 'omitnan'); 
    %c = polyfit(x,filtdata,S.polynom);
    %y_est = polyval(c,x);
    datasetNeu(:,k)=datasetNeu(:,k)-polyval(polyfit(x,medfilt1(y,S.medfiltC, 'omitnan'),S.polynom),x);
    
    %actual neuropil correction
    dataset(:,k)=dataset(:,k)-datasetNeu(:,k);   
    
    edges=-1:0.02:6; 
    histoX=-0.99:0.02:5.99;
    histoY=histcounts(dataset(:,k),edges); 
    f = fit(transpose(histoX),transpose(histoY),'gauss1');
    dataset(:,k)=dataset(:,k)-f.b1;
            
end

%% To visualize the removed artifacts, corrected bleaching, and subtraction of neuropil signal.
for i = 1:1
    figure();
    hold on
    plot(dataset(:,i))
    legend(["F-Fneu"]);
    title("PMT shutter closing artifacts removed and bleaching corrected");
end

%% Signal filtering
startframe=1:400:size(dataset,1);
for k=1:size(dataset,2) % Loops through each neuron.
    
    datasetSm(:,k)=medfilt1(dataset(:,k),S.medfiltSm, 'omitnan'); % medfilt1
    % is a third-order median filter which smooths the voltage signal.
    
    % I'm guessing that in case the filter somehow lowered values below 
    % threshold, this resets those values to NaN.
    dataset(thres,k)=NaN;
    datasetSm(thres,k)=NaN; 
    
%     for i=1:length(startframe)
%         dataset(startframe(i):startframe(i)+2,k)=NaN; 
%         datasetSm(startframe(i):startframe(i)+2,k)=NaN; 
%     end

end

% --------------------
% This smoothing definitely helps, but doesn't completely smooth. Why not
% continue smoothing? Is this just goot enough to fix large random
% fluctuations?
%  - This was done when finding the maximum fluorecence in a bin (S.Rmax).
%  Makes the resulting Rmax value less noise-dependent.
% ----------------------

%% To visualize median-filtered signal.
for i = 1:1
    figure();
    hold on
    plot(datasetSm(:,i))
    legend(["F-Fneu"]);
    title("Smoothed Signal");
end

%%
% Updating variables in S.
S.dataset=dataset;
S.datasetSm=datasetSm; 
S.datasetNeu=datasetNeu; 

% split up df/f data into individual laps 

S.frametimesperlap=cell(1,(length(S.startlevels)-2));
S.dataperlap=cell(1,size(dataset,2));

for k=1:size(dataset,2)  % Loop each cell
    for j=S.lapsused % S.lapsused seems to just be an array of length (# laps)
        % with values equal to indices. So this is looping over each lap.
        
        % Seperates frame timings into laps.
        C=S.frametiming>S.startlevels(j)&S.frametiming<(S.startlevels(j+1))-1; 
        % Saves frame times that fall within each lap.
        S.frametimesperlap{1,j}=S.frametiming(C)-S.startlevels(j);

        % Does th same thing as above; seoerates frame timings into laps.
        D=find(S.frametiming>S.startlevels(j)&S.frametiming<(S.startlevels(j+1))-1); 
        % Saves the unsmoothed fl data into a form seperated by laps.
        S.dataperlap{1,k}{1,j}=dataset(D,k); 
        % Saves the smoothed fl data into lap form.
        S.dataperlapSm{1,k}{1,j} = datasetSm(D,k);
    end
end

% for k=1:size(dataset,2)
%     plot(S.frametiming(:,1),dataset(:,k)); 
%     title(k)
%     pause
%     close all
% end
% -------------------------
% Why is the x-axis on the order of 10^7 here? Why not the normal 10^4 like
% with frames or something smaller for seconds?
% ---------------------------

% clearvars -except S control dataset*

%%
%
% Right now transients from Ca peaks are being included in the heat map.
% This happens sometimes as the mouse starts running, a Ca spike that
% happened just before the mouse started running has its peak excluded (as
% we would like), but the transient activity is still included because it
% lingers into the running time. She would like to exclude these
% transients. 
%
% close all;

% C is a logical array (0s and 1s) with 0 where the velocity is less than
% or equal to 4 cm/s and 1 otherwise. This is later used to leave out all
% fluorescence frames during non-running periods.
C=S.wsALL(:,2)>0.05; %>4cm/sec
S.R=zeros(50, (length(S.startlevelsR)-2), size(dataset,2)); S.R(:,:,:)=NaN; 
S.Rmax=S.R; 

% ------------------------------
% What is S.R and S.Rmax? It seems like S.R and S.Rmax are versions of dF/F...
% - S.R is the mean fluorescence per bin, S.Rmax is the max fluorescnece
% per bin. (50 bins).
% ------------------------------


for k=1:size(dataset,2) % Loops through each cell.
    dataset1=zeros(size(S.wsALL,1), 1); dataset1(:,:)=NaN; 
    dataset1S=dataset1;

    % I believe this is showing how to lay the suite2p frame data (72000
    % data points) over the wavesurfer frame data (3.8x10^7 data points). 
    dataset1(S.frametiming(:,1),1)=dataset(:,k); 
    dataset1S(S.frametiming(:,1),1)=datasetSm(:,k); 
    
    %only running data
    datasetR=dataset1(C);
    datasetRS=dataset1S(C); 
   
    for j=S.lapsused%lapnumber%exclude lap between lap1 and 2 which is too short
        y1=S.startlevelsR(j); y2=(S.startlevelsR(j+1))-1; 
        cg=[]; cg=datasetR(y1:y2,1);  
        cg2=[]; cg2=datasetRS(y1:y2,1); 

        % This section loops through each bin, assigning S.R the mean
        % fluorescence value in that bin. Note this is only done for
        % times when the velocity is above 4 cm/s.
        for binnumber=1:50   
               if (S.threshold((binnumber+1),j)-S.threshold(binnumber,j))>0               
                S.R(binnumber,j,k)=mean(cg(S.threshold(binnumber,j):(S.threshold((binnumber+1),j))-1,1),'omitnan');
                S.Rmax(binnumber,j,k)=max(cg2(S.threshold(binnumber,j):(S.threshold((binnumber+1),j))-1,1));
               end
        end
    end 
end

% % Plots the dF/F fluorescence across laps for each neuron. Also saves
% each plot as a .pdf. 
% for k=1:size(dataset,2)
%     imagesc(transpose(S.Rmax(:,:,k)),'AlphaData',~isnan(transpose(S.R(:,:,k))))
%     colormap jet
%     caxis([0,4])
%     %caxis([0 0.3])
%     colorbar
%     ax=gca;
%     ax.TickDir='out';
%     test=['neuron' num2str(k)];
%     title(test);
%     fileName=[S.s2p_parentfolder, filesep, 'Cell dFF Plots', filesep, 'neuronM', num2str(k), '.pdf'];
%     set(gcf,'Units','inches');
%     set(gcf, 'Position', [0 0 7 10]);
%     saveas(gcf, fileName, 'pdf');
% end

S.MeanR=mean(S.R, 2, 'omitnan');
S.MeanRmax=mean(S.Rmax, 2, 'omitnan');

% clearvars -except S control dataset datasetSm

%%
% Plots to visualize average fluorescence for each neuron.
% for k=2:size(dataset,2)
%     subplot(4,1,[1,2])
%     imagesc(transpose(S.R(:,:,k)),'AlphaData',~isnan(transpose(S.Rmax(:,:,k))))
%     colormap jet
%     caxis([0,2])
%     ax=gca;
%     ax.TickDir='out';
%     test=['neuron' num2str(k)];
%     title(test);
%     subplot(4,1,3)
%     plot(S.MeanR(:,:,k), 'k', 'Linewidth',2)
%     ylabel("MeanR (mean fluorescence?)")
%     subplot(4,1,4)
%     plot(S.MeanRunning, 'r', 'Linewidth',2)
%     ylabel("Mean Running")
%     xlabel("Bins?")
%     ylim([0 inf])
%     pause
% end
% 
%% for k=transpose(S.placecells)
%     plot(S.frametiming(:,1)/S.acq,dataset(:,k)); 
%     title(k)
%     pause
%     close all
% end
% 
%%
% i=0; 
% figure; 
% for k=[12, 18]
%     plot(S.frametiming(:,1)/S.acq,dataset(:,k)-i, 'k', 'Linewidth', 1.1); hold on; 
%     plot(0.00005:0.0001:length(S.wsALL)/S.acq, (S.wsALL(:,6)/10000)+3);
%     i=i+3;
%     xlim([1760 2060])
% end

end