function [S, control] = BRim0_MS(S, data_dir)
% Load and edit wavesurfer data (behavioral data).
%
% Inputs:
% data_dir - file directory to folder that contains wavesurfer and suite2p
% data for a given recording session.
% S - structure to hold behavioral data.
%
% Outputs:
% S - structure containing all relevent processed behavioral data.
% control - 
%
% Note: This function must be run before BRim1_MS.

S.data_dir = data_dir;
S.ws_parentfolder = [data_dir, filesep, 'wavesurfer'];

if ~isfolder(data_dir)
        errorMessage = sprintf('Error: The following folder does not exist:\n%s',data_dir);
        uiwait(warndlg(errorMessage));
        return;  
end

if ~isfolder(S.ws_parentfolder)
        errorMessage = sprintf('Error: The following folder does not exist:\n%s', S.ws_parentfolder);
        uiwait(warndlg(errorMessage));
        return;  
end


d=dir([S.ws_parentfolder, filesep, '*.h5']); 
S.fileNames={d(:).name};

wsCell=cell(1,length(d)); 
S.wsALL=[];
S.startmarkerChange=0;

% Load ws data:
for i=1:length(d)
    startlevels=[];
    data=ws.loadDataFile([S.ws_parentfolder, filesep, S.fileNames{i}]);
    S.acq=data.header.AcquisitionSampleRate;
    names=fieldnames(data);
    ws1=[];
    ws1(:,2)=data.(names{2}).analogScans(:,1);
    ws1(:,4)=data.(names{2}).analogScans(:,3);
    ws1(:,5)=data.(names{2}).analogScans(:,4);
    ws1(:,1)=data.(names{2}).analogScans(:,5);
    ws1(:,3)=data.(names{2}).analogScans(:,6);
         
%      if i==1 
%          S.startmarkerChange=1; %1==true
%          ws1(1,1)=0; ws1(2:10,1)=5;     
%      end
          
    if i==1 %|| i==2 || i==3
        ws1(:,5)=0;
    end

    idxl = ws1(:,1)>=2;
    idxl(1) = 0;
    idx = find(idxl);
    yest = ws1(idx-1,1)<2; 
    startlevels=idx(yest)-1; 
    
    if i>1
        wsCell{1,i}=ws1(:,1:5);
    else        
        wsCell{1,i}=ws1(startlevels(1):length(ws1),1:5);
    end
    
    S.wsALL=vertcat(S.wsALL, (wsCell{1,i})); 
end


% wavesurfer data is loaded into S.wsALL. Columns contain information as
% below: 
% --------------------------------------------
% Ask for more information on what each column holds, including
% the units used. 
% --------------------------------------------
% All units are volts
% S.acq is the data acquisition rate
% S.wsALL(:,1) = lap start
% S.wsALL(:,2) = running
% S.wsALL(:,3) = valve opening
% S.wsALL(:,4) = licking 
% S.wsALL(:,5) = frame clock 
% S.wsALL(:,6) = position, integrated signal from running


% Edit ws data:

S.wsALL(:,2)=S.wsALL(:,2)-1.254;
S.wsALL(4.313e6:4.315e6,5)=0;

% add startmarkers here manually 
%S.wsALL(x:x+10,1)=5; 

% number of recorded frames

idxl = S.wsALL(:,5)>=2; % When the microscope images, the frameclock signal goes up to about 5,
% so selecting frameclock values above 2 sorts out times at which the
% imaging took place. 
idxl(1) = 0;
idx = find(idxl);
yest = S.wsALL(idx-1,5)<2; % I believe this finds the start times of each imaging frame
% by looking for times where the i-1 frameclock value is less than 2 and
% the i value is above 2. 
S.frametiming=idx(yest); % start times (in wavesurfer recording steps) of each imaging frame.

length(S.frametiming) % Shows the number of imaging frames.

% This does somethign similar to above, but for laps. Gets the start times
% of each lap.
idxl = S.wsALL(:,1)>=2;
idxl(1) = 0;
idx = find(idxl);
yest = S.wsALL(idx-1,1)<2; 
S.startlevels=idx(yest)-1;
S.startlevels(length(S.startlevels)+1,1)=length(S.wsALL(:,1))+1; 

%remove licks shorter than 20 ms

S.wsALL(1:5,4)=0; S.wsALL(length(S.wsALL),4)=0; 
idxl = S.wsALL(:,4)>=2; 
idxl(1) = 0;
idx = find(idxl);
yest = S.wsALL(idx-1,4)<2;
lickdetec(:,1)=idx(yest); 

idxl = S.wsALL(:,4)<2; 
idxl(1) = 0;
idx = find(idxl);
yest = S.wsALL(idx-1,4)>2;
lickdetec(:,2)=idx(yest)-1; 

lickdetec(:,3)=(lickdetec(:,2)-lickdetec(:,1)) *1/S.acq; 

for i=1:length(lickdetec)
    if lickdetec(i,3)<0.02
        S.wsALL(lickdetec(i,1):lickdetec(i,2),4)=0; 
    end
end

S.wsALL(1:5,3)=0; S.wsALL(length(S.wsALL),3)=0; 

idxl = S.wsALL(:,3)>=2; 
idxl(1) = 0;
idx = find(idxl);
yest = S.wsALL(idx-1,3)<2;
valvedetec(:,1)=idx(yest); 

idxl = S.wsALL(:,3)<2; 
idxl(1) = 0;
idx = find(idxl);
yest = S.wsALL(idx-1,3)>2;
valvedetec(:,2)=idx(yest)-1; 

valvedetec(:,3)=(valvedetec(:,2)-valvedetec(:,1)) *1/S.acq; 

for i=1:length(valvedetec)
    if valvedetec(i,3)<0.010
        S.wsALL(valvedetec(i,1):valvedetec(i,2),3)=0; 
    end
end

figure; plot(valvedetec(:,3)); ylim([0 0.2]); 
title('duration of valve openings'); 
ylabel('valve opening time (s)'); xlabel ('lap number'); 


% divide up into individual laps and calculate position signal from analog running signal 

S.wsL=cell(1,size(S.wsALL,2)); % wsL is the same information as wsALL, just
% divided up into laps.


% ------------------------------
% What is S.startlevels? - points in time at which laps start, by
% thresholding the lap start signal. Units are in acquisition points
% ------------------------------

for k=1:size(S.wsALL,2)
    for i=1:(length(S.startlevels)-1)
        S.wsL{1,k}{:,i}=S.wsALL(S.startlevels(i):S.startlevels(i+1)-1,k); 
    end
end

S.maxpos=zeros(length(S.startlevels),1); 

position=S.wsL{1,2};
for i=1:(length(S.startlevels)-1)
    position{1,i}=cumtrapz(position{1,i}); % Integration over velocity to get position.
    S.maxpos(i,1)=max(position{1,i}); 
end

S.wsL{1,size(S.wsALL,2)+1}=position; 

interim=[]; 
for k=1:length(S.startlevels)-1
    interim=vertcat(interim, (position{1,k})); 
end

S.wsALL=horzcat(S.wsALL, interim); %S.wsALL(:,8) is position

control=[S.maxpos S.startlevels];

figure; plot(control(:,1));
title('lap length');
xlabel('lap number'); ylabel('lap length (integrated V)');


S.framelength=zeros(length(S.frametiming),3); S.framelength(:,:)=NaN; 
idxl = S.wsALL(:,5)>=2; 
idxl(1) = 0;
idx = find(idxl);
yest = S.wsALL(idx-1,5)<2;
S.framelength(:,1)=idx(yest); % Sets the first column to be the imaging 
% frame start times (in ws units). The same as how S.frametimes
% is defined. 

% This does the inverse of above. First finds all ws times i where there was
% no imaging, then finds those times i-1 where there was imaging. This
% sorts out the end times of imaging frames. 
idxl = S.wsALL(:,5)<2; 
idxl(1) = 0;
idx = find(idxl);
yest = S.wsALL(idx-1,5)>2;
S.framelength(:,2)=idx(yest)-1; % stores end times here.
S.framelength(:,3)=(S.framelength(:,2)-S.framelength(:,1)) *1/S.acq; % Finds
% the differences between each start and end time. Should always be the
% same as each frame should take the same amount of time to capture.

figure; plot(S.framelength(:,3),'o'); ylim([0 0.05]);
ylabel('frameduration (s)'); xlabel('frame number');
title('frame duration')

%adjust frametimes to times with more running
S.frametimingOriginal=S.frametiming;
for k=1:length(S.frametiming)
    if S.wsALL(S.frametiming(k,1),2)<=0.01
        c=[]; 
        c=find(S.wsALL(S.frametiming(k,1):S.framelength(k,2),2)>0.01, 1,'first');
            if isempty(c)==0 && c<(0.032*S.acq) %length of 1 frame
                S.frametiming(k,1)=S.frametiming(k,1)+c-1;
            end
    end
end

% clearvars -except S control


%%
% --------------------------
% What is the purpose of this section?
% ----------------------------

% close all
C=S.wsALL(:,2)>0.05; % 0.05 V = 4 cm/sec

allstartR=S.wsALL(C,1); 
allstartR(1,1)=0;allstartR(2:10,1)=5;
allstartR(3617167+73024:3617167+73024+10,1)=5;
idxl = allstartR>=2;
idxl(1) = 0; idx = find(idxl);
yest = allstartR(idx-1,1)<2; 
S.startlevelsR=idx(yest)-1;
S.startlevelsR(length(S.startlevelsR)+1,1)=length(allstartR)+1;

S.wsR=cell(1,size(S.wsALL,2));

for k=1:size(S.wsALL,2)
    cg=S.wsALL(C,k);
    for i=1:(length(S.startlevelsR)-1)
        S.wsR{1,k}{:,i}=cg(S.startlevelsR(i):S.startlevelsR(i+1)-1,1); 
    end
end


S.maxposR=zeros(length(S.startlevelsR),1); 
for i=1:(length(S.startlevelsR)-1)
    cg=S.wsR{1,6}{1,i};
    cg=cg(2:length(cg),1);
    S.maxposR(i,1)=max(cg); 
end

control=[control S.maxposR S.startlevelsR];

S.lapsused=[1:(length(S.startlevelsR)-2)];

%is number of startlevelsR and startlevels the same? 

% clearvars -except S control

%%
%for each lap calculate bin dimensions

S.threshold=zeros(50,length(S.startlevelsR)-2); % must have 51 rows at the end
 
for i=1:(length(S.startlevelsR)-2)
    
    cg=S.wsR{1,6}{1,i};
    
    for binnumber=1:50
        thres=binnumber * (S.maxposR(i,1)/50); 
        idxl = cg>=thres;
        idxl(1) = 0;
        idx = find(idxl, 1);
        yest = cg(idx-1,1)<thres; 
        if idx(yest)>0
            S.threshold(binnumber,i)=idx(yest);
        else
            S.threshold(binnumber,i)=numel(cg);
        end
    end
end

interim=zeros(1,(length(S.startlevelsR))-2); interim(:,:)=1; 
S.threshold=vertcat(interim, S.threshold);

%behavioral parameters

% licking 

S.L=zeros(50,(length(S.startlevelsR)-2)); S.L(:,:)=NaN; 
S.Lrate=zeros(50,(length(S.startlevelsR)-2)); S.Lrate(:,:)=NaN; 
S.Ltick=cell(1,length(S.startlevelsR)-2); 

for j=S.lapsused
    cg=S.wsR{1,4}{1,j}; 
    
    for binnumber=1:50
        cgg=cg(S.threshold(binnumber,j):S.threshold(binnumber+1,j)-1,1);
        idxl = cgg>=2; idxl(1) = 0;
        idx = find(idxl);
        yest = cgg(idx-1,1)<2; 
        
        if numel(idx(yest))>0
            S.L(binnumber,j)=numel(idx(yest)); 
            S.Lrate(binnumber,j)=numel(idx(yest))/(length(cgg)/S.acq); 
        else
            S.L(binnumber,j)=0;    
            S.Lrate(binnumber,j)=0;  
        end
    end
    
    %licks as a function of position; for each lick determine position;     
    idxl = cg>=2; 
    idxl(1) = 0;
    idx = find(idxl);
    yest = cg(idx-1,1)<2; 
    licktiming=idx(yest);
    
    position=S.wsR{1,6}{1,j};
    S.Ltick{1,j}=position(licktiming); 
end

S.MeanL=mean(S.L, 2, 'omitnan');
S.MeanLrate=mean(S.Lrate, 2, 'omitnan');

%valve opening

S.V=zeros(50,(length(S.startlevelsR)-2)); S.V(:,:)=nan; 
S.Vtick=cell(1,length(S.startlevelsR)-2); 

for j=S.lapsused
    y1=S.startlevelsR(j); y2=(S.startlevelsR(j+1))-1; 
    cg=S.wsR{1,3}{1,j}; 
    for binnumber=1:50
        cgg=cg(S.threshold(binnumber,j):S.threshold(binnumber+1,j)-1,1);
        idxl = cgg>=2; 
        idxl(1) = 0;
        idx = find(idxl);
        yest = cgg(idx-1,1)<2; 
        
        if numel(idx(yest))>0
            S.V(binnumber,j)=numel(idx(yest)); 
        else
            S.V(binnumber,j)=0;    
        end
    end
    
    %licks as a function of position; for each lick determine position;
    idxl = cg>=2; 
    idxl(1) = 0;
    idx = find(idxl);
    yest = cg(idx-1,1)<2; 
    valvetiming=idx(yest);
    
    position=S.wsR{1,6}{1,j};
    S.Vtick{1,j}=position(valvetiming); 
end

S.SumValve=sum(S.V, 2, 'omitnan');
S.MeanValve=mean(S.V, 2, 'omitnan');

% velocity profile 

S.Running=zeros(50,(length(S.startlevels)-2)); %50 x number of laps
S.Running(:,:)=nan; 

for j=S.lapsused % Loop through laps.
    cg=S.wsR{1,2}{:,j}; % wavesurfer velocity data for that lap.
    for binnumber=1:50 % Loops bins.
        % The ws data has way more data points in it than the imaging data,
        % so I think this portion is for binning the ws data. 
        cgg=cg(S.threshold(binnumber,j):S.threshold(binnumber+1,j)-1,1);
        % To reduce the size of the ws data, the mean within each bin is
        % taken of the velocity, and that's what becomes S.Running.
        S.Running(binnumber, j)=mean(cgg,'omitnan');
    end
end

S.MeanRunning=mean(S.Running,2,'omitnan');

figure; 
subplot(3,2,1)
imagesc(transpose(S.Running)*80); 
ax=gca;
set(ax, 'XTick', [0.5 13 25.5 38 50.5]); 
set(ax, 'XTickLabel',({'0', '45', '90','135','180'}));
ylabel('# lap'); 
colorbar
caxis([0 inf])
title ('velocity (cm/sec)')

subplot(3,2,2)
plot(S.MeanRunning*80, 'b', 'Linewidth',2); 
ax=gca;
set(ax, 'XTick', [1 6.25 12.5 18.75 25 31.25 37.5 43.75 50]); 
set(ax, 'XTickLabel',({'0', '', '45', '', '90','','135','','180'}));
xlabel('position (cm)'); ylabel('velocity (cm/sec)'); 
ylim([0 inf]);

subplot(3,2,3); 
xxtick=S.Vtick; 
for j=S.lapsused
    xxtick{1,j}=-j; 
    if isempty(S.Vtick{1,j})<1
       plot(S.Vtick{1,j}, xxtick{1,j}, 'r+'); 
    end
    ylim([(((length(S.startlevelsR)-2)+1)*-1) 1]); 
    MarkerSize=1';
    hold on; 
end
xlabel('position'); ylabel('# lap'); 
xlim([0 max(S.wsALL(:,6))])

subplot(3,2,4)
plot(S.SumValve, 'r', 'Linewidth',2); 
ax=gca;
set(ax, 'XTick', [1 6.25 12.5 18.75 25 31.25 37.5 43.75 50]); 
set(ax, 'XTickLabel',({'0', '', '45', '', '90','','135','','180'}));
xlabel('position (cm)'); ylabel('number of valve openings'); 
ylim([0 inf]);

subplot(3,2,5); 
xtick=S.Ltick; 
for j=S.lapsused
    xtick{1,j}=-j; 
    if isempty(S.Ltick{1,j})<1%length(S.Ltick{1,j})>0
       plot(S.Ltick{1,j}, xtick{1,j}, 'k+'); 
    end
    ylim([(((length(S.startlevelsR)-2)+1)*-1) 1]); 
    MarkerSize=1';
    hold on; 
end
xlim([0 max(S.wsALL(:,6))])
xlabel('position'); ylabel('# lap'); box off; 

subplot(3,2,6)
plot(S.MeanLrate, 'k', 'Linewidth',2); 
ax=gca;
set(ax, 'XTick', [1 6.25 12.5 18.75 25 31.25 37.5 43.75 50]); 
set(ax, 'XTickLabel',({'0', '', '45', '', '90','','135','','180'}));
xlabel('position (cm)'); ylabel('lick rate (Hz)'); 
ylim([0 inf]);

% clearvars -except S control

end