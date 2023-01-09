%function [] = print_standing_lengths(file_dir)
% Prints the amount of time in seconds spent not running (mouse velocity =
% 0) and prints the amount of time spent not running only when preceeded and
% followed by a minimum 5 seconds of standing still. This is to help
% identify good candidate data sets for looking for sharp-wave ripples
% (SWRs) which tend to occur only during extensive periods of rest, and not
% immediately before or after running. **** Shouldn't I be looking for
% long stationary periods coincident with imaging?
%
% Inputs:
% file_dir - the directory contianing the wavesurfer data files to be looped through.
% Each file should contain the wavesurfer data for one experiment.
%
% Note: This function is meant to be run separately from master_script.

file_dir = 'D:\Michael Satchell Rotation\all wavesurfer files';

dir_struct = dir([file_dir, filesep, '*.h5']); % Directory structure of subfolders.
% Makes sure to only include files that end in .h5.



filenames = {dir_struct(:).name};

for i = 1:size(dir_struct, 1) % Loops through folders
    
    ws_data=ws.loadDataFile([file_dir, filesep, filenames{i}]);
    
    acq_rate = ws_data.header.AcquisitionSampleRate;
    vel = ws_data.(names{2}).analogScans(:,1);
    frame_clock = ws_data.(names{2}).analogScans(:,5);

    vel = vel - 1.254; % This is done in BRim0 and i'm not sure why, but likely 
    % to shift all the data so that multiplying by 80 gives the correct
    % velocity in cm/s.

    vel = vel.*80; % shifts velocity units to cm/s.




end


%end