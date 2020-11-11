% This script creates motion plots for participants who survive a threshold
% of mean framewise displacement > 0.5mm. The scale for translation and 
% rotation are set to the total min and max value in that plane (i.e. the
% min x / y /z movement across all participants)
%
% Author: Rory Boyle rorytboyle@gmail.com
% Date: 09/04/2020

%% 1) Read in files
cd('A:\motion_files_thresholded');
motionFiles=dir('rp*.txt');

destFolder = 'A:\motion_plots';

%% 2) Get min and max values in each axis (x,y,z,pitch,roll,yaw) to set limits
all_mvtparaam = [];
all_mvtparaam_diffs = [];

% Loop through files and collate all values
for i = 1:length(motionFiles)
    % get ppt code
    subjCode = motionFiles(i).name(5:11);
    
    % create subdir in destination folder for ppt
%     mkdir([destFolder filesep subjCode]);
    
    % get movement parameters file and load
    file = [motionFiles(i).folder filesep motionFiles(i).name];
    mvtparaam=load(file);
    
    % get differences
    mvtparaam_diff = diff(mvtparaam);
    
    % add to all arrays
    all_mvtparaam = [all_mvtparaam; mvtparaam];
    all_mvtparaam_diffs = [all_mvtparaam_diffs; mvtparaam_diff];
end

% Get min + max values for translation and rotation, and for differences
max_trans = max(max(all_mvtparaam(:,1:3)));
min_trans = min(min(all_mvtparaam(:,1:3)));
max_rot = max(max(all_mvtparaam(:,4:6)));
min_rot = min(min(all_mvtparaam(:,4:6)));

max_trans_diff = max(max(all_mvtparaam_diffs(:,1:3)));
min_trans_diff = min(min(all_mvtparaam_diffs(:,1:3)));
max_rot_diff = max(max(all_mvtparaam_diffs(:,4:6)));
min_rot_diff = min(min(all_mvtparaam_diffs(:,4:6)));

%% 3) Plot
% loop through each file and create plots

for i = 1:length(motionFiles)
    % get ppt code
    subjCode = motionFiles(i).name(5:11);
    
    % create subdir in destination folder for ppt
%     mkdir([destFolder filesep subjCode]);
    
    % get movement parameters file and load
    file = [motionFiles(i).folder filesep motionFiles(i).name];
    mvtparaam=load(file);

    % plot movement in x y z 
    subplot1 = subplot(2, 2, 1);
    plot(mvtparaam(:,1:3));
    ylim(subplot1, [min_trans max_trans]);
    title(['Translation - ' subjCode]);

    % plot first derivative of movement in x y z
    subplot2 = subplot(2, 2, 2);
    plot(diff(mvtparaam(:,1:3)));
    ylim(subplot2, [min_trans_diff max_trans_diff]);
    title(['Derivative of Translation - ' subjCode]);
    legend({'x', 'y', 'z'}, 'Location', 'eastoutside');
    legend('boxoff');

    % plot movement in pitch roll yaw 
    subplot3 = subplot(2, 2, 3);
    plot(mvtparaam(:,4:6));
    ylim(subplot3, [min_rot max_rot]);
    title(['Rotation - ' subjCode]);
    
    % plot first derivative of movement in pitch roll yaw
    subplot4 = subplot(2, 2, 4);
    plot(diff(mvtparaam(:,4:6)));
    ylim(subplot4, [min_rot_diff max_rot_diff]);
    title(['Derivative of Rotation - ' subjCode]);
    legend({'pitch', 'roll', 'yaw'}, 'Location', 'eastoutside');
    legend('boxoff');
    
    savefig([destFolder filesep subjCode]);
    close(gcf);
 end