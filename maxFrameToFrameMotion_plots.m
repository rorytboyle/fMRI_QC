% This script plots a histogram of the top 5 max values for the
% derivatives of rotation and translation, and then saves the values and
% their indices in a .mat file

% Author: Rory Boyle rorytboyle@gmail.com
% Date: 06/04/2020

%CHANGE THE PATH TO THE FOLDER YOU WANT
cd('A:\motion_files');
motionFiles=dir('rp*.txt');

destFolder = 'A:\frame_to_frame_motion';

% loop through each file
%%%%%%%%%%%%%%%%%%%% START LOOP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
transAxes = {'x', 'y', 'z'};
rotAxes = {'pitch', 'roll', 'yaw'};

headers = {'x', 'y', 'z', 'x_ix', 'y_ix', 'z_ix', 'pitch', 'roll',...
    'yaw', 'pitch_ix', 'roll_ix', 'yaw_ix'};

for i = 1:length(motionFiles)
    % get ppt code
    subjCode = motionFiles(i).name(5:11);
    
    % get movement parameters file and load
    file = [motionFiles(i).folder filesep motionFiles(i).name];
    mvtparaam=load(file);
    
    % get max (top 5) values for derivatives of translation  (x y z)
    transDerivs = diff(mvtparaam(:,1:3));
    [top5_transDerivs, top5_Trans_ix] = maxk(transDerivs, 5);
       
    subplot(2, 3, [1:3]);
    hist(maxk(top5_transDerivs(:), 5));
    title(["Top 5 translation derivatives - " + subjCode]);
    
    for j = 4:6
        subplot(2, 3, j);
        hist(top5_transDerivs(:, j-3));
        title([transAxes{j-3} + "-axis"]);
    end
    
    savefig([destFolder filesep subjCode '_translation']);
    close(gcf);
    
    % get max (top 5) values for derivatives of rotation (pitch roll yaw)
    rotDerivs = diff(mvtparaam(:,4:6));
    [top5_rotDerivs, top5_Rot_ix] = maxk(rotDerivs, 5);
    
    subplot(2, 3, [1:3]);
    hist(maxk(top5_rotDerivs(:), 5));
    title(["Top 5 rotation derivatives - " + subjCode]);
    
    for j = 4:6
        subplot(2, 3, j);
        hist(top5_rotDerivs(:, j-3));
        title([rotAxes{j-3}]);
    end    

    savefig([destFolder filesep subjCode '_rotation']);
    close(gcf);
    
    % add data to cell array and save
    top5data = [top5_transDerivs, top5_Trans_ix, top5_rotDerivs, top5_Rot_ix];
    top5cell = [headers; num2cell(top5data)];
    
    save([destFolder filesep subjCode '_maxDerivatives.mat'], 'top5cell');
end
    
