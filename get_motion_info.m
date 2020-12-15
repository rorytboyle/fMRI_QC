function [motion_info, Power_exclusions, all_fd_arrays, ...
    all_euler_arrays, all_dist_arrays]=...
    get_motion_info(rootDir, numVols, save_output, save_excel, save_dir)
% This function calculates mean framewise displacement based on SPM
% realigment parameters (listed in rootDir)
%
% INPUT:
% rootDir =             (string) directory containing realignment
%                       parameter files generated during SPM preprocessing.
%                       Code assumes that all realignment files are named
%                       following the the format: 
%                       rp_aSUBID_TASKNAME_00001.txt where 'SUBID' =
%                       participant code/id and 'TASKNAME' = name of fMRI
%                       task/scan e.g. rp_a2029_PaperFold_00001.txt
% numVols =             (double) Number of volumes in each fMRI scan with 
%                       realignment parameters listed in rootDir. 
% save_output =         (boolean) flag (true/false) to save output to .mat
%                       file (if true, will save 
%                       motion_info to rootDir)
% save_excel =          (boolean) flag (true/false) to save motion_info 
%                       output as excel file (if true, will save 
%                       motion_info to rootDir)
% save_dir =            (string) directory to save output to
%
% OUTPUT:
% motion_info =         (cell array) number of participant * 7 with columns
%                       in following order: 
%                       'subcode' = participant ids
%                       'mean_FWdisplacement = mean FWD
%                       'mean_distance' = mean absolute
%                       distance/displacement
%                       'peak_distance' = max absolute
%                       distance/displacement
%                       'mean_rotation' = mean rotation
%                       'num_moves_0.1mm' = total number of movements with
%                       displacements > 0.1 mm
%                       'num_moves_0.5mm' = total number of movements with
%                       displacements > 0.5 mm
%                       'Percent_bad' = % of volumes with absolute
%                       displacements > 0.5 mm
% Power_exclusions =    (cell array) n * 1 containing participant codes / 
%                       subids of participants who have over 25% of volumes
%                       with absolute displacements greater than 0.5mm.
%                       based on Power et al. (2012) criterion
% all_fd_arrays =       (double) number of volumes * number of ppts 
%                       containing framewise displacements. Calculated 
%                       as outlined by Power et al. (2012) NeuroImage
% all_euler_arrays =    (double) number of volumes * number of ppts 
%                       containing mean rotations. Calculated as Euler
%                       angles.
% all_dist_arrays =     (double) number of volumes * number of ppts
%                       containing mean translations. Calculated as square
%                       root of squared difference in X, Y, and Z
%                       axis between each volume
%
% Primary authors: Kathy Ruddy & Daniel Wooley
% Other authors: Rob Whelan & Rory Boyle
%                (edits to generalise code, convert to function,
%                add further outputs and comments/documentation)
% Contact: rorytboyle@gmail.com
% Date: 14/12/2020
%
% Power et al. (2012) Neuroimage: https://doi.org/10.1016/j.neuroimage.2011.10.018
% Example usage
% [motion_info, Power_exclusions, all_fd_arrays, all_euler_arrays,...
%    all_dist_arrays] = headMotion_generic('A:\realignment_parameters',...
%   430, true, true, 'A:\motion_info');

%% 1) Prep data and preallocate arrays
% Get realignment parameter text files
cd(rootDir);
directory_contents=dir('rp*.txt');

% set thresholds for displacements (in mm)
moveThresh1 = 0.1;
moveThresh2 = 0.5;

% preallocate arrays
final_subcodes=cell(length(directory_contents),1);
meanFWdisplacement_Output=NaN(1,length(directory_contents));
meanDistOutput = NaN(1,length(directory_contents));
peakDistOutput = NaN(1,length(directory_contents));
numMoves1Output = NaN(1,length(directory_contents));
eulerMeanOutput = NaN(1,length(directory_contents));
numMoves2Output = NaN(1,length(directory_contents));
numScanRegOutput = NaN(1,length(directory_contents));
percent_badOutput = NaN(1,length(directory_contents));

all_fd_arrays = NaN(numVols-1, length(directory_contents));
all_euler_arrays = NaN(numVols-1, length(directory_contents));
all_dist_arrays = NaN(numVols-1, length(directory_contents));

%% 2) Loop through realignment parameters and calculate motion info
for ii = 1:length(directory_contents)
    clc;disp(ii);
    % get subject =
    subjcode= extractBefore(extractAfter(directory_contents(ii).name,...
        'rp_a'), '_');
    % load motion parameters file
    headmove_filename=[directory_contents(ii).folder filesep directory_contents(ii).name];
    file=load(headmove_filename);
    % Check numVols is same as number of volumes/frames in realignment
    % parameters files
    if numVols ~= length(file)
        error(['numVols is not equal to the number of volumes/frames in'...
            ' the realignment parameters file'])
    end
    % create paramter variables
    x = file(:,1);
    y = file(:,2);
    z = file(:,3);
    p = file(:,4);
    r = file(:,5);
    yw = file(:,6);

    % find the difference between succesive volumes
    xDiff = diff(x);
    yDiff = diff(y);
    zDiff = diff(z);
    pDiff = diff(p);
    rDiff = diff(r);
    ywDiff = diff(yw);

    % square the difference
    xSq = xDiff.^2;
    ySq = yDiff.^2;
    zSq = zDiff.^2;
    pSq = pDiff.^2;
    rSq = rDiff.^2;
    ywSq = ywDiff.^2;

    % calculate translations
    distArray = sqrt(xSq+ySq+zSq);
    all_dist_arrays(:,ii) = distArray;
    
    % find mean and peak displacement
    meanDist = mean(distArray);
    peakDist = max(distArray);

    % number of translations > moveThresh
    if isempty(find(distArray>moveThresh1))
        numMoves1 = 0;
    else
        numMoves1 = size(find(distArray>moveThresh1),1);
    end

    % calculate rotations (euler angle)
    % euler angle = arccos((cos(phi)cos(theta)+cos(phi)cos(psi)+cos(theta)cos(psi)+sin(phi)sin(psi)sin(theta)-1)./2)
    eulerArray = acos((cos(abs(ywDiff)).*cos(abs(pDiff))+cos(abs(ywDiff)).*cos(abs(rDiff))+cos(abs(pDiff)).*cos(abs(rDiff))+sin(abs(ywDiff)).*sin(abs(rDiff)).*sin(abs(pDiff))-1)./2);
    all_euler_arrays(:,ii) = eulerArray;
    
    % find mean rotation
    eulerMean = mean(eulerArray);

    % calculate framewise displacement - powers et al. 2012 NeuroImage
    fdArray = abs(xDiff) + abs(yDiff) + abs(zDiff) + abs(100*pi*((ywDiff.*(180/pi))./360)) + abs(100*pi*((pDiff.*(180/pi))./360)) + abs(100*pi*((rDiff.*(180/pi))./360));
    all_fd_arrays(:,ii) = fdArray;

    % find mean framewise displacement
    mean_framewise_displacement=nanmean(fdArray);
    
    % number of framewise displacements > moveThresh
    if isempty(find(fdArray>moveThresh2))
        numMoves2 = 0;
        %         R = file(:,:);
        numScanReg = 0;
    else
        numMoves2 = size(find(fdArray>moveThresh2),1);

        % create regressor and save new regressor file
        move2reg = zeros(size(x,1),1);

        for ll = 1:size(fdArray,1)
            if fdArray(ll) > moveThresh2
                if ll == 1
                    move2reg(ll:ll+2) = 1;
                elseif (ll > 1) && (ll < (size(fdArray,1)))
                    move2reg(ll-1:ll+2) = 1;
                elseif ll == (size(fdArray,1))
                    move2reg(ll-1:ll+1) = 1;
                end
            end
        end

        %         R = [file(:,:), move2reg];
        numScanReg = size(find(move2reg == 1),1);

    end

    %calculate the percentage of all scans that would be lost using the
    %Power 2012 criterion        
    percent_bad=(numMoves2/numVols)*100;

    %make big vectors with the outputs for all ppts
    meanDistOutput(ii) = meanDist;
    peakDistOutput(ii) = peakDist;
    numMoves1Output(ii) = numMoves1;
    eulerMeanOutput(ii) = eulerMean;
    numMoves2Output(ii) = numMoves2;
    numScanRegOutput(ii) = numScanReg;
    percent_badOutput(ii) = percent_bad;

    meanFWdisplacement_Output(ii) = mean_framewise_displacement;

    final_subcodes{ii}=subjcode;
end

%% 3) Find ppts who exceed Power criterion
% get indices
power_exclusions_ix = find(percent_badOutput>25);
% get participant codes/ids
Power_exclusions = final_subcodes(power_exclusions_ix);

%% 4) Return output
total_header={};
total_header=[total_header, 'subcode', 'mean_FWdisplacement',...
    'mean_distance','peak_distance','mean_rotation',...
    'num_moves_0.1mm','num_moves_0.5mm','Percent_bad'];
merged_data=[final_subcodes,num2cell(meanFWdisplacement_Output'),...
    num2cell(meanDistOutput'), num2cell(peakDistOutput'), ...
    num2cell(eulerMeanOutput'), num2cell(numMoves1Output'),...
    num2cell(numMoves2Output'),num2cell(percent_badOutput')];

% Add header/column names to data
motion_info=[total_header;merged_data];

% save output if specified
if save_output == true
    save([save_dir filesep 'all_motion_info.mat'], 'motion_info',...
        'Power_exclusions', 'all_fd_arrays', 'all_euler_arrays',...
        'all_dist_arrays');
end

% Write to excel file if specified
if save_excel == true
      xlswrite([save_dir filesep 'head_motion.xlsx'],motion_info);
end

