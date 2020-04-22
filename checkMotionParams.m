function [meanmove, maxmove, moveprc, ppts_to_remove]=...
    checkMotionParams(data, prct, subids)
% This function dentifies participants to remove based on percentile of
% frame-to-frame motion in translation and rotation axes.
%
% INPUT:
% data =            Timepoint * Parameters * Participant matrix containing 
%                   frame to frame motion parameters (i.e. realignment 
%                   parameters in x, y, z, pitch, roll, and yaw directions
%                   and their framewise differences) in each timepoint/vol
%                   for each participant. Should be a timepoint * 12 * 
%                   participant matrix.
% prct =            Percentile at which to establish cut off
% subids =          1 * Participant cell array containing subids in same
%                   order as in data
%
% OUTPUT:
% meanmove =        Participant * 12 array containing mean value for motion
%                   parameters at each timepoint/vol.
% maxmove =         Participant * 12 array containing max value for motion
%                   parameters at each timepoint/vol.
%                   frame motion in 
% moveprc =         1 * 12 array containing value of each motion parameter 
%                   at specified percentile.
% ppts_to_remove =  1 * n cell array containing subids of participants with
%                   motion parameter values beyond the percentile
%                   threshold.

% Author: Rob Whelan
% Edits and comments: Rory Boyle rorytboyle@gmail.com
% Date: 19/04/2020

%% 1) Calculate mean and max motion values and the values at specific percentile
close all;

meanmove=squeeze(mean(abs(data),1))';%get the average of absolute movement for all subs, transposed b/c easier to see
maxmove=squeeze(max(abs(data),[],1))';%get the average of absolute movement for all subs, transposed b/c easier to see
moveprc=squeeze(prctile(maxmove,[prct],1)); %get values at specific percentile

%% 2) Make bar charts to show percentiles and max movement values
% Quick Bar Chart of Percentile for X Y Z Derivatives 
figure;
bar(moveprc([7:9]));    
title('prcs x y z deriv');

% Quick Bar Chart of Percentile for Pitch Roll Yaw Derivatives
figure;
bar(moveprc([10:12]));    
title('prcs p r y deriv');

% Quick Bar Chart of Max Movement for X Y Z Derivatives
figure;
for n=7:9
    subplot(1,3,n-6);
    hist(maxmove(:,n));    
end
title('max move deriv x y z');

% Quick Bar Chart of Max Movement for Pitch Roll Yaw Derivatives
figure;
for n=10:12
    subplot(1,3,n-9);
    hist(maxmove(:,n));    
end
title('max move deriv p r y');

%% 3) Get indices of participants to remove
data_prct=data;%to start
ppts_ix=[];
for n=1:12  % loop through each motion parameter
    for m=1:size(data, 1)  % loop through each timepoint
        for s=1:size(data, 3)  % loop through each participant
            if data_prct(m,n,s)>moveprc(n)
                data_prct(:,n,s)=NaN;
                ppts_ix=[ppts_ix s];
            end
        end
    end
end
ppts_ix=unique(ppts_ix);

%% 4) Return list of participant ids/codes
ppts_to_remove = subids(ppts_ix);