function [qc_array, qc_score] = spm_regChecker(files, compare_to, spm_path)
% Loops through .nii files and allows user store a rating of and short 
% note/comment on the file based on either image orientation, 
% co-registration quality, or normalisation quality. Useful for visual QC
% / inspection of images. Type 'end' to complete ratings. (Note: user
% must type 'end' twice in a row to end function). % Saves a .mat
% file 'qc_info.mat' to current directory containing: qc_array (filenames
% and ratings (Bad/Unsure/Good)), score (numerical rating values =0: bad,
% 1: unsure, 2: good.), and last_file_rated (name of the last file rated).
% CAUTION: THIS WILL OVERWRITE ANY PREVIOUSLY SAVED 'qc_info.mat' FILES IN
% CURRENT DIRECTORY. MAKE SURE TO RENAME PREVIOUSLY SAVED FILES BEFORE 
% CALLING FUNCTION.
%
% INPUT:
% files =           (string) file path for directory containing .nii files   
% compare_to =      (string) SPM template image against which you want to
%                   compare. Must be either: 'single_subj_T1', 'avg152T1',
%                   'avg152T2', 'avg305T1', or 'avg152PD'
% spm_path =        (string) file path for SPM directory. Can access this
%                   by typing "which spm" in MATLAB command window. Note:
%                   should be in format:
%                   'C:\Program Files\MATLAB\R2018B\toolbox\spm12' (i.e.
%                   don't provide the link to actual SPM installation).
%
% OUTPUT:
% qc_array =        (cell) files checked * 3 array containing filenames,
%                   corresponding ratings (Bad/Unsure/Good), and
%                   note/comments on file.
% score =           (double) files checked * 1 array containing numerical
%                   ratings. Values:0 = Bad, 1 = Unsure, 2 = Good
% last_file_rated = (string) name of the last file rated
%
% Example usage:    [qc_array, qc_score] = spm_regChecker(...
%                   'C:\files_to_qc ', 'avg152T1', ...
%                   'C:\Program Files\MATLAB\R2018B\toolbox\spm12')
%
% Author: Rory Boyle
% Contact: rorytboyle@gmail.com
% Date: 17/12/2020
%
%% 1) Check 'files' contains .nii files
if iscell(files) == 1
    ismember(files(1), '.nii');
    files_str = char(files(1));
    if strcmp(files_str(end-3:end), '.nii') == 0 % if files don't end in .nii stop function
        fprintf('You have loaded in non .nii files or your cell array of files does not contain any .nii files.\nPlease call function again and load in only .nii files\n')
        return 
    % check if filenames contain a valid file path 
    elseif isfolder(files_str(1:max(strfind(files_str, filesep))-1)) == 0
        fprintf('Your filenames do not contain a valid filepath. \nPlease call function again and load in files with a valid path within their name. \nAlternatively, load in a folder containing only T1 images\n')
        return
    end
end

%% 2) Read in .nii files
% read in filenames from a directory, if files = a cell array of image names,
% then skip and proceed 
if iscell(files) == 0
    file_struct  = dir(files); % read in files
    file_struct = file_struct(3:end); % remove first two rows in struct ("." and "..")
    clear files
    files = cell(length(file_struct),1);
    for i = 1:length(file_struct)
        files{i} = [char(file_struct(i).folder) filesep char(file_struct(i).name)];
        files_str = char(files(i));
    end
    if strcmp(files_str(end-3:end), '.nii') == 0
        fprintf('\nYour folder contains non .nii files. Please load in a folder containing only .nii files.\n') %%% CURRENTLY NOT WORKING PROPERLY
        return
    end
end

%% 3) Get reference file 
ref_image = [spm_path filesep 'canonical' filesep compare_to '.nii']

% check reference file exists/is correctly specified
if exist(ref_image) ~= 2
    error('Reference image (compare_to) is not correctly specified or does not exist in spm_path you have provided')
    return
end
%% 4) Loop through images and call spm CheckReg
% Initialise arrays
qc_array = cell(length(files), 1);
qc_score = zeros(size(files));

user_ends = 0; % variable to break out of while loop

for i = 1:length(files)
    % Get file
    qc_array(i,1) = files(i);
    % Call CheckReg
    tic; spm_check_registration(files{i}, ref_image); toc
    while user_ends == 0
            x = input('\n To end program & save info, type end \n Rate image. Enter 0 for bad, Enter 1 for unsure, Enter 2 for good: \n', 's');
            y = input('\n Add note on file: \n', 's');
           % if y is empty (i.e. no note provided), add empty string
           if isempty(y)
               y = ' ';
           end
           % if image is bad (0) --> add 'Bad' to cell array and 0 to score array
            if strcmp(x, '0')
               qc_array(i, 2) = {'Bad'};
               qc_array(i, 3) = cellstr(y);
               qc_score(i) = 0;
               last_file_rated = files(i);
               save('qc_info', 'qc_array', 'qc_score', 'last_file_rated');
               break
           % if user is unsure about image (0) --> add 'Unsure' to cell array and 1
           % to score array
            elseif strcmp(x, '1')
               qc_array(i, 2) = {'Unsure'};
               qc_array(i, 3) = cellstr(y);
               qc_score(i) = 1;
               last_file_rated = files(i);
               save('qc_info', 'qc_array', 'qc_score', 'last_file_rated');
               break
           % if image is good (2) --> add 'Good' to cell array and 2 to score array
            elseif strcmp(x, '2')
               qc_array(i, 2) = {'Good'};
               qc_array(i, 3) = cellstr(y);
               qc_score(i) = 2;
               last_file_rated = files(i);
               save('qc_info', 'qc_array', 'qc_score', 'last_file_rated');
               break
           % else request input again
            elseif strcmp(x, 'end') == 0
               fprintf('\n Invalid input - Please enter a number between 0-2\n');
            elseif strcmp(x, 'end')
                %%% return qc_array and score for all files up until i
               qc_array = qc_array(1:i-1, :);
               qc_score = qc_score(1:i-1);
               if i>1
                   last_file_rated = files(i-1);
               end
               user_ends = 1;
               break
           break
           end
       end
   if user_ends == 1
       break
   end
   % save variables to qc_info.mat in current directory
if user_ends == 0 
    save('qc_info', 'qc_array', 'qc_score');
elseif user_ends == 1 && i > 1
    save('qc_info', 'qc_array', 'qc_score', 'last_file_rated');
end
end
