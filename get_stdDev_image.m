%-----------------------------------------------------------------------
% This script loops through preprocessed images, gets the residualised 
% volumes and calculates a standard deviation image
% Author: Rory Boyle rorytboyle@gmail.com
% Date: 03/04/2020
%-----------------------------------------------------------------------
clear all
%% Read in folders and files
folders = {'A:\batch1_preprocessed', 'A:\batch2_preprocessed',...
           'A:\batch3_preprocessed', 'A:\batch4_preprocessed',...
           'A:\batch5_preprocessed'}
       
% create cell array of all dir names (with full paths) i.e. 1 dir per ppt
ppt_dir = {};
for folder = 1:length(folders) % loop through folders
    fileStruct = dir(folders{folder});
    for file = 1:length(fileStruct) % loop through files and append
        fileName = [fileStruct(file).folder filesep fileStruct(file).name];
        if fileName(end) ~= '.' % don't add . and .. to list of filenames
            ppt_dir = [ppt_dir, fileName];
        end
    end
end

       
% get list of Res_*.nii files for each ppt, get subid, calculate mean image
for i = 1:length(ppt_dir)
    subDir = ppt_dir{i};
    subid = extractAfter(subDir, "preprocessed\");
    
    % list out residualised images i.e. *.Res_*.nii files
    cd(subDir)
    imgs = dir('Res_*.nii');
    img_list = {};
    for img = 1:length(imgs)
        imgName = [imgs(img).folder filesep imgs(img).name ',1'];
        img_list = [img_list; imgName];
    end
    
    %% Run imcalc and calculate mean image   
    spm('defaults', 'fmri');
    spm_jobman('initcfg');
    
    matlabbatch{1}.spm.util.imcalc.input = img_list;
    matlabbatch{1}.spm.util.imcalc.output = [subDir filesep subid '_stdDev'];
    matlabbatch{1}.spm.util.imcalc.outdir = {''};
    matlabbatch{1}.spm.util.imcalc.expression = 'std(X)';
    matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
    matlabbatch{1}.spm.util.imcalc.options.dmtx = 1;
    matlabbatch{1}.spm.util.imcalc.options.mask = 0;
    matlabbatch{1}.spm.util.imcalc.options.interp = 1;
    matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
    
    spm_jobman('run',matlabbatch);
    clear matlabbatch
end       

