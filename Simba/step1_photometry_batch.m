%% Blackbox Batch Processing Code - Centralized Output %%
% Processing TDT fiber photometry data for keypoint moseq blackbox
% Saves all results to a single dated folder in the Parent Directory

%% Setup and Directory Selection %%
clear all;
clc;

% --- CONFIGURATION ---
% Update these paths if you move your SDK or AUC folders
SDKpath = 'C:\Users\Arlene\Documents\Margolis Lab\Fiber photometry\FP codes\TDTMatlabSDK-master';
aurocpath = 'C:\Users\Arlene\Documents\Margolis Lab\Fiber photometry\FP codes\AUC';
addpath(genpath(SDKpath));
addpath(genpath(aurocpath));

% Select the PARENT directory containing individual animal folders
fprintf('Please select the PARENT directory containing the animal subfolders.\n');
parentDir = uigetdir(pwd, 'Select Parent Directory');
if parentDir == 0
    error('No directory selected.');
end

% --- CREATE CENTRAL RESULTS FOLDER ---
dateStr = datestr(now, 'yyyymmdd'); 
resultsFolderName = ['Cleaned_photometry_', dateStr];
resultsPath = fullfile(parentDir, resultsFolderName);

if ~exist(resultsPath, 'dir')
    mkdir(resultsPath);
    fprintf('Created new results folder: %s\n', resultsPath);
else
    fprintf('Saving to existing results folder: %s\n', resultsPath);
end

% Filter for animal subfolders
allContents = dir(parentDir);
dirFlags = [allContents.isdir];
subFolders = allContents(dirFlags);
subFolders = subFolders(~ismember({subFolders.name}, {'.', '..', resultsFolderName}));

fprintf('Found %d animal folders. Starting batch process...\n\n', length(subFolders));

successCount = 0;
failCount = 0;

%% Batch Loop %%
for k = 1:length(subFolders)
    
    animalID = subFolders(k).name;
    currentSubPath = fullfile(parentDir, animalID);
    
    fprintf('--------------------------------------------------\n');
    fprintf('Processing (%d/%d): %s\n', k, length(subFolders), animalID);
    
    try
        %% 1. Import Data %%
        Data = TDTbin2mat(currentSubPath);
        
        if isempty(Data.streams) || ~isfield(Data.epocs, 'PC0_')
            error('Missing TDT streams or PC0_ (Blackbox TTL) markers.');
        end
        
        % Extract raw streams starting from Sample 1
        raw_G = Data.streams.G__B.data.'; 
        raw_I = Data.streams.IsoB.data.';   
        Fs = Data.streams.G__B.fs;
        Ts_raw = ((0:(numel(raw_G) - 1)) / Fs)'; 
        bb = Data.epocs.PC0_.onset; 

        %% 2. The Cut (Synchronization) %%
        % Slice the data to the Blackbox TTL window BEFORE any math
        % This removes technical initialization surges
        startIdx = find(Ts_raw >= bb(1), 1, 'first');
        endIdx = find(Ts_raw <= bb(end), 1, 'last');
        
        sliced_G = raw_G(startIdx:endIdx);
        sliced_I = raw_I(startIdx:endIdx);
        synced_Ts = Ts_raw(startIdx:endIdx) - bb(1); % Time 0 = Video Start

        %% 3. Preprocessing and Normalization %%
        % Linear detrend to remove slow baseline drift
        detG = detrend(sliced_G);
        detI = detrend(sliced_I);
        
        % Robust Z-score calculation (Median/MAD method)
        zGCaMP = (detG - median(detG)) / mad(detG, 1);
        zIso = (detI - median(detI)) / mad(detI, 1);
        
        % Movement-corrected signal
        synced_Z = zGCaMP - zIso;

        %% 4. Artifact Removal %%
        z_mean = mean(synced_Z);
        z_std = std(synced_Z);
        
        % Thresholds: +4 SD for spikes, -2 SD for drops
        pos_thresh = z_mean + 4 * z_std;
        neg_thresh = z_mean - 2 * z_std;
        
        remove_idx = (synced_Z > pos_thresh) | (synced_Z < neg_thresh);
        cleaned_Z = synced_Z;
        cleaned_Z(remove_idx) = NaN;
        
        % Interpolate missing values and apply smoothing filter
        cleaned_Z = fillmissing(cleaned_Z, 'linear');
        cleaned_Z = medfilt1(cleaned_Z, 5); 

        %% 5. Export Results %%
        % Save MAT file for MATLAB workflows
        export_data = struct();
        export_data.timestamps = synced_Ts;           
        export_data.cleaned_Z = cleaned_Z;            
        export_data.ZGCaMP = zGCaMP;         
        export_data.ZIso = zIso;              
        export_data.sampling_rate = Fs;               
        export_data.session_start_time = bb(1);      
        export_data.session_end_time = bb(end);      
        
        matFileName = [animalID, '_cleaned.mat'];
        save(fullfile(resultsPath, matFileName), 'export_data');
        
        % Save CSV for Python Ground Truth Validation
        csvFileName = [animalID, '_cleaned.csv'];
        photometry_table = table(synced_Ts, cleaned_Z, zGCaMP, zIso, ...
            'VariableNames', {'time_sec', 'cleaned_Z', 'ZGCaMP', 'ZIso'});
        writetable(photometry_table, fullfile(resultsPath, csvFileName));
        
        fprintf('SUCCESS: Saved %s\n', csvFileName);
        successCount = successCount + 1;

    catch ME
        fprintf('ERROR processing %s: %s\n', animalID, ME.message);
        failCount = failCount + 1;
    end
end

%% Final Summary %%
fprintf('--------------------------------------------------\n');
fprintf('Batch Processing Complete.\n');
fprintf('Results Location: %s\n', resultsPath);
fprintf('Total Successful: %d | Total Failed: %d\n', successCount, failCount);