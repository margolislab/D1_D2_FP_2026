%% Final Statistical Recovery Script
% This script reads your auROC Excel file and generates manuscript-ready stats.
% It ignores any notes or extra text you have added to the bottom of the sheet.

clear; clc; close all;

%% 1. SELECT DATA FILE
[file, path] = uigetfile('*.xlsx', 'Select your Analysis Excel File');
if isequal(file,0); return; end
data = readtable(fullfile(path, file));

%% 2. DATA CLEANING (Remove notes/footer text)
% Convert key columns to strings and trim whitespace to ensure matches work
data.Condition = strtrim(string(data.Condition));
data.Timepoint = strtrim(string(data.Timepoint));

% Remove any rows where Condition or Timepoint is "missing" or empty
% This effectively deletes your notes at the bottom of the Excel sheet
valid_rows = ~ismissing(data.Condition) & data.Condition ~= "" & ...
             ~ismissing(data.Timepoint) & data.Timepoint ~= "";
data = data(valid_rows, :);

% Ensure auROC is numeric (in case text in the column forced it to a cell/string)
if ~isnumeric(data.auROC)
    data.auROC = str2double(string(data.auROC));
end

% Remove any rows where auROC is not a number (NaN)
data(isnan(data.auROC), :) = [];

%% 3. DEFINE HEDGES' G FUNCTION
% Formula includes the correction factor for small sample sizes
calc_g = @(x1, x2) (mean(x1) - mean(x2)) / ...
    sqrt(((length(x1)-1)*var(x1) + (length(x2)-1)*var(x2)) / (length(x1)+length(x2)-2)) * ...
    (1 - (3 / (4*(length(x1)+length(x2)) - 9)));

%% 4. STATISTICAL CALCULATION & REPORTING
groups = {'Sham', 'SNI'};
fprintf('\n======================================================\n');
fprintf('   FULL STATISTICAL REPORT: Baseline2 vs Week3\n');
fprintf('======================================================\n');

for i = 1:length(groups)
    cond = groups{i};
    
    % Extract auROC values for this specific condition
    bl_vals  = data.auROC(strcmpi(data.Condition, cond) & strcmpi(data.Timepoint, 'Baseline2'));
    w3_vals = data.auROC(strcmpi(data.Condition, cond) & strcmpi(data.Timepoint, 'Week3'));
    
    n_bl = length(bl_vals);
    n_w3 = length(w3_vals);
    
    % Only run if we have enough data points
    if n_bl > 1 && n_w3 > 1
        % --- Kolmogorov-Smirnov Test (Distribution Shift) ---
        [~, p_ks, D] = kstest2(bl_vals, w3_vals);
        
        % --- Variance Test (F-test for heterogeneity) ---
        % The 4th output 'stats_var' contains df and the F-statistic
        [~, p_var, ~, stats_var] = vartest2(bl_vals, w3_vals);
        
        % --- Hedges' g (Effect Size) ---
        g = calc_g(bl_vals, w3_vals);
        
        % --- Print Results in Manuscript-Ready Format ---
        fprintf('\n>> CONDITION: %s\n', upper(cond));
        fprintf('   Sample Sizes:  n(BL)=%d, n(W3)=%d\n', n_bl, n_w3);
        
        % KS Result
        fprintf('   KS Test:       D(%d,%d) = %.4f, p = %.4f\n', n_bl, n_w3, D, p_ks);
        
        % Variance Result
        fprintf('   Var Test:      F(%d,%d) = %.4f, p = %.4f\n', stats_var.df1, stats_var.df2, stats_var.fstat, p_var);
        
        % Effect Size Result
        fprintf('   Hedges'' g:     g = %.4f\n', g);
        
        % Combined string for quick copy-pasting into Word
        fprintf('   [Summary]: (KS-D=%.3f, p=%.3f; Var-F=%.2f, p=%.3f; g=%.2f)\n', D, p_ks, stats_var.fstat, p_var, g);
    else
        fprintf('\n>> CONDITION: %s\n', upper(cond));
        fprintf('   Insufficient data for comparison (n_bl=%d, n_w3=%d)\n', n_bl, n_w3);
    end
end
fprintf('======================================================\n');