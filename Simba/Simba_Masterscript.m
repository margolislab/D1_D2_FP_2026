%% Master SimBA-Photometry Pipeline: PSTH, auROC, and Manuscript Stats
clear all; clc; close all;
rng(1); % For reproducibility
warning('off', 'MATLAB:print:ContentTypeImageSuggested');

%% 1. CONFIGURATION & PARAMETERS
target_genotype   = 'A2A';       % 'D1' or 'A2A'
target_behavior   = 'locomotion'; % 'grooming', 'rearing', or 'locomotion'
target_timepoints = {'Baseline2', 'Week12'};
target_condition  = {'SNI', 'Sham'}; 

% Bout & PSTH Windows
min_bout_dur  = 1.0;  
max_bout_dur  = 20.0; 
t_pre         = 2.0;   
t_post        = 15.0; 

% auROC Permutation Settings
num_shuffles  = 10000; 
min_shift_sec = 20;
tail_buffer   = 0.5; 

% Export Path
base_export_path = 'E:\Behavior\Blackbox_analysis\simba_DLS_2\Simba_aligned_photometry';

%% 2. COLOR SETUP (Fake Transparency Trick for Illustrator)
if strcmpi(target_genotype, 'D1')
    c_BL_mean = [0.99, 0.57, 0.45]; c_BL_sem = [0.99, 0.73, 0.63]; % Reds
    c_W12_mean = [0.65, 0.06, 0.08]; c_W12_sem = [0.87, 0.18, 0.15]; 
else % A2A
    c_BL_mean = [0.62, 0.79, 0.88]; c_BL_sem = [0.78, 0.86, 0.94]; % Blues
    c_W12_mean = [0.03, 0.32, 0.61]; c_W12_sem = [0.19, 0.51, 0.74]; 
end
sham_BL_m = [0.74, 0.74, 0.74]; sham_BL_s = [0.85, 0.85, 0.85]; % Grays
sham_W12_m = [0.15, 0.15, 0.15]; sham_W12_s = [0.39, 0.39, 0.39]; % Black

%% 3. LOAD DATA & SETUP DIRECTORIES
[simba_file, simba_path] = uigetfile({'*.xlsx';'*.csv'}, 'Select SimBA Summary');
if isequal(simba_file,0), return; end
simba_data = readtable(fullfile(simba_path, simba_file));

root_dir = uigetdir(pwd, 'Select folder containing the .mat photometry files');
if isequal(root_dir,0), return; end

% Setup Export Directory
master_folder_name = ['auROC_PSTH_', target_behavior, '_', target_genotype];
export_dir = fullfile(base_export_path, master_folder_name, 'Baseline_vs_Week12');
if ~exist(export_dir, 'dir'), mkdir(export_dir); end

% Filter Behavior Table
idx = ismember(strtrim(string(simba_data.Condition)), target_condition) & ...
      ismember(strtrim(string(simba_data.Timepoint)), target_timepoints) & ... 
      strcmpi(strtrim(string(simba_data.Genotype)), target_genotype) & ...
      strcmpi(strtrim(string(simba_data.Event)), target_behavior) & ...
      simba_data.Bout_time >= min_bout_dur & ...
      simba_data.Bout_time <= max_bout_dur;
filtered_data = simba_data(idx, :);
unique_combos = unique(filtered_data(:, {'Subject_ID', 'Timepoint'}), 'rows');

% Initialize Buckets
roc_results = table(); 
sham_BL_snippets = []; sham_W12_snippets = [];
sni_BL_snippets  = []; sni_W12_snippets  = [];

%% 4. PROCESSING LOOP (Extraction & auROC)
for i = 1:height(unique_combos)
    this_subject = char(strtrim(string(unique_combos.Subject_ID(i))));
    this_tp      = char(strtrim(string(unique_combos.Timepoint(i))));
    
    fprintf('Processing %d/%d: %s [%s]\n', i, height(unique_combos), this_subject, this_tp);
    
    % Find specific .mat file
    mat_search = dir(fullfile(root_dir, '**', ['*', this_subject, '*.mat']));
    match_idx = contains({mat_search.folder}, this_tp, 'IgnoreCase', true) | ...
                contains({mat_search.name}, this_tp, 'IgnoreCase', true);
    if ~any(match_idx), continue; end 
    
    target_file = mat_search(find(match_idx, 1));
    load(fullfile(target_file.folder, target_file.name)); 
    Z = double(export_data.cleaned_Z); ts = export_data.timestamps; Fs = export_data.sampling_rate;
    
    subject_bouts = filtered_data(strcmpi(filtered_data.Subject_ID, this_subject) & ...
                                  strcmpi(filtered_data.Timepoint, this_tp), :);
    this_cond = string(unique(subject_bouts.Condition));
    num_bouts = height(subject_bouts);
    if num_bouts == 0, continue; end
    
    % --- auROC Calculation ---
    bout_mask = zeros(size(Z));
    for b = 1:num_bouts
        s_idx = round(subject_bouts.Start_time(b) * Fs); 
        e_idx = round((subject_bouts.End_time(b) + tail_buffer) * Fs); 
        s_idx = max(1, s_idx); e_idx = min(length(Z), e_idx);
        bout_mask(s_idx:e_idx) = 1;
    end
    
    sig_in = Z(bout_mask == 1); sig_out = Z(bout_mask == 0);
    [~, ~, stats] = ranksum(sig_in, sig_out);
    U = stats.ranksum - (length(sig_in)*(length(sig_in)+1)/2);
    obs_auROC = U / (length(sig_in) * length(sig_out));
    
    fprintf('  -> Running %d permutations...', num_shuffles);
    shuff_auCs = zeros(num_shuffles, 1);
    for s = 1:num_shuffles
        if mod(s, 1000) == 0, fprintf('.'); end
        shift = randi([round(min_shift_sec*Fs), length(Z)-round(min_shift_sec*Fs)]);
        s_mask = circshift(bout_mask, shift);
        s_in = Z(s_mask==1); s_out = Z(s_mask==0);
        [~, ~, s_stats] = ranksum(s_in, s_out);
        s_U = s_stats.ranksum - (length(s_in)*(length(s_in)+1)/2);
        shuff_auCs(s) = s_U / (length(s_in) * length(s_out));
    end
    p_val = sum(abs(shuff_auCs - 0.5) >= abs(obs_auROC - 0.5)) / num_shuffles;
    fprintf(' Done.\n');
    
    roc_results = [roc_results; table(string(this_subject), string(this_cond), string(this_tp), ...
        num_bouts, obs_auROC, p_val, 'VariableNames', {'SubjectID', 'Condition', 'Timepoint', 'Bouts', 'auROC', 'pVal'})];
        
    % --- PSTH Snippet Extraction ---
    expected_length = round(t_pre * Fs) + round(t_post * Fs) + 1;
    sub_snippets = [];
    for b = 1:num_bouts
        [~, onset_idx] = min(abs(ts - subject_bouts.Start_time(b)));
        s_idx = onset_idx - round(t_pre * Fs);
        e_idx = s_idx + expected_length - 1;
        if s_idx > 0 && e_idx <= length(Z)
            snip = Z(s_idx:e_idx); 
            snip = snip - mean(snip(1:round(t_pre*Fs))); 
            sub_snippets = [sub_snippets; snip(:)'];
        end
    end
    
    if isempty(sub_snippets), continue; end
    if strcmpi(this_cond, 'Sham')
        if strcmpi(this_tp, 'Baseline2'), sham_BL_snippets = [sham_BL_snippets; sub_snippets];
        elseif strcmpi(this_tp, 'Week12'), sham_W12_snippets = [sham_W12_snippets; sub_snippets]; end
    else 
        if strcmpi(this_tp, 'Baseline2'), sni_BL_snippets = [sni_BL_snippets; sub_snippets];
        elseif strcmpi(this_tp, 'Week12'), sni_W12_snippets = [sni_W12_snippets; sub_snippets]; end
    end
end

%% 5. FIGURE 1: DOWNSAMPLED VECTOR PSTH
t_ax_raw = (0:expected_length-1)/Fs - t_pre;
ds_factor = max(1, round(Fs / 100)); % Downsample to ~100Hz for clean vector paths
time_axis = t_ax_raw(1:ds_factor:end);

f_psth = figure('Color', 'w', 'Renderer', 'painters', 'Position', [100 100 1000 400]);
% Sham Panel
subplot(1,2,1); hold on;
if ~isempty(sham_BL_snippets)
    data_bl = sham_BL_snippets(:, 1:ds_factor:end);
    m = mean(data_bl, 1); s = std(data_bl, 0, 1)/sqrt(size(data_bl,1));
    patch([time_axis, fliplr(time_axis)], [m+s, fliplr(m-s)], sham_BL_s, 'EdgeColor', 'none', 'FaceAlpha', 1);
    plot(time_axis, m, 'Color', sham_BL_m, 'LineWidth', 2);
end
if ~isempty(sham_W12_snippets)
    data_w12 = sham_W12_snippets(:, 1:ds_factor:end);
    m = mean(data_w12, 1); s = std(data_w12, 0, 1)/sqrt(size(data_w12,1));
    patch([time_axis, fliplr(time_axis)], [m+s, fliplr(m-s)], sham_W12_s, 'EdgeColor', 'none', 'FaceAlpha', 1);
    plot(time_axis, m, 'Color', sham_W12_m, 'LineWidth', 2);
end
line([0 0], ylim, 'Color', 'k', 'LineStyle', '--'); title('Sham: BL vs W12'); set(gca, 'TickDir', 'out', 'Box', 'off');

% SNI Panel
subplot(1,2,2); hold on;
if ~isempty(sni_BL_snippets)
    data_bl = sni_BL_snippets(:, 1:ds_factor:end);
    m = mean(data_bl, 1); s = std(data_bl, 0, 1)/sqrt(size(data_bl,1));
    patch([time_axis, fliplr(time_axis)], [m+s, fliplr(m-s)], c_BL_sem, 'EdgeColor', 'none', 'FaceAlpha', 1);
    plot(time_axis, m, 'Color', c_BL_mean, 'LineWidth', 2);
end
if ~isempty(sni_W12_snippets)
    data_w12 = sni_W12_snippets(:, 1:ds_factor:end);
    m = mean(data_w12, 1); s = std(data_w12, 0, 1)/sqrt(size(data_w12,1));
    patch([time_axis, fliplr(time_axis)], [m+s, fliplr(m-s)], c_W12_sem, 'EdgeColor', 'none', 'FaceAlpha', 1);
    plot(time_axis, m, 'Color', c_W12_mean, 'LineWidth', 2);
end
line([0 0], ylim, 'Color', 'k', 'LineStyle', '--'); title(['SNI ' target_genotype ': BL vs W12']); set(gca, 'TickDir', 'out', 'Box', 'off');

%% 6. FIGURE 2: AUROC HISTOGRAMS
edges = 0:0.05:1;
h_sham_bl  = roc_results.auROC(roc_results.Condition == "Sham" & roc_results.Timepoint == "Baseline2");
h_sham_w12 = roc_results.auROC(roc_results.Condition == "Sham" & roc_results.Timepoint == "Week12");
h_sni_bl   = roc_results.auROC(roc_results.Condition == "SNI"  & roc_results.Timepoint == "Baseline2");
h_sni_w12  = roc_results.auROC(roc_results.Condition == "SNI"  & roc_results.Timepoint == "Week12");

f_hist = figure('Color', 'w', 'Position', [150 150 1000 400]);
subplot(1,2,1); hold on;
if ~isempty(h_sham_bl), histogram(h_sham_bl, edges, 'Normalization', 'probability', 'FaceColor', sham_BL_m, 'FaceAlpha', 0.4); end
if ~isempty(h_sham_w12), histogram(h_sham_w12, edges, 'Normalization', 'probability', 'FaceColor', sham_W12_m, 'FaceAlpha', 0.6); end
xline(0.5, '--k', 'Chance'); title('Sham auROC Distribution');

subplot(1,2,2); hold on;
if ~isempty(h_sni_bl), histogram(h_sni_bl, edges, 'Normalization', 'probability', 'FaceColor', c_BL_mean, 'FaceAlpha', 0.4); end
if ~isempty(h_sni_w12), histogram(h_sni_w12, edges, 'Normalization', 'probability', 'FaceColor', c_W12_mean, 'FaceAlpha', 0.6); end
xline(0.5, '--k', 'Chance'); title(['SNI (' target_genotype ') auROC Distribution']);

%% 7. EXPORT DATA & FIGURES
% Export Figures
exportgraphics(f_psth, fullfile(export_dir, 'Illustrator_PSTH_Vector.pdf'), 'ContentType', 'vector');
print(f_psth, fullfile(export_dir, 'Illustrator_PSTH_Vector.svg'), '-dsvg', '-vector');
exportgraphics(f_hist, fullfile(export_dir, 'Histogram_Distributions.pdf'), 'ContentType', 'vector');

% Export Stats & Prism Tables
writetable(roc_results, fullfile(export_dir, 'Final_Stats_Individual_Bouts.xlsx'));

p_sham_bl  = histcounts(h_sham_bl, edges, 'Normalization', 'probability')';
p_sham_w12 = histcounts(h_sham_w12, edges, 'Normalization', 'probability')';
p_sni_bl   = histcounts(h_sni_bl, edges, 'Normalization', 'probability')';
p_sni_w12  = histcounts(h_sni_w12, edges, 'Normalization', 'probability')';
hist_table = table(edges(1:end-1)', p_sham_bl, p_sham_w12, p_sni_bl, p_sni_w12, ...
    'VariableNames', {'Bin_Start', 'Sham_BL_Prob', 'Sham_W12_Prob', 'SNI_BL_Prob', 'SNI_W12_Prob'});
writetable(hist_table, fullfile(export_dir, 'Histogram_Data_for_Prism.xlsx'));

%% 8. IN-MEMORY STATISTICAL RECOVERY (Manuscript Output)
% Hedges' g inline function
calc_g = @(x1, x2) (mean(x1) - mean(x2)) / ...
    sqrt(((length(x1)-1)*var(x1) + (length(x2)-1)*var(x2)) / (length(x1)+length(x2)-2)) * ...
    (1 - (3 / (4*(length(x1)+length(x2)) - 9)));

fprintf('\n======================================================\n');
fprintf('   FULL STATISTICAL REPORT: Baseline2 vs Week12\n');
fprintf('======================================================\n');

groups = {'Sham', 'SNI'};
for i = 1:length(groups)
    cond = groups{i};
    bl_vals  = roc_results.auROC(roc_results.Condition == cond & roc_results.Timepoint == "Baseline2");
    w12_vals = roc_results.auROC(roc_results.Condition == cond & roc_results.Timepoint == "Week12");
    
    bl_vals(isnan(bl_vals)) = [];
    w12_vals(isnan(w12_vals)) = [];
    
    n_bl = length(bl_vals);
    n_w12 = length(w12_vals);
    
    fprintf('\n>> CONDITION: %s\n', upper(cond));
    if n_bl > 1 && n_w12 > 1
        [~, p_ks, D] = kstest2(bl_vals, w12_vals);
        [~, p_var, ~, stats_var] = vartest2(bl_vals, w12_vals);
        g = calc_g(bl_vals, w12_vals);
        
        fprintf('   Sample Sizes:  n(BL)=%d, n(W12)=%d\n', n_bl, n_w12);
        fprintf('   KS Test:       D(%d,%d) = %.4f, p = %.4f\n', n_bl, n_w12, D, p_ks);
        fprintf('   Var Test:      F(%d,%d) = %.4f, p = %.4f\n', stats_var.df1, stats_var.df2, stats_var.fstat, p_var);
        fprintf('   Hedges'' g:     g = %.4f\n', g);
        fprintf('   [Summary]: (KS-D=%.3f, p=%.3f; Var-F=%.2f, p=%.3f; g=%.2f)\n', D, p_ks, stats_var.fstat, p_var, g);
    else
        fprintf('   Insufficient data for comparison (n_bl=%d, n_w12=%d)\n', n_bl, n_w12);
    end
end
fprintf('======================================================\n');
fprintf('Export Complete: Saved to %s\n', export_dir);