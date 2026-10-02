%% VAME State Analysis: Baseline 2 vs. Week 3 (Integrated Stats & Export)
clear all; clc; close all;

% Reproducibility
rng(1); 
warning('off', 'MATLAB:print:ContentTypeImageSuggested');

%% 1. TARGET SELECTION
target_timepoints = {'Baseline2', 'Week3'}; 
target_condition  = {'SNI', 'Sham'}; 
target_genotype   = 'A2A';      % Change to 'D1' as needed
target_community  = 1;          % Change to 5 for Locomotion

% --- PSTH WINDOWS ---
t_pre = 2.0;   
t_post = 8.0; 

% --- BOUT FILTERS ---
min_bout_s = 1.0; 

% --- AUROC SETTINGS ---
num_shuffles = 10000; 
min_shift_sec = 20;

% --- SYSTEM SETTINGS ---
video_fps = 45;
photo_fs  = 1017.25;
vame_offset_frames = 22; 

% --- PATHS (Updated as per your local structure) ---
base_export_path = 'E:\Behavior\Blackbox_analysis\VAME_Projects\VAME_Photometry_Results\VAME_Baseline_vs_Week3';
meta_path  = 'E:\Behavior\Blackbox_analysis\VAME_Projects\Metadata.xlsx';
vame_root  = 'E:\Behavior\Blackbox_analysis\VAME_Projects\SNI_Longitudinal_Full-Jul15-2026\results';
photo_root = 'E:\Behavior\Blackbox_analysis\photometry';

%% 2. SETUP DIRECTORIES & METADATA
export_dir = fullfile(base_export_path, ['Community_', num2str(target_community)], target_genotype);
if ~exist(export_dir, 'dir'), mkdir(export_dir); end

meta_data = readtable(meta_path);
idx = ismember(strtrim(string(meta_data.Condition)), target_condition) & ...
      ismember(strtrim(string(meta_data.Timepoint)), target_timepoints) & ... 
      strcmpi(strtrim(string(meta_data.Genotype)), target_genotype);
filtered_meta = meta_data(idx, :);

% Fixed table initialization
roc_results = table('Size',[0 5], ...
    'VariableTypes',{'string','string','string','double','double'}, ...
    'VariableNames',{'SubjectID','Condition','Timepoint','auROC','pVal'});

sham_BL_snippets = []; sham_W3_snippets = [];
sni_BL_snippets  = []; sni_W3_snippets  = [];

%% 3. PROCESSING LOOP
processed_count = 0;
for i = 1:height(filtered_meta)
    this_sub = char(filtered_meta.Subject_ID(i));
    this_tp  = char(filtered_meta.Timepoint(i));
    this_cond = char(filtered_meta.Condition(i));
    session_name = char(filtered_meta.session(i));
    
    fprintf('Processing %s | %s [%s]...\n', this_sub, this_tp, this_cond);

    % Load Photometry
    photo_file = fullfile(photo_root, this_tp, [this_sub, '*_cleaned.mat']);
    d = dir(photo_file); if isempty(d), d = dir(fullfile(photo_root, '**', [this_sub, '*_cleaned.mat'])); end
    if isempty(d), continue; end
    
    match_tp = contains({d.folder}, this_tp, 'IgnoreCase', true) | contains({d.name}, this_tp, 'IgnoreCase', true);
    if ~any(match_tp), continue; end
    d = d(find(match_tp, 1));
    
    load(fullfile(d.folder, d.name)); 
    Z = double(export_data.cleaned_Z); Fs = export_data.sampling_rate;
    
    % Load VAME
    vame_file = fullfile(vame_root, session_name, 'VAME', 'hmm-15', 'community', ['community_label_', session_name, '.mat']);
    if ~exist(vame_file, 'file'), continue; end
    v_data = load(vame_file); labels = v_data.community_labels;
    
    ratio = Fs / video_fps;
    state_mask = zeros(size(Z));
    sub_snippets = [];
    
    onsets = find(labels(2:end) == target_community & labels(1:end-1) ~= target_community) + 1;
    
    for b = 1:length(onsets)
        curr = onsets(b);
        lat = labels(curr:end); dur_f = find(lat ~= target_community, 1, 'first') - 1;
        if isempty(dur_f), dur_f = length(lat); end
        if (dur_f / video_fps) < min_bout_s, continue; end
        
        center = round((curr + vame_offset_frames) * ratio);
        bout_end_idx = center + round((dur_f/video_fps) * Fs);
        
        if center > 0 && bout_end_idx <= length(Z)
            state_mask(center : min(length(Z), bout_end_idx)) = 1;
            s_idx = center - round(t_pre * Fs);
            e_idx = center + round(t_post * Fs);
            if s_idx > 0 && e_idx <= length(Z)
                snip = Z(s_idx:e_idx);
                snip = snip - mean(snip(1:round(t_pre*Fs))); 
                sub_snippets = [sub_snippets; snip(:)'];
            end
        end
    end
    
    if isempty(sub_snippets) || sum(state_mask) == 0, continue; end
    
    % auROC & Permutation
    sig_in = Z(state_mask == 1); sig_out = Z(state_mask == 0);
    [~, ~, stats] = ranksum(sig_in, sig_out);
    U = stats.ranksum - (length(sig_in)*(length(sig_in)+1)/2);
    obs_auROC = U / (length(sig_in) * length(sig_out));
    
    shuff_auCs = zeros(num_shuffles, 1);
    for s = 1:num_shuffles
        shift = randi([round(min_shift_sec*Fs), length(Z)-round(min_shift_sec*Fs)]);
        s_mask = circshift(state_mask, shift);
        s_in = Z(s_mask==1); s_out = Z(s_mask==0);
        [~, ~, s_stats] = ranksum(s_in, s_out);
        s_U = s_stats.ranksum - (length(s_in)*(length(s_in)+1)/2);
        shuff_auCs(s) = s_U / (length(s_in) * length(s_out));
    end
    p_val = sum(abs(shuff_auCs - 0.5) >= abs(obs_auROC - 0.5)) / num_shuffles;
    
    roc_results = [roc_results; {string(this_sub), string(this_cond), string(this_tp), obs_auROC, p_val}];
    processed_count = processed_count + 1;

    if strcmpi(this_cond, 'Sham')
        if strcmpi(this_tp, 'Baseline2'), sham_BL_snippets = [sham_BL_snippets; mean(sub_snippets, 1)];
        else, sham_W3_snippets = [sham_W3_snippets; mean(sub_snippets, 1)]; end
    else
        if strcmpi(this_tp, 'Baseline2'), sni_BL_snippets = [sni_BL_snippets; mean(sub_snippets, 1)];
        else, sni_W3_snippets = [sni_W3_snippets; mean(sub_snippets, 1)]; end
    end
end

%% 4. PLOTTING
if processed_count == 0, error('No data found.'); end

if strcmpi(target_genotype, 'D1')
    c_BL_m = [0.99, 0.57, 0.45]; c_W3_m = [0.65, 0.06, 0.08]; 
else
    c_BL_m = [0.62, 0.79, 0.88]; c_W3_m = [0.03, 0.32, 0.61]; 
end
sham_BL_m = [0.74, 0.74, 0.74]; sham_W3_m = [0.15, 0.15, 0.15];

f_final = figure('Color', 'w', 'Position', [50 50 1100 800]);
time_axis = linspace(-t_pre, t_post, size(sham_BL_snippets, 2));

subplot(2,2,1); hold on; title('Sham: BL2 vs Week 3');
if ~isempty(sham_BL_snippets), fill_plot(time_axis, sham_BL_snippets, sham_BL_m, 'BL2'); end
if ~isempty(sham_W3_snippets), fill_plot(time_axis, sham_W3_snippets, sham_W3_m, 'W3'); end
line([0 0], ylim, 'Color', 'k', 'LineStyle', '--'); legend;

subplot(2,2,2); hold on; title(['SNI (', target_genotype, '): BL2 vs Week 3']);
if ~isempty(sni_BL_snippets), fill_plot(time_axis, sni_BL_snippets, c_BL_m, 'BL2'); end
if ~isempty(sni_W3_snippets), fill_plot(time_axis, sni_W3_snippets, c_W3_m, 'W3'); end
line([0 0], ylim, 'Color', 'k', 'LineStyle', '--'); legend;

% Histograms
edges = 0:0.05:1;
subplot(2,2,3); hold on; title('Sham auROC Distribution');
h_sham_bl = roc_results.auROC(roc_results.Condition == "Sham" & roc_results.Timepoint == "Baseline2");
h_sham_w3 = roc_results.auROC(roc_results.Condition == "Sham" & roc_results.Timepoint == "Week3");
if ~isempty(h_sham_bl), histogram(h_sham_bl, edges, 'Normalization', 'probability', 'FaceColor', sham_BL_m, 'FaceAlpha', 0.4); end
if ~isempty(h_sham_w3), histogram(h_sham_w3, edges, 'Normalization', 'probability', 'FaceColor', sham_W3_m, 'FaceAlpha', 0.6); end
xline(0.5, '--k');

subplot(2,2,4); hold on; title(['SNI (', target_genotype, ') auROC Distribution']);
h_sni_bl = roc_results.auROC(roc_results.Condition == "SNI" & roc_results.Timepoint == "Baseline2");
h_sni_w3 = roc_results.auROC(roc_results.Condition == "SNI" & roc_results.Timepoint == "Week3");
if ~isempty(h_sni_bl), histogram(h_sni_bl, edges, 'Normalization', 'probability', 'FaceColor', c_BL_m, 'FaceAlpha', 0.4); end
if ~isempty(h_sni_w3), histogram(h_sni_w3, edges, 'Normalization', 'probability', 'FaceColor', c_W3_m, 'FaceAlpha', 0.6); end
xline(0.5, '--k');

%% 5. INTEGRATED STATS (Your requested block)
fprintf('\n--- auROC STATISTICAL COMPARISONS (BL2 vs W3) ---\n');
calc_g = @(x1, x2) (mean(x1) - mean(x2)) / sqrt(((length(x1)-1)*var(x1) + (length(x2)-1)*var(x2)) / (length(x1)+length(x2)-2)) * (1 - (3 / (4*(length(x1)+length(x2)) - 9)));

if ~isempty(h_sni_bl) && ~isempty(h_sni_w3)
    [~, p_ks] = kstest2(h_sni_bl, h_sni_w3); [~, p_var] = vartest2(h_sni_bl, h_sni_w3); g = calc_g(h_sni_bl, h_sni_w3);
    fprintf('SNI (%s): KS-p=%.4f, Var-p=%.4f, Hedges-g=%.4f\n', target_genotype, p_ks, p_var, g);
end
if ~isempty(h_sham_bl) && ~isempty(h_sham_w3)
    [~, p_ks] = kstest2(h_sham_bl, h_sham_w3); [~, p_var] = vartest2(h_sham_bl, h_sham_w3); g = calc_g(h_sham_bl, h_sham_w3);
    fprintf('Sham: KS-p=%.4f, Var-p=%.4f, Hedges-g=%.4f\n', p_ks, p_var, g);
end

%% 6. EXPORT (Integrated Save & Prism block)
save_fn = sprintf('%s_C%d_Analysis_BL2_vs_W3', target_genotype, target_community);
exportgraphics(f_final, fullfile(export_dir, [save_fn, '.pdf']), 'ContentType', 'vector');
writetable(roc_results, fullfile(export_dir, [save_fn, '_IndividualStats.xlsx']));

p_sbl = histcounts(h_sham_bl, edges, 'Normalization', 'probability')';
p_sw3 = histcounts(h_sham_w3, edges, 'Normalization', 'probability')';
p_nbl = histcounts(h_sni_bl, edges, 'Normalization', 'probability')';
p_nw3 = histcounts(h_sni_w3, edges, 'Normalization', 'probability')';

hist_table = table(edges(1:end-1)', p_sbl, p_sw3, p_nbl, p_nw3, ...
    'VariableNames', {'Bin_Start', 'Sham_BL2', 'Sham_W3', 'SNI_BL2', 'SNI_W3'});
writetable(hist_table, fullfile(export_dir, [save_fn, '_Histogram_Data_Prism.xlsx']));

fprintf('✅ Community %d Analysis Saved.\n', target_community);

% Helper function for fill
function fill_plot(t, data, col, lbl)
    m = mean(data, 1); s = std(data,0,1)/sqrt(size(data,1));
    fill([t, fliplr(t)], [m+s, fliplr(m-s)], col, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    plot(t, m, 'Color', col, 'LineWidth', 2, 'DisplayName', lbl);
end