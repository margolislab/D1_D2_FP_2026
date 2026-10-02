%% VAME-Photometry: Sham vs. SNI Comparison Suite
clear all; clc; close all;

%% 1. TARGET SELECTION
target_timepoint = 'Week3';    % Choose the timepoint to compare
target_genotype  = 'A2A';      % 'D1' or 'A2A'
target_community = 3;          % 5=Locomotion, 3=Guarding
exclude_subjects = {'woody', 'shaggy'}; 

% --- BOUT FILTERS & WINDOWS ---
min_bout_s = 1.5;    
purity_s   = 2.0;    
t_pre      = 5;      
t_post     = 8;      % Visual window set to 4s as requested
math_win_s = 8;      % Statistical analysis window set to 4s

% --- SYSTEM SETTINGS ---
video_fps  = 45;
photo_fs   = 1017.25;
vame_offset_frames = 22; 
save_dir = 'E:\Behavior\Blackbox_analysis\VAME_Comparison_Results';
if ~exist(save_dir, 'dir'), mkdir(save_dir); end

%% 2. PATHS & METADATA
meta_path = 'E:\Behavior\Blackbox_analysis\Metadata.xlsx';
vame_root = 'E:\Behavior\VAME_Projects\SNI_Longitudinal_Full-Jul15-2026\results';
photo_root = 'E:\Behavior\Blackbox_analysis\photometry';

meta_data = readtable(meta_path);
% Filter only for Timepoint and Genotype (keep both Sham and SNI)
idx = strcmpi(meta_data.Timepoint, target_timepoint) & ...
      strcmpi(meta_data.Genotype, target_genotype);
group_metadata = meta_data(idx, :);
unique_subjects = unique(group_metadata.Subject_ID);

%% 3. PROCESSING LOOP
% Storage for plotting and stats
data_store = struct('Sham_Traces', [], 'SNI_Traces', [], 'Sham_Stats', [], 'SNI_Stats', [], 'Sham_Corr', [], 'SNI_Corr', []);

for i = 1:length(unique_subjects)
    this_sub = char(unique_subjects{i});
    if ismember(lower(this_sub), lower(exclude_subjects)), continue; end
    
    sub_info = group_metadata(strcmpi(group_metadata.Subject_ID, this_sub), :);
    this_cond = char(sub_info.Condition(1)); % 'Sham' or 'SNI'
    current_sub_bouts = [];
    
    for s = 1:height(sub_info)
        session_name = char(sub_info.session(s));
        
        % Load Data
        photo_file = fullfile(photo_root, target_timepoint, [this_sub, '*_cleaned.mat']);
        d = dir(photo_file); if isempty(d), d = dir(fullfile(photo_root, '**', [this_sub, '*_cleaned.mat'])); end
        if isempty(d), continue; end
        load(fullfile(d(1).folder, d(1).name)); 
        z_signal = export_data.cleaned_Z;
        
        vame_file = fullfile(vame_root, session_name, 'VAME', 'hmm-15', 'community', ['community_label_', session_name, '.mat']);
        if ~exist(vame_file, 'file'), continue; end
        v_data = load(vame_file); labels = v_data.community_labels;
        
        onsets = find(labels(2:end) == target_community & labels(1:end-1) ~= target_community) + 1;
        ratio = photo_fs / video_fps;
        
        for b = 1:length(onsets)
            curr = onsets(b);
            if curr <= round(purity_s*video_fps), continue; end
            if any(labels(curr-round(purity_s*video_fps):curr-1) == target_community), continue; end
            lat = labels(curr:end); dur_f = find(lat ~= target_community, 1, 'first') - 1;
            if isempty(dur_f), dur_f = length(lat); end
            bout_dur_s = dur_f / video_fps;
            if bout_dur_s < min_bout_s, continue; end
            
            center = round((curr + vame_offset_frames) * ratio);
            start_idx = center - round(t_pre * photo_fs);
            end_idx   = center + round(t_post * photo_fs);
            
            if start_idx > 0 && end_idx <= length(z_signal)
                snippet = z_signal(start_idx : end_idx);
                base_idx = 1 : round(2*photo_fs); % -5 to -3s
                bout_z = (snippet - nanmean(snippet(base_idx))) / nanstd(snippet(base_idx));
                
                % Store individual bout for Correlation
                resp_idx = round(t_pre*photo_fs) : round((t_pre + math_win_s)*photo_fs);
                bout_peak = max(bout_z(resp_idx));
                if strcmpi(this_cond, 'Sham'), data_store.Sham_Corr = [data_store.Sham_Corr; bout_dur_s, bout_peak];
                else, data_store.SNI_Corr = [data_store.SNI_Corr; bout_dur_s, bout_peak]; end
                
                current_sub_bouts = [current_sub_bouts; bout_z(:)'];
            end
        end
    end
    
    % --- ANIMAL LEVEL STATS ---
    if ~isempty(current_sub_bouts)
        sub_mean = nanmean(current_sub_bouts, 1);
        t_vec = linspace(-t_pre, t_post, length(sub_mean));
        post_idx = t_vec >= 0 & t_vec <= math_win_s;
        base_idx = t_vec >= -5 & t_vec <= -3;
        
        peak_val = max(sub_mean(post_idx));
        auc_val  = trapz(t_vec(post_idx), sub_mean(post_idx));
        [~,~,~,auroc_val] = perfcurve([ones(sum(post_idx),1); zeros(sum(base_idx),1)], [sub_mean(post_idx)'; sub_mean(base_idx)'], 1);
        
        res_row = {this_sub, peak_val, auc_val, auroc_val};
        if strcmpi(this_cond, 'Sham')
            data_store.Sham_Traces = [data_store.Sham_Traces; sub_mean];
            data_store.Sham_Stats = [data_store.Sham_Stats; res_row];
        else
            data_store.SNI_Traces = [data_store.SNI_Traces; sub_mean];
            data_store.SNI_Stats = [data_store.SNI_Stats; res_row];
        end
    end
end

%% 4. COLORS & SETTINGS
% Sham is always Black. SNI depends on genotype.
shamCol = [0.2 0.2 0.2]; 
if strcmpi(target_genotype, 'D1'), sniCol = [0.8 0 0]; % Red
else, sniCol = [0 0.45 0.74]; end % Blue

base_fn = sprintf('VAME_C%d_%s_%s_COMPARE', target_community, target_genotype, target_timepoint);

%% 5. FIGURE 1: OVERLAY PETH
f1 = figure('Color','w', 'Units', 'inches', 'Position', [2 2 6 4]); hold on;
groups = {'Sham', 'SNI'};
traces = {data_store.Sham_Traces, data_store.SNI_Traces};
colors = {shamCol, sniCol};

for g = 1:2
    if isempty(traces{g}), continue; end
    m = nanmean(traces{g}, 1);
    s = nanstd(traces{g}, 0, 1) / sqrt(size(traces{g}, 1));
    fill([t_vec, fliplr(t_vec)], [m+s, fliplr(m-s)], colors{g}, 'FaceAlpha', 0.15, 'EdgeColor', 'none');
    plot(t_vec, m, 'Color', colors{g}, 'LineWidth', 2.5, 'DisplayName', groups{g});
end
line([0 0], ylim, 'Color', [0.5 0.5 0.5], 'LineStyle', '--');
xlabel('Time from Onset (s)'); ylabel('Z-Score'); legend('Location','best');
title(sprintf('Community %d | %s | Sham vs SNI', target_community, target_timepoint));
grid on; set(gca, 'Box', 'off');

%% 6. FIGURE 2: SIDE-BY-SIDE STATS
f2 = figure('Color','w','Position', [100 100 1000 400]);
metrics = {'Peak_Z', 'AUC', 'auROC'};
for k = 1:3
    subplot(1,3,k); hold on;
    % Prepare data for boxplot
    d_sham = cell2mat(data_store.Sham_Stats(:, k+1));
    d_sni = cell2mat(data_store.SNI_Stats(:, k+1));
    
    % Use a grouped boxplot approach
    boxplot([d_sham; d_sni], [zeros(size(d_sham)); ones(size(d_sni))], 'Colors', 'k', 'Symbol', '');
    
    % Overlay Sham dots
    scatter(1+(rand(length(d_sham),1)-0.5)*0.1, d_sham, 40, 'filled', 'MarkerFaceColor', shamCol, 'MarkerFaceAlpha', 0.6);
    % Overlay SNI dots
    scatter(2+(rand(length(d_sni),1)-0.5)*0.1, d_sni, 40, 'filled', 'MarkerFaceColor', sniCol, 'MarkerFaceAlpha', 0.6);
    
    set(gca, 'XTick', [1 2], 'XTickLabel', {'Sham', 'SNI'});
    title(metrics{k}); ylabel('Value'); grid on;
end

%% 7. FIGURE 3: DUAL CORRELATION
f3 = figure('Color','w'); hold on;
for g = 1:2
    if g==1, d = data_store.Sham_Corr; c = shamCol; lbl = 'Sham';
    else, d = data_store.SNI_Corr; c = sniCol; lbl = 'SNI'; end
    if isempty(d), continue; end
    scatter(d(:,1), d(:,2), 30, 'filled', 'MarkerFaceColor', c, 'MarkerFaceAlpha', 0.3, 'DisplayName', lbl);
    h = lsline; set(h(1), 'Color', c, 'LineWidth', 2);
end
xlabel('Bout Duration (s)'); ylabel('Peak Z-Score'); title('Bout Correlation'); legend; grid on;

%% 8. SAVE
save_list = {f1, f2, f3}; names = {'_PETH_Compare', '_Stats_Compare', '_Corr_Compare'};
for i = 1:3
    fn = fullfile(save_dir, [base_fn, names{i}]);
    exportgraphics(save_list{i}, [fn, '.pdf'], 'ContentType', 'vector');
    saveas(save_list{i}, [fn, '.svg']);
end