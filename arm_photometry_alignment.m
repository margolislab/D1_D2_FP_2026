%% ARM Mechanostimulation - Photometry TTL Alignment & Metrics
% Detrends, Z-scores, and aligns TDT photometry data to ARM TTL pulses.
% Calculates sliding-window AUROC and calcium decay kinetics (Tau, slopes).

clear; clc; close all;
warning('off', 'MATLAB:print:ContentTypeImageSuggested');

%% 1. CONFIGURATION & SETUP
% Define Paths (Update these if moving the repository)
SDKpath = 'C:\Users\Arlene\Documents\Margolis Lab\Fiber photometry\FP codes\TDTMatlabSDK-master';
aurocpath = 'C:\Users\Arlene\Documents\Margolis Lab\Fiber photometry\FP codes\AUC';
addpath(genpath(SDKpath));
addpath(genpath(aurocpath));

% Select Data Directory
target_path = uigetdir(pwd, 'Select the TDT Data Folder');
if target_path == 0, return; end
cd(target_path);

% User Prompts
EventsFlag = str2double(cell2mat(inputdlg('Does this file contain event marks? If yes, input 1: ')));
if isnan(EventsFlag), EventsFlag = 0; end

%% 2. IMPORT DATA & DETREND
fprintf('Importing TDT data...\n');
Data = TDTbin2mat(target_path);

% Start from sample 8000 to bypass initialization noise and warm-up transients
start_idx = 8000;
GCaMP = Data.streams.G__B.data(start_idx:end).'; 
Iso   = Data.streams.IsoB.data(start_idx:end).'; 
Fs    = Data.streams.IsoB.fs; 
Ts    = ((start_idx:(start_idx + numel(GCaMP) - 1)) / Fs)'; 

% Line of Best Fit Detrending
fprintf('Detrending signals...\n');
list = {GCaMP, Iso};
DetrendSignals = cell(1, 2);
for i = 1:length(list)
    BestFitG = polyfit(Ts, list{i}, 2);
    TLG = (BestFitG(1)*(Ts.^2)) + (BestFitG(2)*Ts) + BestFitG(3);
    DetrendSignals{i} = list{i} - TLG;
end

%% 3. Z-SCORING & ARTIFACT REMOVAL
% Robust MAD Z-scoring
medianBaseG = median(DetrendSignals{1});
madBaseG    = mean(mad(DetrendSignals{1}));
medianBaseI = median(DetrendSignals{2});
madBaseI    = mean(mad(DetrendSignals{2}));

ZmadGCaMP = (DetrendSignals{1} - medianBaseG) / madBaseG;
ZmadIso   = (DetrendSignals{2} - medianBaseI) / madBaseI;

% Movement-corrected signal
Z = ZmadGCaMP - ZmadIso;

% Artifact Removal (+4 SD / -2 SD thresholds)
z_mean = mean(Z);
z_std  = std(Z);
remove_idx = (Z > (z_mean + 4 * z_std)) | (Z < (z_mean - 2 * z_std));

cleaned_Z = Z;
cleaned_Z(remove_idx) = NaN;
cleaned_Z = fillmissing(cleaned_Z, 'linear');
cleaned_Z = medfilt1(cleaned_Z, 5);
Z = cleaned_Z; % Overwrite with cleaned signal

%% 4. TTL ALIGNMENT & EXTRACTION
if EventsFlag == 1
    fprintf('Aligning to ARM TTL events...\n');
    Events2 = Data.epocs.PC1_.onset; % Adjust TTL port if necessary
    
    Seconds = 10;
    Before = round(Seconds * Fs);
    After = round(Seconds * Fs);
    
    % Preallocate alignment matrices
    num_events = length(Events2) - 2; 
    ZAlign2 = NaN(num_events, Before + After + 1);
    ZAlign2_sp = NaN(num_events, After);
    
    row_idx = 1;
    for i = 2:length(Events2)-1
        [~, a] = min(abs(Ts - Events2(i)));
        
        % Stimulus Window
        if a - Before > 0 && a + After <= length(Z)
            ZAlign2(row_idx, :) = Z(a - Before : a + After);
        end
        
        % Spontaneous Window (10s to 20s post-event)
        if a + After*2 <= length(Z)
            ZAlign2_sp(row_idx, :) = Z(a + After + 1 : a + After*2);
        end
        row_idx = row_idx + 1;
    end
    
    % Group Means & Exports
    ZAlignMean2 = nanmean(ZAlign2, 1);
    ZAlignMean2_sp = nanmean(ZAlign2_sp, 1);
    
    AnimalTrace = ZAlignMean2;  
    AnimalTrials = ZAlign2;     
    time_vector = linspace(-Seconds, Seconds, Before + After + 1);
    save('Animal_Traces.mat', 'AnimalTrace', 'AnimalTrials', 'time_vector');
    
    % --- PLOTTING ---
    time_vector_sp = linspace(Seconds, Seconds*2, size(ZAlignMean2_sp, 2));
    f1 = figure('Name', 'Signal Alignment', 'Color', 'w');
    
    subplot(2,2,1); imagesc(time_vector, 1:size(ZAlign2,1), ZAlign2);
    caxis([-0.5 3]); xline(0, 'r', 'LineWidth', 1); title('All Trials Raster');
    
    subplot(2,2,3); plot(time_vector, ZAlignMean2);
    xlim([-Seconds, Seconds]); ylim([-3, 4]); xline(0, 'r'); title('Mean Z-Scored Signal');
    
    subplot(2,2,2); imagesc(time_vector_sp, 1:size(ZAlign2_sp,1), ZAlign2_sp);
    caxis([-0.5 3]); title('Spontaneous Raster');
    
    subplot(2,2,4); plot(time_vector_sp, ZAlignMean2_sp);
    xlim([Seconds, Seconds*2]); ylim([-3, 4]); title('Spontaneous Signal');
    
    saveas(f1, 'ZAlign_summary.svg');
    
    %% 5. SLIDING WINDOW AUROC ANALYSIS
    fprintf('Calculating Sliding Window AUROC...\n');
    SlidingSteps = round(0.05 * Fs); % 50 ms slide
    Windowsize   = round(0.2 * Fs);  % 200 ms integration window
    StartValue   = Before + round(0.1 * Fs); 
    EndValue     = Before + round(2.0 * Fs); 
    InititalBasalWindow = Before - round(1 * Fs); 
    EndBasalWindow = Before; 
    
    WinTime = []; AreaUC = []; confidenceP = [];
    Wincounter = 1;
    
    for yy = StartValue:SlidingSteps:EndValue
        if yy + Windowsize > size(ZAlign2, 2), break; end
        
        WinTime(Wincounter, 1) = (yy - Before) / Fs;  
        y_basal = mean(ZAlign2(:, InititalBasalWindow:EndBasalWindow), 2);
        y_event = mean(ZAlign2(:, yy:yy+Windowsize), 2);
        
        y_comb = [y_basal; y_event]; 
        categories2 = [zeros(size(y_basal,1),1); ones(size(y_event,1),1)]; 
        
        [A, ~] = auc([categories2, y_comb], 0.05, 'hanley');
        AreaUC(1, Wincounter) = A; 
        confidenceP(1, Wincounter) = auc_bootstrap([categories2, y_comb]); 
        
        Wincounter = Wincounter + 1; 
    end
    
    f2 = figure('Name', 'AUROC', 'Color', 'w');
    plot(WinTime, AreaUC, 'LineWidth', 2);
    xlabel('Time from stimulus (s)'); ylabel('AUROC'); title('AUROC Over Time'); grid on;
    saveas(f2, 'AUROC_plot.svg');
    
    %% 6. METRICS EXTRACTION (AUROC & Calcium Kinetics)
    fprintf('Extracting Kinetics and Peak Metrics...\n');
    % AUROC Metrics
    [peak_auc, peak_idx] = max(AreaUC);
    time_of_peak_auc = WinTime(peak_idx);
    
    sig_idx = find(confidenceP < 0.05, 1, 'first');
    latency_to_sig_auc = ifelse(isempty(sig_idx), NaN, WinTime(sig_idx));
    
    above_thresh = AreaUC > 0.6;
    duration_above_thresh = sum(above_thresh) * SlidingSteps / Fs; 
    auc_integral = trapz(WinTime, AreaUC);
    
    writetable(table(peak_auc, time_of_peak_auc, latency_to_sig_auc, duration_above_thresh, auc_integral), 'AUROC_metrics.csv');
    
    % Trace-Based Calcium Metrics
    zero_index = find(time_vector >= 0, 1, 'first');
    post_trace = ZAlignMean2(zero_index:end);
    post_time  = time_vector(zero_index:end);
    
    [max_amp, rel_peak_idx] = max(post_trace);
    peak_latency = post_time(rel_peak_idx);
    trace_auc = trapz(post_time, post_trace);
    
    threshold = 2;
    above_thresh_idx = find(post_trace > threshold);
    peak_duration = ifelse(isempty(above_thresh_idx), 0, (above_thresh_idx(end) - above_thresh_idx(1)) / Fs);
    
    % Decay Metrics
    decay_trace = post_trace(rel_peak_idx:end);
    xData_s = (0:(length(decay_trace)-1))' ./ Fs;
    
    % Exponential Fit
    try
        ft = fittype('A*exp(-x/tau)+C', 'independent','x','coefficients',{'A','tau','C'});
        opts = fitoptions(ft);
        opts.StartPoint = [double(max(decay_trace)), 1, double(min(decay_trace))]; 
        opts.Lower = [0, 0, -Inf];
        [fitresult, gof] = fit(xData_s, double(decay_trace(:)), ft, opts); 
        decay_tau = fitresult.tau; fit_A = fitresult.A; fit_C = fitresult.C; Rsq = gof.rsquare;
    catch
        fitresult = []; decay_tau = NaN; fit_A = NaN; fit_C = NaN; Rsq = NaN;
    end
    
    % Half Decay
    if ~isempty(fitresult) && ~isnan(fit_A)
        half_value = fit_C + 0.5 * (decay_trace(1) - fit_C);
    else
        half_value = 0.5 * decay_trace(1);
    end
    half_idx = find(decay_trace <= half_value, 1, 'first');
    decay_half_time = ifelse(isempty(half_idx), NaN, (half_idx - 1) / Fs);
    
    % 90-70 Robust Slope
    idx_90 = find(decay_trace <= decay_trace(1) * 0.9, 1, 'first');
    idx_70 = find(decay_trace <= decay_trace(1) * 0.7, 1, 'first');
    if ~isempty(idx_90) && ~isempty(idx_70) && idx_70 > idx_90
        slope_90_70 = (decay_trace(idx_90) - decay_trace(idx_70)) / ((idx_70 - idx_90) / Fs);
    else
        slope_90_70 = NaN;
    end
    
    writetable(table(max_amp, trace_auc, peak_latency, peak_duration, decay_tau, decay_half_time, slope_90_70), 'ZScore_CalciumTrace_Metrics.csv');
    
    % Plot Final Fit
    f3 = figure('Name', 'Decay Fit', 'Color', 'w'); hold on;
    plot(xData_s, decay_trace, 'ko', 'MarkerSize', 4);
    if ~isempty(fitresult)
        plot(xData_s, fit_A * exp(-xData_s / decay_tau) + fit_C, 'r-', 'LineWidth', 2);
        yline(fit_C, 'b--');
        title(sprintf('Decay Fit: \\tau = %.3f s, R^2 = %.3f', decay_tau, Rsq));
    end
    xlabel('Time from Peak (s)'); ylabel('Z-Score'); grid on;
    saveas(f3, 'Decay_Fit.svg');
    fprintf('Processing Complete. Metrics saved to target folder.\n');
else
    figure('Color', 'w'); plot(Z); title('Cleaned Z-Score Trace (No Events)');
end

% Inline helper function for ternary operations
function out = ifelse(cond, true_val, false_val)
    if cond, out = true_val; else, out = false_val; end
end