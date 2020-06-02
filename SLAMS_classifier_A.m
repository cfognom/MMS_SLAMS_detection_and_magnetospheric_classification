classdef SLAMS_classifier_A < handle
    %SLAMS_CLASSIFIER_A Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        SW_classifier
    end
    
    methods
        function obj = SLAMS_classifier_A()
            %SLAMS_FINDER1 Construct an instance of this class
            %   Detailed explanation goes here

        end

        function ts = load_ts(obj, tint, ts)

            if nargin < 3 || ~isfield(ts, 'b')
                ts.b = load_tints_MMS('mms1_fgm_srvy_l2', 'mms1_fgm_b_gse_srvy_l2', tint);
            end
            if ~isfield(ts, 'b_abs')
                ts.b_abs = irf_abs(ts.b);
            end
        end

        function obj = train(obj)

            X = readmatrix(['cache/', 'traindata1.txt']);
            % train_data2 = readmatrix([cache_path, 'traindata2.txt']);

            X = [irf_abs(X(:,2:4), 1), irf_abs(X(:,5:7), 1), X(:,8)];
            % X = [irf_abs(X(:,2:4), 1), X(:,5) - (irf_abs(X(:,6:7), 1)), X(:,8)];
            % X = [irf_abs(X(:,2:4), 1), irf_dot(X(:,5:7), [1, 0, 0]), X(:,8)];
            % X = [irf_abs(X(:,5:7), 1), irf_abs(X(:,11:13), 1), X(:,15)];
            % X = [irf_abs(X(:,2:4), 1), irf_abs(X(:,8:10), 1), X(:,14)];
            [n, ~] = size(X);

            % Normalize data
            m = mean(X);
            s = std(X);
            X_norm = (X - m)./s;

            figure;
            scatter3(X(:,1), X(:,2), X(:,3), 2)
            xlabel('abs(b)');
            ylabel('abs(bulkv_i)');
            zlabel('n_i')

            y = dbscan(X_norm, 0.3, int32(n/60));

            figure;
            scatter3(X(:,1), X(:,2), X(:,3), 2, y)
            xlabel('abs(b)');
            ylabel('abs(bulkv_i)');
            zlabel('n_i')

            scal = 20;
            model1 = fitcsvm(X(y == 1, :), nonzeros(y == 1), 'Kernelscale', scal, 'Standardize', true, 'KernelFunction', 'rbf');
            model2 = fitcsvm(X(y == 2, :), nonzeros(y == 2), 'Kernelscale', scal, 'Standardize', true, 'KernelFunction', 'rbf');
            model3 = fitcsvm(X(y == -1, :), nonzeros(y == -1), 'Kernelscale', scal, 'Standardize', true, 'KernelFunction', 'rbf', 'OutlierFraction', 0.1);

            obj.SW_classifier = @classifierFunction;

            y = obj.SW_classifier(X);

            figure;
            cVec = [1 0 0; 0 1 0; 0 0 1; 1 0 1];
            % c = y == 1:max(y);
            scatter3(X(:,1), X(:,2), X(:,3), 2, cVec(y,:))
            xlabel('abs(b)');
            ylabel('abs(bulkv_i)');
            zlabel('n_i')

            % dj

            % hold on;
            % scatter3(C(:,1), C(:,2), C(:,3), 20, 'k', 'filled')
            % view(3)

            function [y, score] = classifierFunction(X, smoothWindowSize)
                [~, score1] = predict(model1, X);
                [~, score2] = predict(model2, X);
                [~, score3] = predict(model3, X);
                if nargin > 1 && smoothWindowSize > 1
                    score1 = movmean(score1, smoothWindowSize);
                    score2 = movmean(score2, smoothWindowSize);
                    score3 = movmean(score3, smoothWindowSize);
                end
                score = [score1, score2, score3];
                [~, y] = max(score, [], 2);
            end
        end
        
        function tints_pred_SLAMS = evaluate(obj, tint, tints_true_SLAMS, params)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            % tint_load = select_tint(tints_find_SLAMS, i);
            % extra_load_time = 60;
            extra_load_time = max(params.bgWindow, params.meanWindow);
            tint_load = [tint(1) + -extra_load_time, tint(2) + extra_load_time];

            % Load fgm
            fgm_filePrefix = 'mms1_fgm_srvy_l2';
            b = load_tints_MMS(fgm_filePrefix, 'mms1_fgm_b_gse_srvy_l2', tint_load);
            b_abs = irf_abs(b);

            if true
                % Load ion-fpi
                ion_fpi_filePrefix = 'mms1_fpi_fast_l2_dis-moms';
                v_ion = load_tints_MMS(ion_fpi_filePrefix, 'mms1_dis_bulkv_gse_fast', tint_load);
                v_ion_abs = irf_abs(v_ion);
                n_ion = load_tints_MMS(ion_fpi_filePrefix, 'mms1_dis_numberdensity_fast', tint_load);

                t = b_abs.time;
                % t = v_ion.time;
                % v_ion = v_ion.resample(b_abs);
                % v_ion = v_ion.resample(b_abs);
                v_ion_abs = v_ion_abs.resample(b_abs);
                n_ion = n_ion.resample(b_abs);

                X = [b_abs.data, v_ion_abs.data, n_ion.data];
                % X = [b_abs.data, v_ion.data(:, 1) - irf_abs(v_ion.data(:, 2:3), 1), n_ion.data];
                % X = [b_abs.data, irf_dot(v_ion.data(:, 1:3), [1, 0, 0]), n_ion.data];
                [y, score] = obj.SW_classifier(X);
                score = TSeries(t, score);
                % score = v_ion.data(:, 1)./(1 + irf_abs(v_ion.data(:, 2:3), 1));
                % y = double(score < -3);
                % score = TSeries(v_ion.time, score);
                % ys = TSeries(v_ion.time, y);

                ts_label = TSeries(b_abs.time, y);
                tints_cell_MSP_SW_MSH = labels2tints(ts_label, [1,2,3]);
                tints_SW = tints_cell_MSP_SW_MSH{2};
                tints_SW = remove_empty_tints(tints_SW);
            end


            % % Method1
            % [b_trimmed_data, idx] = rmoutliers(b_abs.data, 'movmean', 20*1e9, 'SamplePoints', double(b_abs.time.epoch), 'ThresholdFactor', 1.2);
            % b_trimmed = TSeries(b_abs.time(~idx), b_trimmed_data);
            % b_bg = moving_window_func(b_trimmed, 20, @(x) mean(x.data, 1), b_abs.time);
            % % b_bg = b_trimmed.resample(b_abs);
            % ts_label = TSeries(b_abs.time, double(b_abs.data >= 2*b_bg.data));

            % % Method2
            % % b_smooth = moving_window_func(b_abs, 5, @(x) mean(x.data, 1));
            % % b_floor = moving_window_func(b_smooth, 30, @(x) min(x.data, [], 1));
            % b_smooth_data = movmean(b_abs.data, 6*1e9, 'SamplePoints', b_abs.time.epoch);
            % b_floor_data = movmin(b_smooth_data, 20*1e9, 'SamplePoints', b_abs.time.epoch);
            % % c = 0.1;
            % % w = 1./(1 + c*(b_smooth_data - b_floor_data));
            % % b_bg_data = b_smooth_data.*w;
            % % b_bg = TSeries(b_abs.time, b_bg_data);
            % b_bg = TSeries(b_abs.time, b_floor_data);
            % ts_label = TSeries(b_abs.time, double(b_abs.data >= 2*b_bg.data));

            % Method3
            % b_compsmooth_data = movmean(b.data, 30*1e9, 1, 'SamplePoints', b.time.epoch);
            % b_bg = TSeries(b_abs.time, irf_abs(b_compsmooth_data, 1));
            % ts_label = TSeries(b_abs.time, double(b_abs.data >= 2*b_bg.data));

            % % Method4
            % b_smooth_data = movmedian(b_abs.data, 30*1e9, 1, 'SamplePoints', b_abs.time.epoch);
            % b_bg = TSeries(b_abs.time, b_smooth_data);
            % ts_label = TSeries(b_abs.time, double(b_abs.data >= 2*b_bg.data));

            % % Method5
            % % b_abs = TSeries(b_abs.time, movmean(b_abs.data, 2*1e9, 1, 'SamplePoints', b_abs.time.epoch));
            % % rising = diff(b_abs.data) > 0
            % b_mean_data = params.meanMult*movmean(b_abs.data, params.meanWindow*1e9, 1, 'SamplePoints', b_abs.time.epoch);
            % b_mean_harm = moving_window_func(b_abs, 30, @(x) harmmean(x.data), b_abs.time);
            % % b_mean_data = b_mean.data;
            % b_mean = TSeries(b_abs.time, b_mean_data);
            % b_bg = TSeries(b_abs.time(b_abs.data < b_mean_data), b_abs.data(b_abs.data < b_mean_data));
            % % b_bg = b_bg.resample(b_abs);
            % % b_bg  = TSeries(b_abs.time, movmean(b_bg.data, 30*1e9, 1, 'SamplePoints', b_abs.time.epoch));
            % % b_norm = TSeries(b_abs.time, b_abs.data - b_mean_data);
            % tints_peaktrough = get_tints_peaktrough(b_abs.data, b_mean_data, b_abs.time);
            % n_tints_peaktrough = length(tints_peaktrough)/2;
            % % means = zeros(1, n_tints_peaktrough);
            % % means_bg = zeros(1, n_tints_peaktrough);
            % % maxs = zeros(1, n_tints_peaktrough);
            % % for i = 1:n_tints_peaktrough
            % %     tint_peaktrough = select_tint(tints_peaktrough, i);
            % %     d = b_abs.tlim(tint_peaktrough).data;
            % %     % d_bg = b_bg.tlim([tint_peaktrough(1) + -60, tint_peaktrough(2) + 60]).data;
            % %     means(i) = mean(d);
            % %     % means_bg(i) = mean(d_bg);
            % %     % maxs(i) = max(d);
            % % end
            % SLAMS_logical = zeros(1, n_tints_peaktrough*2);
            % for i = 2:2:n_tints_peaktrough-1
            %     % m_bg = mean([means(i - 1), means(i + 1)]);
            %     tint_peaktrough = select_tint(tints_peaktrough, i);
            %     d = b_abs.tlim(tint_peaktrough).data;
            %     d_bg = b_bg.tlim([tint_peaktrough(1) + -params.bgWindow/2, tint_peaktrough(2) + params.bgWindow/2]).data;
            %     m_bg = mean(d_bg);
            %     m_p = mean(d);
            %     if m_p > params.SLAMSthresh*m_bg
            %         SLAMS_logical(i*2 - 1:i*2) = 1;
            %     end
            % end
            % tints_pred_SLAMS = tints_peaktrough(SLAMS_logical == 1);

            % ts_label = TSeries(b_abs.time, double(b_abs.data > 2*b_mean_harm.data));
            % ts_label = TSeries(b_abs.time, double(b_abs.data >= 2*b_bg.data));
            % b_bg = TSeries(tints_peaktrough, reshape([means; means], 1, [])');
            % ends_peak = find(dif < 0);
            % ends_valley = find(dif > 0);
            % starts_peak = ends_valley + 1;
            % starts_valley = ends_peak + 1;
            % if b_abs.data(1,:) > b_mean_data(1,:)
            %     starts_peak = [1, starts_peak];
            % else
            %     starts_valley = [1, starts_valley];
            % end
            % if b_abs.data(end,:) > b_mean_data(end,:)
            %     ends_peak = [ends_valley, length(b_mean_data)];
            % else
            %     ends_valley = [ends_valley, length(b_mean_data)];
            % end
            % tints_peak = b_abs.time(sort([starts_peak, ends_peak]));
            % tints_valley = b_abs.time(sort([starts_valley, ends_valley]));


            % b_mean = TSeries(b_abs.time, b_mean_data);
            % b_bg = TSeries(b_abs.time(b_abs.data < b_mean.data), b_abs.data(b_abs.data < b_mean.data));
            % b_bg_smooth = moving_window_func(b_bg, 30, @(x) mean(x.data), b_abs.time);
            % b_min_data = movmin(b_abs.data, 30*1e9, 1, 'SamplePoints', b_abs.time.epoch);
            % b_min = TSeries(b_abs.time, b_min_data);
            % b_norm = b_abs - b_min;
            % b_var_data = movvar(b_norm.data, 30*1e9, 1, 'SamplePoints', b_norm.time.epoch);
            % b_var = TSeries(b_norm.time, b_var_data);

            % Method6
            % b_bg = rmoutliers(b_abs.data, 'movmean', 30*1e9, 'SamplePoints', double(b_abs.time.epoch), 'ThresholdFactor', 1.2);
            % n_points = length(b_abs.time);
            % b_bg = [b_abs.data(1)];
            % t = b_abs.time(1);
            % for i = 1:n_points
            %     avg = moving_window_func(TSeries(t, b_bg), 30, @(x) mean(x.data), t(end));
            % end
            % b_smooth_data = b_abs.data;
            % b_smooth = TSeries(b_abs.time, b_smooth_data);
            % b_mean_data = movmean(b_smooth_data, 30*1e9, 1, 'SamplePoints', b_abs.time.epoch);
            % b_mean = TSeries(b_abs.time, b_mean_data);
            % % b_var = TSeries(b_abs.time, (b_abs.data - b_mean_data).^2);
            % b_std_data = movstd(b_abs.data, 30*1e9, 1, 'SamplePoints', b_abs.time.epoch);
            % b_std = TSeries(b_abs.time, b_std_data);
            % tints_peaktrough = get_tints_peaktrough(b_smooth_data, b_mean_data, b_abs.time);
            % b_harm = moving_window_func(b_smooth, 30, @(x) harmmean(x.data), b_abs.time);
            % ts_label = TSeries(b_abs.time, double((b_abs.data - b_mean.data) > 2*b_std.data));
            % tints_cell_NOTSLAMS_SLAMS = labels2tints(ts_label, [0,1]);
            % tints_pred_SLAMS = tints_cell_NOTSLAMS_SLAMS{2};
            % disc_array = discretize_tints(tints_peaktrough, tints_pred_SLAMS);
            % idx = 2*find(disc_array);
            % idx = sort([idx - 1, idx]);
            % tints_pred_SLAMS = tints_peaktrough(idx);
            FT = fft(b_abs.data, [], 1);
            N = fix(length(FT)/2);
            % N_clear = 3;
            % FT(1:N_clear) = 0;
            % FT(end - N_clear + 1:end) = 0;
            N_clear = N - 60;
            FT = fftshift(FT);
            FT(1:N_clear) = 0;
            FT(end - N_clear + 1:end) = 0;
            FT = ifftshift(FT);
            
            % L = b_abs.time(end) - b_abs.time(1);
            % P2 = abs(FT/L);
            % P1 = P2(1:L/2 + 1);
            % P1(2:end - 1) = 2*P1(2:end - 1);
            % f = 0.0625*(0:(L/2))/L;
            % figure;
            % plot(f, P1);
            % figure;
            % plot(b_abs.time.epoch, b_abs.data);
            % figure;
            b_fft_low = TSeries(b_abs.time, ifft(FT, 'symmetric'));

            % FT = fft(b_abs.data, [], 1);
            % N_clear = 100;  
            % FT(1:N_clear) = 0;
            % FT(end - N_clear + 1:end) = 0;
            % b_fft_high = TSeries(b_abs.time, ifft(FT, 'symmetric'));

            b_fft_1der = diff_ts(b_fft_low);
            b_fft_2der = diff_ts(b_fft_1der);

            % b_mean_m = 1.2*TSeries(b_abs.time, movmean(b_fft_low.data, 60*1e9, 1, 'SamplePoints', b_abs.time.epoch));
            b_mean = moving_window_func(b_fft_low, 60, @(x) harmmean(x.data), b_abs.time);
            % b_mean = TSeries(b_abs.time, movmedian(b_fft_low.data, 30*1e9, 1, 'SamplePoints', b_abs.time.epoch));
            % b_mean = TSeries(b_abs.time, irf_abs(movmean(b.data, 30*1e9, 1, 'SamplePoints', b_abs.time.epoch), 1));
            % b_stdev = TSeries(b_abs.time, movstd(b_abs.data, 30*1e9, 1, 'SamplePoints', b_abs.time.epoch));
            % b_dev = b_abs - b_mean;
            b_label = TSeries(b_abs.time, double(b_fft_low.data > 1.8*b_mean.data));
            % b_label = TSeries(b_abs.time, double(b_abs.data > 1.8*b_mean.data));

            % [b_label_test, b_bg] = test_onepatatime(b_abs.data, 320);
            % b_label_test = TSeries(b_abs.time, b_label_test');
            % b_bg = TSeries(b_abs.time, b_bg');

            tints_cell_NOTSLAMS_SLAMS = labels2tints(b_label, [0,1]);
            tints_pred_SLAMS = tints_cell_NOTSLAMS_SLAMS{2};

            % tints_peaktrough = get_tints_peaktrough(b_fft_low.data, b_mean_m.data, b_fft_low.time);
            % disc_array = discretize_tints(tints_peaktrough, tints_pred_SLAMS);
            % idx = 2*find(disc_array);
            % idx = sort([idx - 1, idx]);
            % tints_pred_SLAMS = tints_peaktrough(idx);

            tints_pred_SLAMS = remove_empty_tints(tints_pred_SLAMS);
            % tints_pred_SLAMS = intersect_tints(tints_pred_SLAMS, tints_SW);
            tints_pred_SLAMS = intersect_tints(tint, tints_pred_SLAMS);
            % tints_pred_SLAMS = remove_short_tints(tints_pred_SLAMS, params.minDur);
            % tints_pred_SLAMS = remove_short_tints(tints_pred_SLAMS, 2);
            if isempty(tints_pred_SLAMS)
                tints_pred_SLAMS = [];
            end

            if true
                plt = modular_plot('title', 'Predicted vs true');
                % plt.lineplot('B_smooth', {b_abs});
                plt.lineplot('B_pred', {b_abs, b_fft_low});
                plt.lineplot('B_pred', {b_mean}, 'color', 'c');
                plt.lineplot('B_pred', {1.8*b_mean}, 'color', 'r', 'ylabel', 'Pred');
                plt.lineplot('B_pred_state', b_label_test)
                plt.lineplot('B_pred_bg', b_bg)
                % plt.lineplot('B_pred', {b_mean_m}, 'color', 'g');
                % plt.lineplot('B_pred', {b_bg}, 'color', 'c');
                % plt.lineplot('B_pred', {2*b_harm}, 'color', 'r', 'ylabel', 'Pred');
                plt.lineplot('B_true_cont', {b, b_abs}, 'ylabel', 'True');
                plt.lineplot('n', n_ion, 'ylabel', 'n_i (cm^{-3})')
                % plt.lineplot('T_ion', {T_ion_para, T_ion_perp}, 'ylabel', 'T_{i} (eV)', 'legend', {'T_{ipara}','T_{iperp}'})
                % plt.spectrogram('E_ion_spectr', {E_ion_spectr, E_ion}, 'ylabel', 'W_i (eV)')
                plt.lineplot('v_ion', {v_ion, v_ion_abs}, 'ylabel', 'v_{i} (kms^{-1})', 'legend', {'v_{ix}','v_{iy}','v_{iz}'})
                plt.lineplot('class_score', score, 'ylabel', 'Score', 'legend', {'Msphere', 'SW','Msheath'}, 'colorOrder', [1 0 0; 0 0.5 0; 0 0 1])
                % plt.lineplot('pos', pos, 'ylabel', 'R_{GSE} (R_E)', 'legend', {'x','y','z'})
                plt.show(tint);
                plt.mark(tints_pred_SLAMS, 'target', 'B_pred', 'intervals', true);
                % plt.mark(tints_pred_SLAMS, 'target', 'B_tmp', 'intervals', true);
                plt.mark(tints_true_SLAMS, 'target', 'B_true_cont', 'intervals', true);
                plt.mark(tints_cell_MSP_SW_MSH, 'target', 'class_score', 'intervals', true);
            end
        end
    end
end


