classdef SLAMS_finder < handle
    %SLAMS_FINDER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        region_classifier
        ts
        last_run_results
        model_region
        tint_load
    end
    
    methods
        function obj = SLAMS_finder()
            %SLAMS_FINDER Construct an instance of this class
            %   Detailed explanation goes here
            obj.ts = struct;
            obj.last_run_results = struct;
            obj.SW_train(true);
        end

        function data = create_region_data(~, E_ion, E_ion_omni, pos)
            E_bins = log(E_ion);
            % E_bins = 1:32;
            tot = sum(E_ion_omni, 2);
            E_mean = sum(E_ion_omni.*E_bins, 2)./tot;
            E_MAD = sum(E_ion_omni.*(abs(E_bins - E_mean)), 2)./tot;
            cosang = pos(:,1)./(sqrt(sum(pos.^2, 2)));
            cosang(isnan(cosang)) = 1;
            % cosang = (cosang/2 + 1/2).^5;
            % dis(isnan(dis)) = 1;
            % v_rel = log(abs(v_ion(:,1))./sqrt(sum(v_ion(:,2:3).^2, 2)));
            data = [E_MAD, E_mean, cosang];
        end

        % function y = cluster_tree(tree, n)
            
        % end

        function SW_train(obj, plot_progress)
            disp('Training SW classifier...')
            % X1 = readmatrix(['cache/', 'traindata1.txt']);
            X4 = readmatrix(['cache/', 'traindata4.txt']);
            % X5 = readmatrix(['cache/', 'traindata5.txt']);
            % train_data2 = readmatrix([cache_path, 'traindata2.txt']);

            % X = [irf_abs(X1(:,2:4), 1), irf_abs(X1(:,5:7), 1), X1(:,8)];
            % X_n = X1(:,8);
            % X_b = irf_abs(X1(:,2:4), 1);
            E_ion = X4(:,37:68);
            % E_ion = X5(:,2:33);
            E_ion_omni = X4(:,5:36);
            % v_ion = X4(:,2:4);
            pos = X4(:,69:71);
            X = obj.create_region_data(E_ion, E_ion_omni, pos);
            % X = [X_e_m, X1(:,1), X_e_v];
            % X = [irf_abs(X(:,2:4), 1), X(:,5) - (irf_abs(X(:,6:7), 1)), X(:,8)];
            % X = [irf_abs(X(:,2:4), 1), irf_dot(X(:,5:7), [1, 0, 0]), X(:,8)];
            % X = [irf_abs(X(:,5:7), 1), irf_abs(X(:,11:13), 1), X(:,15)];
            % X = [irf_abs(X(:,2:4), 1), irf_abs(X(:,8:10), 1), X(:,14)];
            [n, ~] = size(X);

            % % Normalize data
            % m = mean(X);
            % s = std(X);
            % X_norm = (X - m)./s;

            if plot_progress
                figure;
                scatter3(X(:,1), X(:,2), X(:,3), 2)
                title('Raw data')
                xlabel('log(E_{ion} bins) MAD');
                ylabel('log(E_{ion} bins) Mean');
                zlabel('cos(angle)')
            end

            % y = dbscan(X_norm, 0.3, int32(n/60));
            % y = dbscan(X_norm, 0.2, int32(n/180));

            % Z = linkage(X_norm, 'single', 'euclidean');
            % % size(Z)
            % figure;
            % dendrogram(Z,1000)
            % df

            % tes = pdist(X_norm);
            % sq = squareform(tes);
            % % sq(sq == 0) = 10000;
            % mindist = mean(sq, 2);
            % seq = 1:10000;
            % figure;
            % plot(seq, mindist)
            % [rmoutlierdist, idx] = rmoutliers(mindist, 'percentiles', [0, 80]);
            % figure;
            % plot(seq(idx == 0), mindist(idx == 0))
            % X_norm = X_norm(idx == 0, :);
            % X = X(idx == 0, :);

            % try
            %     % df
            %     load('model.mat', 'model')
            % catch
            %     rng(90)
            %     model = fitgmdist(X, 3, 'CovarianceType', 'full', 'SharedCovariance', false, 'Replicates', 1, 'Start', 'randSample', 'Options', statset('Display', 'final', 'MaxIter', 1000,  'TolFun', 1e-12));
            %     save('model.mat', 'model')
            % end
            % y = cluster(model, X);

            % [y, ~, ~, ~, mahaldist] = cluster(model, X);
            % outliers = all(sqrt(mahaldist) > 2, 2)
            % unique(outliers)
            % X_reduced = X(~outliers, :);
            % rng(90)
            % model = fitgmdist(X_reduced, 3, 'CovarianceType', 'full', 'SharedCovariance', false, 'Replicates', 1, 'Start', 'randSample', 'Options', statset('Display', 'final', 'MaxIter', 1000,  'TolFun', 1e-12));
            % y_reduced = cluster(model, X_reduced);
            % y = cluster(model, X);

            % Hierarchical clustering.
            tree = linkage(X, 'average', 'mahalanobis');
            % % y = cluster(Z, 'Cutoff', 0.1, 'Criterion', 'distance');

            % Lowest number of clusters that separates MSH, SW and MSP.
            n_clust = 10;
            y = cluster(tree, 'MaxClust', n_clust);
            
            if plot_progress
                figure;
                scatter3(X(:,1), X(:,2), X(:,3), 2, y)
                title('Hierarchical clusters')
                xlabel('log(E_{ion} bins) MAD');
                ylabel('log(E_{ion} bins) Mean');
                zlabel('cos(alpha)')
            end

            % Get mu, sigma and p of each cluster.
            mu = zeros(n_clust, 3);
            sig = zeros(3, 3, n_clust);
            p = zeros(n_clust, 1);

            for i = 1:n_clust
                logi = y == i;
                X_clust = X(logi, :);
                mu(i, :) = mean(X_clust, 1);
                sig(:, :, i) = cov(X_clust);
                p(i) = nnz(logi)/n;
            end

            % If cluster is small (p is low) remove it.
            idx_small = p < 0.03;
            mu(idx_small, :) = [];
            sig(:, :, idx_small) = [];
            p(idx_small) = [];

            % Create GMM of remaining clusters.
            model = gmdistribution(mu, sig, p);
            y = cluster(model, X);

            if plot_progress
                figure;
                cVec = [1 0 0; 0 0.8 0; 0 0 1; 0 1 1; 1 0 1; 0 0 0];
                scatter3(X(:,1), X(:,2), X(:,3), 2, cVec(y, :))
                title('GMM clusters')
                xlabel('log(E_{ion} bins) MAD');
                ylabel('log(E_{ion} bins) Mean');
                zlabel('cos(alpha)')
            end

            % df

            % merge_map = {[1, 3], [5, 6], [2, 4]};
            % y = merge_clusters(y, merge_map);
            % merge_map = {[2, 5], 1, [3, 4, 6, 7, 8, 9, 10]};
            % merge_map = {[2, 5], 1, [3, 4, 6, 7, 8, 9, 10, 11]};

            % Create a map to merge gaussians which belong to the same class.
            merge_map = {[2, 3], 1, [4, 5, 6]};
            y = merge_clusters(y, merge_map);

            % ys = unique(y);
            % counts = hist(y, ys);
            % [~, idx] = sort(counts, 'descend');
            % rang = 1:2;
            % ys = ys(idx(rang));
            % idx = ismember(y, ys);
            % X_trim = X(idx, :);
            % y_trim = y(idx);

            % y(y == 1200) = 1;
            % y(y == 1201) = 2;
            % y(y == 1203) = 3;
            % y(y == 1205) = 4;
            % y(idx == 0) = -1;
            % y(idx) = rng;
            % unique(y_trim)
            % y = cluster(Z, 'Maxclust', 4);

            % figure;
            % histogram(y, length(unique(y)))
            % y = spectralcluster(X_norm, 5)
            % unique(y)

            if plot_progress % && false
                figure;
                cVec = [1 0 0; 0 0.8 0; 0 0 1; 0 1 1; 1 0 1];
                scatter3(X(:,1), X(:,2), X(:,3), 2, cVec(y, :))
                title('Merged clusters')
                xlabel('log(E_{ion} bins) MAD');
                ylabel('log(E_{ion} bins) Mean');
                zlabel('cos(alpha)')
            end

            % df

            % model1 = fitcsvm(X(y == 1, :), nonzeros(y == 1), 'Kernelscale', scal, 'Standardize', true, 'KernelFunction', 'rbf');
            % model2 = fitcsvm(X(y == 2, :), nonzeros(y == 2), 'Kernelscale', scal, 'Standardize', true, 'KernelFunction', 'rbf');
            % model3 = fitcsvm(X(y == -1, :), nonzeros(y == -1), 'Kernelscale', scal, 'Standardize', true, 'KernelFunction', 'rbf', 'OutlierFraction', 0.1);
            % sw_idx = y_trim == 1200;
            % numClasses = length(unique(y));
            % model = cell(numClasses, 1);
            % % rng(100);
            % for i = 1:numClasses
            %     y_tmp = y == i;
            %     unique(y_tmp)
            %     % if i == 3
            %     %     model{i} = fitcsvm(X, y_tmp, 'BoxConstraint', 0.03, 'Kernelscale', 1, 'Standardize', true, 'KernelFunction', 'rbf', 'OutlierFraction', 0, 'Weights', y_tmp + 0.01);
            %     % else BoxConstraint = 0.03
            %         model{i} = fitcsvm(X(y_tmp, :), y_tmp(y_tmp), 'BoxConstraint', 1, 'Kernelscale', 1, 'Standardize', true, 'KernelFunction', 'rbf', 'OutlierFraction', 0.10, 'ScoreTransform', @(x) (0.01*x./(abs(0.01*x) + 1) + 1)/2);
            %     % end
            %     % model{i} = fitPosterior(model{i});

            %     % X_tmp = X(y_tmp, :);
            %     % model{i} = struct('mu', mean(X_tmp, 1), 'sigma', cov(X_tmp));
            % end

            obj.region_classifier = @classifierFunction;
            obj.model_region = model;

            y = obj.region_classifier(X);

            % [t, probs, mahaldist, logpdf, nlogL] = obj.region_classifier([0.5, 7, 0])
            % unique(y)

            if plot_progress
                % Create boundary points
                n_bound = 100;
                % [b_x, b_y, b_z] = meshgrid(linspace(4, 10, n_bound), linspace(-400, 800, n_bound), linspace(-1, 1, n_bound));
                [b_x, b_y, b_z] = meshgrid(linspace(0, 3, n_bound), linspace(4, 10, n_bound), linspace(-1, 1, n_bound));
                X_bound = [reshape(b_x, [], 1, 1), reshape(b_y, [], 1, 1), reshape(b_z, [], 1, 1)];
                [~, prob_test, ~, logpdf, ~] = obj.region_classifier(X_bound);
                % % prob_test(1:10,:)
                % dif_probs = abs([ ...
                %     prob_test(:,2) - prob_test(:,3), ...
                %     prob_test(:,1) - prob_test(:,3), ...
                %     prob_test(:,1) - prob_test(:,2) ...
                % ]);
                % small_dif = (dif_probs < 0.05);
                % [~, idx] = min(prob_test, [], 2, 'linear');
                % bound = small_dif(idx);
                % bound = logpdf > -5 & logpdf < -4;
                bound = any(prob_test > 0.45 & prob_test < 0.55, 2);
                % bound = any(prob_test > 0.01 & prob_test < 0.05, 2);
                % bound = any(prob_test > 0.85 & prob_test < 0.95, 2);
                % bound = any(prob_test > 0.15 & prob_test < 0.25, 2);
                % bound = any(prob_test > -0.05 & prob_test < 0.05, 2);
                % bound = all(prob_test > 0.00 & prob_test < 0.05, 2);
                X_bound = X_bound(bound,:);

                figure;
                cVec = [1 0 0; 0 0.8 0; 0 0 1; 0 1 1; 1 0 1];
                % c = y == 1:max(y);
                scatter3(X(:,1), X(:,2), X(:,3), 2, cVec(y, :))
                hold on
                scatter3(X_bound(:,1), X_bound(:,2), X_bound(:,3), 2, [0, 0, 0])
                title('Classified raw data')
                xlabel('log(Energy) MAD');
                ylabel('log(Energy) Mean');
                zlabel('cos(alpha)')
            end

            % ff

            % hold on;
            % scatter3(C(:,1), C(:,2), C(:,3), 20, 'k', 'filled')
            % view(3)

            function [y, probs, mahaldist, logpdf, nlogL] = classifierFunction(X)
                [n, ~] = size(X);
                l = length(merge_map);
                probs = zeros(n, l);
                mahaldist = zeros(n, l);
                [y, nlogL, P, logpdf, d2] = cluster(model, X);
                % post = posterior(model, X);
                % [y,nlogL,P,logpdf,d2] = cluster(model, X);
                % probs = pdf(model, X);
                % probs = logpdf;
                % prop = zeros(n, l);
                % mahaldists = sqrt(mahal(model, X));
                % % idx = zeros(n, l)
                % % model.ComponentProportion
                % % unique(prop)
                % mahaldists = reshape(mahaldists, [], 1);
                % tmp_probs = 2*normcdf(-mahaldists, 0, 1);
                % tmp_probs = reshape(tmp_probs, [], model.NumComponents);
                % tmp_probs = tmp_probs.*(model.ComponentProportion);
                % model.ComponentProportion
                for j = 1:l
                    probs(:,j) = sum(P(:, merge_map{j}), 2);
                    mahaldist(:,j) = min(d2(:, merge_map{j}), [], 2);
                    % idx(:,i) = model.ComponentProportion(merge_map{i}(idx));
                end
                mahaldist = sqrt(mahaldist);
                % for j = 1:numClasses
                %     [~, prob] = predict(model{j}, X);
                %     % prob = sum(post(:, correction{j}), 2);
                %     % prob = mvnpdf(X, model{j}.mu, model{j}.sigma);
                %     % size(prob)
                %     % prob = mvksdensity(model{j}, X, 'Bandwidth', 0.01, 'Kernel', 'epanechnikov');
                %     if nargin > 1 && smoothWindowSize > 1
                %         prob = movmean(prob, smoothWindowSize*1e9, 1, 'SamplePoints', t.epoch);
                %         % score2 = movmean(score2, smoothWindowSize*1e9, 1, 'SamplePoints', t.epoch);
                %         % score3 = movmean(score3, smoothWindowSize*1e9, 1, 'SamplePoints', t.epoch);
                %     end
                %     % prob = (score + model{j}.Nu).*((model{j}.KernelParameters.Scale/pi)^(3/2)/sum(model{j}.Alpha));
                %     % prob = score;
                %     probs(:, j) = prob;
                %     % score = [score1, score2, score3];
                %     % score = score1;
                % end
                % % probs(1:10,:)
                % [~, y] = min(mahaldist, [], 2);
                % % y = cluster(model, X);
                % % probs = posterior(model, X);
                % % for i = 1:l
                % %     probs(:,i) = max(probs(:, merge_map{i}), [], 2);
                % %     % idx(:,i) = model.ComponentProportion(merge_map{i}(idx));
                % % end
                y = merge_clusters(y, merge_map);
            end

            function y = merge_clusters(y, map)
                l = length(map);
                idx = cell(1, l);
                for k = 1:l
                    idx{k} = ismember(y, map{k});
                end
                for k = 1:l
                    y(idx{k}) = k;
                end
            end
        end

        function region_classify(obj, r)
            obj.load_ts({'E_ion', 'E_ion_omni', 'pos'});

            t = obj.ts.E_ion_omni.time;
            % % t = v_ion.time;
            % % v_ion = v_ion.resample(b_abs);
            % % v_ion = v_ion.resample(b_abs);
            % v_ion_abs = obj.ts.v_ion_abs.resample(obj.ts.b_abs);
            % n_ion = obj.ts.n_ion.resample(obj.ts.b_abs);
            pos = obj.ts.pos.resample(obj.ts.E_ion_omni);

            % X = [obj.ts.b_abs.data, v_ion_abs.data, n_ion.data];

            X = obj.create_region_data(obj.ts.E_ion.data, obj.ts.E_ion_omni.data, pos.data);
            % size(X)

            % X = [b_abs.data, v_ion.data(:, 1) - irf_abs(v_ion.data(:, 2:3), 1), n_ion.data];
            % X = [b_abs.data, irf_dot(v_ion.data(:, 1:3), [1, 0, 0]), n_ion.data];
            % [y, prob] = obj.region_classifier(X, r.SW_smoothing, t);
            % n = length(r.smoothWindowSize) + 1;
            [y, probs, mahaldist, logpdf, nlogL] = obj.region_classifier(X);
            % if r.smoothWindowSize > 0
            %     probs = movmean(probs, smoothWindowSize*1e9, 1, 'SamplePoints', t.epoch);
            %     logpdf = movmean(logpdf, smoothWindowSize*1e9, 1, 'SamplePoints', t.epoch);
            %     mahaldist = movmean(mahaldist, smoothWindowSize*1e9, 1, 'SamplePoints', t.epoch);
            % end
            probs = TSeries(t, probs);
            logpdf = TSeries(t, logpdf);
            mahaldist = TSeries(t, mahaldist);
            obj.last_run_results.probs = probs;
            obj.last_run_results.logpdf = logpdf;
            obj.last_run_results.mahaldist = mahaldist;
            % region = struct('nlogL', nlogL, )
            % prob = v_ion.data(:, 1)./(1 + irf_abs(v_ion.data(:, 2:3), 1));
            % y = double(prob < -3);
            % prob = TSeries(v_ion.time, prob);
            % ys = TSeries(v_ion.time, y);

            ts_label = TSeries(t, y);
            % ts_label = moving_window_func(ts_label, 60, @(x) majorityVote(x.data), t);
            tints_cell_region = labels2tints(ts_label, 1:obj.model_region.NumComponents);
            obj.last_run_results.tints_cell_region = tints_cell_region;
            % tints_SW = tints_cell_region{2};
            % % tints_SW = remove_short_tints(tints_SW, 60);
            % tints_SW = remove_empty_tints(tints_SW);

            obj.last_run_results.SW_plotter = @plotter;
            function plotter(plt)
                plt.lineplot('class_probs', mahaldist, 'ylabel', 'Posterior probability', 'legend', {'MSP', 'SW', 'MSH'}, 'colorOrder', [1 0 0; 0 0.8 0; 0 0 1])
            end

            function v = majorityVote(x)
                [~, v] = max(accumarray(x,1));
            end
        end
        
        function SLAMS = SLAMS_classify(obj, r)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            obj.load_ts({'b', 'b_abs'})

            t = obj.ts.b_abs.time;

            switch r.SLAMS_B_bg_method
            case 'mean'
                obj.ts.b_bg = TSeries(t, movmean(obj.ts.b_abs.data, r.SLAMS_B_bg_window*1e9, 1, 'SamplePoints', obj.ts.b_abs.time.epoch));
            case 'median'
                obj.ts.b_bg = TSeries(t, movmedian(obj.ts.b_abs.data, r.SLAMS_B_bg_window*1e9, 1, 'SamplePoints', obj.ts.b_abs.time.epoch));
            case 'harmmean'
                obj.ts.b_bg = moving_window_func(obj.ts.b_abs, r.SLAMS_B_bg_window, @(x) harmmean(x.data), t);
            end

            b_upper_thresh = r.SLAMS_Threshold*obj.ts.b_bg;
            b_lower_thresh = ((r.SLAMS_Threshold - 1)/2 + 1)*obj.ts.b_bg;

            logic_above = obj.ts.b_abs.data > b_upper_thresh.data;
            logic_under = obj.ts.b_abs.data < b_lower_thresh.data;
            logic_SLAMS = merge_SLAMS(logic_above, logic_under);

            SLAMS_label = TSeries(t, double(logic_SLAMS));
            
            tints_cell_NOTSLAMS_SLAMS = labels2tints(SLAMS_label, [0,1]);
            tints_SLAMS = tints_cell_NOTSLAMS_SLAMS{2};

            SLAMS = tints_SLAMS;
            
            % tints_SLAMS = merge_tints(tints_SLAMS, r.SLAMS_DeltaMaxMerge);
            % tints_SLAMS = remove_short_tints(tints_SLAMS, r.SLAMS_MinDur);
            % tints_SLAMS = remove_empty_tints(tints_SLAMS);

            % SLAMS = struct('startUTC', tint_SLAMS(1:2:end).utc, 'endUTC', tint_SLAMS(2:2:end).utc)
            
            obj.last_run_results.SLAMS_plotter = @plotter;
            function plotter(plt, name, arg)
                plt.lineplot(name, obj.ts.b_bg, 'color', 'c');
                % plt.lineplot(name, b_trigger, 'color', 'g');
                plt.lineplot(name, b_lower_thresh, 'color', 'm');
                plt.lineplot(name, b_upper_thresh, 'color', 'r', arg{:});
            end

            function logic_above = merge_SLAMS(logic_above, logic_under)
                if all(logic_above == 0)
                    return
                end
                logic_above = logic_above';
                logic_under = logic_under';
            
                index_above = find(logic_above);
                index_under = find(logic_under);
                l = length(index_under);
                tmp = index_under' > index_above;
                tmp = diff(tmp, 1);
                tmp = mod(find(tmp), l - 1);
                tmp(tmp == 0) = l - 1;
                tmp = unique(tmp);
                tmp = sort([tmp; tmp + 1]);
                tmp = [1; tmp; l];
                tmp = index_under(tmp);
                for i = 1:2:length(tmp)
                    logic_under(tmp(i):tmp(i+1)) = 1;
                end
                logic_above = ~logic_under';
            end
        end

        function tints_SLAMS = SLAMS_classify_experimental(obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            obj.load_ts({'b', 'b_abs'})

            FT = fft(obj.ts.b_abs.data, [], 1);
            N = fix(length(FT)/2);
            N_clear = N - 60;
            FT = fftshift(FT);
            FT(1:N_clear) = 0;
            FT(end - N_clear + 1:end) = 0;
            FT = ifftshift(FT);
            
            b_smooth = TSeries(obj.ts.b_abs.time, ifft(FT, 'symmetric'));

            [B, TF] = rmoutliers(b_smooth.data, 'movmean', 60*1e9, 'ThresholdFactor', 1.5, 'SamplePoints', double(obj.ts.b_abs.time.epoch));
            b_o = TSeries(obj.ts.b_abs.time(~TF), B);
            b_mean = TSeries(obj.ts.b_abs.time, movmean(obj.ts.b_abs.data, 60*1e9, 1, 'SamplePoints', obj.ts.b_abs.time.epoch));
            % b_mean = moving_window_func(obj.ts.b_abs, 30, @(x) harmmean(x.data), obj.ts.b_abs.time);
            % b_std = TSeries(obj.ts.b_abs.time, movstd(obj.ts.b_abs.data, 30*1e9, 1, 'SamplePoints', obj.ts.b_abs.time.epoch));
            % b_dev = TSeries(obj.ts.b_abs.time, obj.ts.b_abs.data - b_mean.data);
            % b_reldev = b_dev/b_mean;
            % [label, bg] = test_onepatatime(obj.ts.b_abs.data, 320);
            % b_label = TSeries(obj.ts.b_abs.time, label');
            % b_bg = TSeries(obj.ts.b_abs.time, bg');

            b_label = TSeries(obj.ts.b_abs.time, double(TF));

            tints_cell_NOTSLAMS_SLAMS = labels2tints(b_label, [0,1]);
            tints_SLAMS = tints_cell_NOTSLAMS_SLAMS{2};

            obj.last_run_results.SLAMS_plotter = @plotter;
            function plotter(plt, name, arg)
                plt.lineplot(name, b_mean, 'color', 'c');
                plt.lineplot(name, b_smooth, 'color', 'g', arg{:});
            end
        end

        function plot_comparison(obj, tint, tints_true_SLAMS)
            obj.load_ts({'b', 'b_abs', 'v_ion', 'v_ion_abs', 'n_ion', 'E_ion_omni', 'E_ion', 'pos'})

            plt = modular_plot('title', 'Predicted vs true');
            if isfield(obj.last_run_results, 'SW_plotter')
                obj.last_run_results.SW_plotter(plt);
            end
            plt.lineplot('B_pred', obj.ts.b_abs);
            if isfield(obj.last_run_results, 'SLAMS_plotter')
                obj.last_run_results.SLAMS_plotter(plt, 'B_pred', {'ylabel', 'Pred'})
            end
            plt.lineplot('B_true', obj.ts.b_abs, 'ylabel', 'True');
            plt.lineplot('B', {obj.ts.b, obj.ts.b_abs}, 'ylabel', 'B_{GSE} (nT)');
            plt.lineplot('n_ion', obj.ts.n_ion, 'ylabel', 'n_i (cm^{-3})')
            % plt.lineplot('T_ion', {T_ion_para, T_ion_perp}, 'ylabel', 'T_{i} (eV)', 'legend', {'T_{ipara}','T_{iperp}'})
            plt.spectrogram('E_ion_omni', {obj.ts.E_ion_omni, obj.ts.E_ion}, 'ylabel', 'W_i (eV)')
            plt.lineplot('v_ion', {obj.ts.v_ion, obj.ts.v_ion_abs}, 'ylabel', 'v_{i} (kms^{-1})', 'legend', {'v_{ix}','v_{iy}','v_{iz}'})
            plt.lineplot('pos', obj.ts.pos, 'ylabel', 'R_{GSE} (R_E)', 'legend', {'x','y','z'})
            plt.show(tint);
            plt.mark(obj.last_run_results.tints_SLAMS, 'target', 'B_pred', 'intervals', true);
            plt.mark(tints_true_SLAMS, 'target', 'B_true', 'intervals', true);
            if isfield(obj.last_run_results, 'SW_plotter')
                plt.mark(obj.last_run_results.tints_cell_region, 'target', 'class_probs', 'intervals', true);
            end
        end

        function plot_prediction(obj, tint)
            obj.load_ts({'b', 'b_abs', 'v_ion', 'v_ion_abs', 'n_ion', 'E_ion_omni', 'E_ion', 'pos'})

            plt = modular_plot('title', 'Predicted SLAMS');
            if isfield(obj.last_run_results, 'SW_plotter')
                obj.last_run_results.SW_plotter(plt);
            end
            plt.lineplot('B_pred', obj.ts.b_abs);
            if isfield(obj.last_run_results, 'SLAMS_plotter')
                obj.last_run_results.SLAMS_plotter(plt, 'B_pred', {'ylabel', 'Pred'})
            end
            plt.lineplot('B', {obj.ts.b, obj.ts.b_abs}, 'ylabel', 'B_{GSE} (nT)');
            plt.lineplot('n_ion', obj.ts.n_ion, 'ylabel', 'n_i (cm^{-3})')
            plt.lineplot('E_s', obj.ts.E_ion_omni, 'ylabel', 'spectr')
            plt.lineplot('E', obj.ts.E_ion, 'ylabel', 'energy')
            % plt.lineplot('T_ion', {T_ion_para, T_ion_perp}, 'ylabel', 'T_{i} (eV)', 'legend', {'T_{ipara}','T_{iperp}'})
            plt.spectrogram('E_ion_omni', {obj.ts.E_ion_omni, obj.ts.E_ion}, 'ylabel', 'W_i (eV)')

            E_bins = obj.ts.E_ion.data;
            % E_bins = (1:32)/32;
            tot = sum(obj.ts.E_ion_omni.data, 2);
            X_e_m = sum(obj.ts.E_ion_omni.data.*E_bins, 2)./tot;
            X_e_v = sum(obj.ts.E_ion_omni.data.*((E_bins - X_e_m).^6), 2)./tot;
            test1 = TSeries(obj.ts.E_ion_omni.time, X_e_m);
            test2 = TSeries(obj.ts.E_ion_omni.time, X_e_v);
            plt.lineplot('test1', test1, 'ylabel', 'test1')
            plt.lineplot('test2', test2, 'ylabel', 'test2')

            plt.lineplot('v_ion', {obj.ts.v_ion, obj.ts.v_ion_abs}, 'ylabel', 'v_{i} (kms^{-1})', 'legend', {'v_{ix}','v_{iy}','v_{iz}'})
            plt.lineplot('pos', obj.ts.pos, 'ylabel', 'R_{GSE} (R_E)', 'legend', {'x','y','z'})
            plt.show(tint);
            plt.mark(obj.last_run_results.tints_SLAMS, 'target', 'B_pred', 'intervals', true);
            if isfield(obj.last_run_results, 'SW_plotter')
                plt.mark(obj.last_run_results.tints_cell_region, 'target', 'class_probs', 'intervals', true);
            end
        end

        function [tints_SLAMS, B_max, B_max_rel] = evaluate(obj, tint, varargin)

            p = inputParser;
            addRequired(p, 'tint');
            addParameter(p, 'Experimental', false)
            addParameter(p, 'Use_SW_classifier', true)
            addParameter(p, 'SW_smoothing', 0)
            addParameter(p, 'ExtraLoadTime', 60)
            addParameter(p, 'SLAMS_B_bg_method', 'median')
            addParameter(p, 'SLAMS_B_bg_window', 60)
            addParameter(p, 'SLAMS_Threshold', 2)
            addParameter(p, 'SLAMS_DeltaMaxMerge', 0.5)
            addParameter(p, 'SLAMS_MinDur', 1)
            addParameter(p, 'Compute_B_max', false)
            parse(p, tint, varargin{:})
            r = p.Results;

            obj.ts = rmfield(obj.ts, fieldnames(obj.ts));
            obj.last_run_results = rmfield(obj.last_run_results, fieldnames(obj.last_run_results));

            obj.tint_load = [r.tint(1) + -r.ExtraLoadTime, r.tint(2) + r.ExtraLoadTime];
            if r.Experimental
                tints_SLAMS = obj.SLAMS_classify_experimental();
            else
                tints_SLAMS = obj.SLAMS_classify(r);
            end
            if r.Use_SW_classifier && ~r.Experimental
                obj.region_classify(r);
                % tints_SLAMS_SW = intersect_tints(tints_SW, tints_SLAMS);
                % tints_SLAMS = intersect_tints(tints_SW, tints_SLAMS);
                % tints_SLAMS = remove_edge_SLAMS(tints_SLAMS_SW, tints_SLAMS);
            end
            tints_SLAMS = intersect_tints(tints_SLAMS, r.tint);
            tints_SLAMS = remove_empty_tints(tints_SLAMS);

            obj.last_run_results.tints_SLAMS = tints_SLAMS;

            if r.Compute_B_max
                [B_max, B_max_rel] = get_B_max_SLAMS(tints_SLAMS, obj.ts.b_abs, obj.ts.b_bg);
            else
                B_max = [];
                B_max_rel = [];
            end

            function tints_SLAMS = remove_edge_SLAMS(tints_SLAMS_SW, tints_SLAMS)
                if isempty(tints_SLAMS_SW)
                    tints_SLAMS = [];
                    return
                end
                mem = ismember(tints_SLAMS_SW.epoch, tints_SLAMS.epoch);
                dif = diff(mem);
                idx = dif(1:2:end) == 0;
                idx = reshape([idx, idx]', [], 1);
                tints_SLAMS = tints_SLAMS_SW(idx);
            end

            function [B_max, B_max_rel] = get_B_max_SLAMS(tints_SLAMS, b_abs, b_bg)
                %GET_B_MAX_SLAMS Summary of this function goes here
                %   Detailed explanation goes here
                if isempty(tints_SLAMS)
                    B_max = [];
                    B_max_rel = [];
                    return
                end
            
                n_SLAMS = length(tints_SLAMS)/2;
                B_max = zeros(n_SLAMS, 1);
                B_max_rel = zeros(n_SLAMS, 1);
                for i = 1:n_SLAMS
                    tint_SLAMS = select_tint(tints_SLAMS, i);
                    b_data_SLAMS = b_abs.tlim(tint_SLAMS).data;
                    b_rel_data_SLAMS = b_data_SLAMS./b_bg.tlim(tint_SLAMS).data;
                    B_max(i) = max(b_data_SLAMS, [], 1);
                    B_max_rel(i) = max(b_rel_data_SLAMS, [], 1);
                end
            end
        end

        function load_ts(obj, cell_string)
            for i = 1:length(cell_string)
                switch cell_string{i}
                case 'b'
                    if ~isfield(obj.ts, 'b')
                        obj.ts.b = load_tints_MMS('mms1_fgm_srvy_l2', 'mms1_fgm_b_gse_srvy_l2', obj.tint_load);
                    end
                case 'b_abs'
                    if ~isfield(obj.ts, 'b_abs')
                        obj.ts.b_abs = irf_abs(obj.ts.b);
                    end
                case 'v_ion'
                    if ~isfield(obj.ts, 'v_ion')
                        obj.ts.v_ion = load_tints_MMS('mms1_fpi_fast_l2_dis-moms', 'mms1_dis_bulkv_gse_fast', obj.tint_load);
                    end
                case 'v_ion_abs'
                    if ~isfield(obj.ts, 'v_ion_abs')
                        obj.ts.v_ion_abs = irf_abs(obj.ts.v_ion);
                    end
                case 'n_ion'
                    if ~isfield(obj.ts, 'n_ion')
                        obj.ts.n_ion = load_tints_MMS('mms1_fpi_fast_l2_dis-moms', 'mms1_dis_numberdensity_fast', obj.tint_load);
                    end
                case 'E_ion_omni'
                    if ~isfield(obj.ts, 'E_ion_omni')
                        obj.ts.E_ion_omni = load_tints_MMS('mms1_fpi_fast_l2_dis-moms', 'mms1_dis_energyspectr_omni_fast', obj.tint_load);
                    end
                case 'E_ion'
                    if ~isfield(obj.ts, 'E_ion')
                        obj.ts.E_ion = load_tints_MMS('mms1_fpi_fast_l2_dis-moms', 'mms1_dis_energy_fast', obj.tint_load);
                    end
                case 'pos'
                    if ~isfield(obj.ts, 'pos')
                        obj.ts.pos = load_tints_MMS('mms1_mec_srvy_l2_ephts04d', 'mms1_mec_r_gse', obj.tint_load)/6378;
                    end
                otherwise
                    error(['''' cell_string{i} ''' loader not implemented.'])
                end
            end
        end
    end
end

