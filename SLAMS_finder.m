classdef SLAMS_finder < handle
    %SLAMS_FINDER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        region_classifier
        ts
        last_run_results
        model_region
        tint_load
        tint
    end
    
    methods
        function obj = SLAMS_finder(varargin)
            %SLAMS_FINDER Construct an instance of this class
            %   Detailed explanation goes here

            p = inputParser;
            addParameter(p, 'Show_train_progress', false)
            parse(p, varargin{:})
            r = p.Results;

            obj.ts = struct;
            obj.last_run_results = struct;
            obj.region_train(r.Show_train_progress);
        end

        function data = create_region_data(~, E_ion_bin_center, E_ion_bin_delta, E_ion_omni, pos, v_ion)
            E_bins_center = log(E_ion_bin_center);
            % E_bins_center = 1:32;
            E_per_bin = E_ion_omni.*E_ion_bin_delta;
            % N_per_bin = E_per_bin./E_bins_center;
            % E_per_bin = E_ion_omni;
            % N_bins = E_per_bin./E_bins_center;
            E_tot = sum(E_per_bin, 2);
            % N_tot = sum(N_per_bin, 2);
            % Echar = E_tot./N_tot;
            E_mean = sum(E_per_bin.*E_bins_center, 2)./E_tot;
            E_MAD = sum(E_per_bin.*(abs(E_bins_center - E_mean)), 2)./E_tot;
            % cosang = pos(:,1)./(sqrt(sum(pos.^2, 2)));
            % cosang(isnan(cosang)) = 1;
            ang = acos(-v_ion(:,1)./sqrt(sum(v_ion.^2, 2)))/pi;
            ang = log(ang);
            % mean(v_ion, 1)
            % size(v_ion)
            % tmp = v_ion(:,1).*(1./(1 + 0.05.*sqrt(sum(v_ion(:,2:3).^2, 2))));
            data = [E_MAD, E_mean, ang];
        end

        function region_train(obj, plot_progress)
            disp('Training region classifier...')
            % X1 = readmatrix(['cache/', 'traindata1.txt']);
            X4 = readmatrix(['cache/', 'traindata4.txt']);
            X5 = readmatrix(['cache/', 'traindata5.txt']);
            % train_data2 = readmatrix([cache_path, 'traindata2.txt']);

            % X = [irf_abs(X1(:,2:4), 1), irf_abs(X1(:,5:7), 1), X1(:,8)];
            % X_n = X1(:,8);
            % X_b = irf_abs(X1(:,2:4), 1);
            E_ion_bin_center = X4(:,37:68);
            E_ion_bin_delta = X5(:,2:33);
            E_ion_omni = X4(:,5:36);
            v_ion = X4(:,2:4);
            pos = X4(:,69:71);
            X = obj.create_region_data(E_ion_bin_center, E_ion_bin_delta, E_ion_omni, pos, v_ion);
            [n, ~] = size(X);

            % Normalize data
            m = mean(X);
            s = std(X);
            X_norm = (X - m)./s;

            if plot_progress
                plot_points(X, 'Raw data')
            end

            % Hierarchical clustering.
            tree = linkage(X_norm, 'average', 'euclidean');
            % y = cluster(Z, 'Cutoff', 0.1, 'Criterion', 'distance');

            % Lowest number of clusters that separates MSH, SW and MSP.
            n_clust = 10;
            y = cluster(tree, 'MaxClust', n_clust);
            % y = dbscan(X_norm, 0.15, 100, 'Distance', 'squaredeuclidean');
            % s = [-0.5, 1.7, 1.2;
            %     -1.5, -1, -1.2;
            %     0.4, -0.6, 0.33;
            %     1, -0.2, -1];
            % y = kmeans(X_norm, 4, 'Start', s);
            % y = cluster(tree,'Cutoff',  1.5, 'Criterion', 'distance');
            
            % try
            %     % df
            %     y = readmatrix('y3.txt');
            % catch
            %     % k = ceil(size(X_norm, 1)^0.5);
            %     % SimMatrix = SNN(X_norm, k);
            %     % SimMatrix = k - SimMatrix;
            %     % Eps = 58;
            %     % MinPts = 70;
            %     % y = Mdbscan(X_norm, MinPts, Eps, SimMatrix)';

            %     Ratio=1.7; 
            %     k = ceil(size(X_norm,1)^0.5);
            %     SimMatrix = SNN(X_norm, k);
            %     SimMatrix = k - SimMatrix;
            %     Eps = 57;
            %     % threshold = 0.5796;
            %     threshold = 0.90;
            %     eta = Eps*Ratio;
            %     y = DRSCAN(X_norm, threshold, Eps, eta, SimMatrix)';

            %     unique(y)
            %     writematrix(y, 'y3.txt')
            % end
            if plot_progress
                plot_points(X, 'Hierarchical clusters', y)
            end
            
            % % merge_map = {[2, 6], 3, [5, 1], 4};
            % merge_map = {[1, 3], [2, 4], [5, 7], [6, 8], [9, 10]};
            merge_map = {[3, 8, 9], [5, 7], [1, 4], [2, 6, 10]};
            % % merge_map = {[2, 3], 5, 1, 6, 4, 7};
            y = merge_clusters(y, merge_map);
            % outlier = ismember(y, [-1]);
            % X_no = X(~outlier, :);
            % y = y(~outlier);
            n_clust = length(unique(y));

            if plot_progress % && false
                plot_points(X, 'Merged clusters', y)
            end

            % Get mu, sigma and p of each cluster.
            mu = zeros(n_clust, 3);
            sig = zeros(3, 3, n_clust);
            p = zeros(n_clust, 1);

            for i = 1:n_clust
                logic = y == i;
                X_clust = X(logic, :);
                mu(i, :) = mean(X_clust, 1);
                sig(:, :, i) = cov(X_clust);
                p(i) = nnz(logic)/n;
            end

            % % If cluster is small (p is low) remove it.
            % max(p)
            % idx_small = p < 0.003;
            % mu(idx_small, :) = [];
            % sig(:, :, idx_small) = [];
            % p(idx_small) = [];
            % p

            % Create GMM of remaining clusters.
            model = gmdistribution(mu, sig, p);
            % model.ComponentProportion
            % options = statset('Display', 'final', 'MaxIter', 1500, 'TolFun', 1e-3);
            % rng(101)
            % % S.mu = [
            % %     0.4, 9.5, -0.5;
            % %     0.2, 7, -3.5;
            % %     0.7, 7.2, -1.5;
            % %     0.9, 7.2, -3.5
            % % ];
            % % S.Sigma = repmat(eye(3)*0.4, [1, 1, 4]);
            % % S.ComponentProportion = [1, 2, 2, 1];
            % under = X(:, 3) < -2.7;
            % model = fitgmdist(X(under, :), 2, 'Start', 'randSample', 'Replicates', 1, 'Options', options);
            % rng(11)
            % model = fitgmdist(X, 4, 'Start', 'randSample', 'Replicates', 1, 'Options', options);
            % model = gmdistribution(S.mu, S.Sigma, S.ComponentProportion);
            % if plot_progress
            %     tmp = 1:n_clust;
            %     tmp = tmp(~idx_small);
            %     logic = ismember(y, tmp);
            %     % unique(y(logic))
            %     y = normalize_y(y(logic));
            %     % unique(y)
            %     plot_points(X(logic, :), 'reduced', y)
            % end
            % df
            y = cluster(model, X);


            if plot_progress
                plot_points(X, 'GMM clusters', y)
            end

            % Create a map to merge gaussians which belong to the same class.
            % merge_map = {[2, 3], 1, [4, 5, 6]};
            % y = merge_clusters(y, merge_map);

            % if plot_progress % && false
            %     plot_points(X, 'Merged clusters', y)
            % end

            obj.region_classifier = @classifierFunction;
            obj.model_region = model;

            y = obj.region_classifier(X);

            % [t, probs, mahaldist, logpdf, nlogL] = obj.region_classifier([0.5, 7, 0])

            if plot_progress
                % Create boundary points
                n_bound = 100;
                % [b_x, b_y, b_z] = meshgrid(linspace(4, 10, n_bound), linspace(-400, 800, n_bound), linspace(-1, 1, n_bound));
                % [b_x, b_y, b_z] = meshgrid(linspace(0, 3, n_bound), linspace(4, 10, n_bound), linspace(-1, 1, n_bound));
                [b_x, b_y, b_z] = meshgrid(linspace(0, 2, n_bound), linspace(5, 10, n_bound), linspace(-8, 1, n_bound));
                X_bound = [reshape(b_x, [], 1, 1), reshape(b_y, [], 1, 1), reshape(b_z, [], 1, 1)];
                [~, prob_test, mahal_test, logpdf] = obj.region_classifier(X_bound);
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
                % bound = any(mahal_test > 6 & prob_test < 7, 2);
                % bound = any(prob_test > 0.01 & prob_test < 0.05, 2);
                % bound = any(prob_test > 0.85 & prob_test < 0.95, 2);
                % bound = any(prob_test > 0.15 & prob_test < 0.25, 2);
                % bound = any(prob_test > -0.05 & prob_test < 0.05, 2);
                % bound = all(prob_test > 0.00 & prob_test < 0.05, 2);
                X_bound = X_bound(bound,:);

                plot_points(X, 'Classified raw data', y, X_bound)
            end

            % df

            function [y, probs, mahaldist, logpdf] = classifierFunction(X)
                [n, ~] = size(X);
                % l = length(merge_map);
                % probs = zeros(n, l);
                % mahaldist = zeros(n, l);
                [y, ~, probs, logpdf, d2] = cluster(model, X);
                % for j = 1:l
                %     probs(:,j) = sum(P(:, merge_map{j}), 2);
                %     mahaldist(:,j) = min(d2(:, merge_map{j}), [], 2);
                %     % idx(:,i) = model.ComponentProportion(merge_map{i}(idx));
                % end
                mahaldist = sqrt(d2);
                % y = merge_clusters(y, merge_map);
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

            function y = normalize_y(y)
                l = length(y);
                uniq = unique(y);
                n_uniq = length(uniq);
                idx = zeros(l, n_uniq);
                for k = 1:n_uniq
                    idx(:, k) = y == uniq(k);
                end
                for k = 1:n_uniq
                    y(logical(idx(:, k))) = k;
                end
            end

            function plot_points(X, name, y, X_bound)
                figure;
                if nargin == 2
                    scatter3(X(:,1), X(:,2), X(:,3), 2)
                else
                    cVec = [1 0 0; 0 0.8 0; 0 0 1; 0 1 1; 1 0 1];
                    [l, ~] = size(cVec);
                    uniq = unique(y);
                    if length(uniq) > l || any(uniq == -1)
                        c = y;
                    else
                        c = cVec(y, :);
                    end
                    scatter3(X(:,1), X(:,2), X(:,3), 2, c)
                    if nargin == 4
                        hold on
                        scatter3(X_bound(:,1), X_bound(:,2), X_bound(:,3), 2, [0, 0, 0])
                    end
                end
                title(name)
                xlabel('MAD log energy flux');
                ylabel('Mean log energy flux');
                zlabel('log v angle')
            end
        end

        function [probs, mahaldist, logpdf] = region_classify(obj, r)
            obj.load_ts({'E_ion_center', 'E_ion_delta', 'E_ion_omni', 'pos', 'v_ion'});

            t = obj.ts.E_ion_omni.time;
            pos = obj.ts.pos.resample(obj.ts.E_ion_omni);

            % X = obj.create_region_data(obj.ts.E_ion_center.data, obj.ts.E_ion_omni.data, pos.data);
            X = obj.create_region_data(obj.ts.E_ion_center.data, obj.ts.E_ion_delta.data, obj.ts.E_ion_omni.data, pos.data, obj.ts.v_ion.data);
            % size(X)

            [y, probs, mahaldist, logpdf] = obj.region_classifier(X);
            % if r.smoothWindowSize > 0
            %     probs = movmean(probs, smoothWindowSize*1e9, 1, 'SamplePoints', t.epoch);
            %     logpdf = movmean(logpdf, smoothWindowSize*1e9, 1, 'SamplePoints', t.epoch);
            %     mahaldist = movmean(mahaldist, smoothWindowSize*1e9, 1, 'SamplePoints', t.epoch);
            % end
            probs = TSeries(t, probs);
            mahaldist = TSeries(t, mahaldist);
            logpdf = TSeries(t, logpdf);
            % region = struct('nlogL', nlogL, )
            % prob = v_ion.data(:, 1)./(1 + irf_abs(v_ion.data(:, 2:3), 1));
            % y = double(prob < -3);
            % prob = TSeries(v_ion.time, prob);
            % ys = TSeries(v_ion.time, y);

            ts_label = TSeries(t, y);
            % ts_label = moving_window_func(ts_label, 60, @(x) majorityVote(x.data), t);
            tints_cell_region = labels2tints(ts_label, 1:5);
            obj.last_run_results.tints_cell_region = tints_cell_region;
            % tints_SW = tints_cell_region{2};
            % % tints_SW = remove_short_tints(tints_SW, 60);
            % tints_SW = remove_empty_tints(tints_SW);

            obj.last_run_results.region_plotter = @plotter;
            function plotter(plt)
                plt.lineplot('class_probs', probs, 'ylabel', 'Posterior', 'legend', {'MSP', 'SW', 'MSH', 'x1', 'x2'}, 'colorOrder', [1 0 0; 0 0.8 0; 0 0 1; 0 1 1; 1 0 1])
            end

            function v = majorityVote(x)
                [~, v] = max(accumarray(x,1));
            end
        end
        
        function [SLAMS, logic_SLAMS] = SLAMS_find(obj, r)
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

            b_upper_thresh = r.SLAMS_threshold*obj.ts.b_bg;
            b_lower_thresh = ((r.SLAMS_threshold - 1)/2 + 1)*obj.ts.b_bg;

            logic_above = obj.ts.b_abs.data > b_upper_thresh.data;
            logic_under = obj.ts.b_abs.data < b_lower_thresh.data;
            logic_SLAMS = merge_SLAMS(logic_above, logic_under);

            SLAMS_label = TSeries(t, double(logic_SLAMS));
            
            tints_cell_NOTSLAMS_SLAMS = labels2tints(SLAMS_label, [0,1]);
            tints_SLAMS = tints_cell_NOTSLAMS_SLAMS{2};
            tints_SLAMS = remove_edge_SLAMS(tints_SLAMS, obj.tint);
            obj.last_run_results.tints_SLAMS = tints_SLAMS;
            if ~isempty(tints_SLAMS)
                logic_SLAMS(t < tints_SLAMS(1) | t > tints_SLAMS(end)) = 0;

                n_SLAMS = length(tints_SLAMS)/2;
                vec = repelem(1, n_SLAMS);
                starts = mat2cell(tints_SLAMS(1:2:end), vec, 1);
                stops = mat2cell(tints_SLAMS(2:2:end), vec, 1);
                SLAMS = struct('start', starts, 'stop', stops);
            else
                logic_SLAMS(:) = 0;

                SLAMS = [];
            end

            
            % tints_SLAMS = merge_tints(tints_SLAMS, r.SLAMS_DeltaMaxMerge);
            % tints_SLAMS = remove_short_tints(tints_SLAMS, r.SLAMS_MinDur);
            % tints_SLAMS = remove_empty_tints(tints_SLAMS);
            
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

            function tints_SLAMS = remove_edge_SLAMS(tints_SLAMS, tint_allowed)
                start_idx = find(tints_SLAMS >= tint_allowed(1), 1);
                if mod(start_idx, 2) == 0
                    start_idx = start_idx + 1;
                end
                stop_idx = find(tints_SLAMS <= tint_allowed(2), 1, 'last');
                if mod(stop_idx, 2) == 1
                    stop_idx = stop_idx - 1;
                end
                if start_idx < stop_idx
                    tints_SLAMS = tints_SLAMS(start_idx:stop_idx);
                else
                    tints_SLAMS = [];
                end
            end
        end

        function plot_comparison(obj, tint, tints_true_SLAMS)
            obj.load_ts({'b', 'b_abs', 'v_ion', 'v_ion_abs', 'n_ion', 'E_ion_omni', 'E_ion_center', 'E_ion_delta', 'pos'})

            plt = modular_plot('title', 'Predicted vs true');
            if isfield(obj.last_run_results, 'region_plotter')
                obj.last_run_results.region_plotter(plt);
            end
            plt.lineplot('B_pred', obj.ts.b_abs);
            if isfield(obj.last_run_results, 'SLAMS_plotter')
                obj.last_run_results.SLAMS_plotter(plt, 'B_pred', {'ylabel', 'Pred'})
            end
            plt.lineplot('B_true', obj.ts.b_abs, 'ylabel', 'True');
            plt.lineplot('B', {obj.ts.b, obj.ts.b_abs}, 'ylabel', 'B_{GSE} (nT)');
            plt.lineplot('n_ion', obj.ts.n_ion, 'ylabel', 'n_i (cm^{-3})')
            % plt.lineplot('T_ion', {T_ion_para, T_ion_perp}, 'ylabel', 'T_{i} (eV)', 'legend', {'T_{ipara}','T_{iperp}'})
            plt.spectrogram('E_ion_omni', {obj.ts.E_ion_omni, obj.ts.E_ion_center}, 'ylabel', 'W_i (eV)')
            % plt.spectrogram('E_ion_omni', {obj.ts.E_ion_omni.*obj.ts.E_ion_delta, obj.ts.E_ion_center}, 'ylabel', 'W_i (eV)')
            plt.lineplot('v_ion', {obj.ts.v_ion, obj.ts.v_ion_abs}, 'ylabel', 'v_{i} (kms^{-1})', 'legend', {'v_{ix}','v_{iy}','v_{iz}'})
            plt.lineplot('pos', obj.ts.pos/6378, 'ylabel', 'R_{GSE} (R_E)', 'legend', {'x','y','z'})
            plt.show(tint);
            plt.mark(obj.last_run_results.tints_SLAMS, 'target', 'B_pred', 'intervals', true);
            plt.mark(tints_true_SLAMS, 'target', 'B_true', 'intervals', true);
            if isfield(obj.last_run_results, 'region_plotter')
                plt.mark(obj.last_run_results.tints_cell_region, 'target', 'class_probs', 'intervals', true);
            end
        end

        function plot_prediction(obj, tint)
            obj.load_ts({'b', 'b_abs', 'v_ion', 'v_ion_abs', 'n_ion', 'E_ion_omni', 'E_ion_center', 'pos'})

            plt = modular_plot('title', 'Predicted SLAMS');
            if isfield(obj.last_run_results, 'region_plotter')
                obj.last_run_results.region_plotter(plt);
            end
            plt.lineplot('B_pred', obj.ts.b_abs);
            if isfield(obj.last_run_results, 'SLAMS_plotter')
                obj.last_run_results.SLAMS_plotter(plt, 'B_pred', {'ylabel', 'Pred'})
            end
            plt.lineplot('B', {obj.ts.b, obj.ts.b_abs}, 'ylabel', 'B_{GSE} (nT)');
            plt.lineplot('n_ion', obj.ts.n_ion, 'ylabel', 'n_i (cm^{-3})')
            % plt.lineplot('E_s', obj.ts.E_ion_omni, 'ylabel', 'spectr')
            % plt.lineplot('E', obj.ts.E_ion_center, 'ylabel', 'energy')
            % plt.lineplot('T_ion', {T_ion_para, T_ion_perp}, 'ylabel', 'T_{i} (eV)', 'legend', {'T_{ipara}','T_{iperp}'})
            plt.spectrogram('E_ion_omni', {obj.ts.E_ion_omni, obj.ts.E_ion_center}, 'ylabel', 'W_i (eV)')

            % E_bins = obj.ts.E_ion_center.data;
            % % E_bins = (1:32)/32;
            % tot = sum(obj.ts.E_ion_omni.data, 2);
            % X_e_m = sum(obj.ts.E_ion_omni.data.*E_bins, 2)./tot;
            % X_e_v = sum(obj.ts.E_ion_omni.data.*((E_bins - X_e_m).^6), 2)./tot;
            % test1 = TSeries(obj.ts.E_ion_omni.time, X_e_m);
            % test2 = TSeries(obj.ts.E_ion_omni.time, X_e_v);
            % plt.lineplot('test1', test1, 'ylabel', 'test1')
            % plt.lineplot('test2', test2, 'ylabel', 'test2')

            plt.lineplot('v_ion', {obj.ts.v_ion, obj.ts.v_ion_abs}, 'ylabel', 'v_{i} (kms^{-1})', 'legend', {'v_{ix}','v_{iy}','v_{iz}'})
            plt.lineplot('pos', obj.ts.pos/6378, 'ylabel', 'R_{GSE} (R_E)', 'legend', {'x','y','z'})
            plt.show(tint);
            plt.mark(obj.last_run_results.tints_SLAMS, 'target', 'B_pred', 'intervals', true);
            if isfield(obj.last_run_results, 'region_plotter')
                plt.mark(obj.last_run_results.tints_cell_region, 'target', 'class_probs', 'intervals', true);
            end
        end

        function newRun(obj, tint, extra_load_time)
            obj.ts = rmfield(obj.ts, fieldnames(obj.ts));
            obj.last_run_results = rmfield(obj.last_run_results, fieldnames(obj.last_run_results));

            obj.tint = tint;
            obj.tint_load = [tint(1) + -extra_load_time, tint(2) + extra_load_time];
        end

        function SLAMS = evaluate(obj, tint, varargin)

            p = inputParser;
            addRequired(p, 'tint');
            addParameter(p, 'Include_B_stats', true)
            addParameter(p, 'Include_region_stats', true)
            addParameter(p, 'Region_time_windows', [])
            addParameter(p, 'Include_GSE_coords', true)
            addParameter(p, 'Extra_load_time', 60)
            addParameter(p, 'SLAMS_B_bg_method', 'median')
            addParameter(p, 'SLAMS_B_bg_window', 60)
            addParameter(p, 'SLAMS_threshold', 2)
            parse(p, tint, varargin{:})
            r = p.Results;

            obj.newRun(r.tint, r.Extra_load_time)

            [SLAMS, logic_SLAMS] = obj.SLAMS_find(r);
            if r.Include_region_stats
                [probs, mahaldist, logpdf] = obj.region_classify(r);
            end
            
            if ~isempty(SLAMS) && false
                n_SLAMS = length(SLAMS);

                if r.Include_GSE_coords || r.Include_region_stats
                    starts = [SLAMS.start];
                    stops = [SLAMS.stop];
                    SLAMS_mid = (starts.epoch + stops.epoch)/2;
                    cell_converter_tool = repelem(1, n_SLAMS);
                    
                    if r.Include_GSE_coords
                        obj.load_ts({'pos'})
                        pos = interp1(double(obj.ts.pos.time.epoch), double(obj.ts.pos.data), double(SLAMS_mid));
                        pos_cell = mat2cell(pos, cell_converter_tool, 3);
                        [SLAMS.pos_GSE] = pos_cell{:};
                    end
                    
                    if r.Include_region_stats
                        t = probs.time.epoch;
                        [~, idx_closest] = min(abs(SLAMS_mid - t'), [], 2);
                        probs_cell = mat2cell(probs.data(idx_closest, :), cell_converter_tool, 3);
                        mahaldist_cell = mat2cell(mahaldist.data(idx_closest, :), cell_converter_tool, 3);
                        logpdf_cell = mat2cell(logpdf.data(idx_closest, :), cell_converter_tool, 1);

                        [SLAMS.region_posterior] = probs_cell{:};
                        [SLAMS.region_mahaldist] = mahaldist_cell{:};
                        [SLAMS.region_logpdf] = logpdf_cell{:};

                        if ~isempty(r.Region_time_windows)
                            % SLAMS.region_windows = r.Region_time_windows;
                            n_windows = length(r.Region_time_windows);
                            dt = int64((1e9*r.Region_time_windows)/2);
                            SLAMS_window_start = SLAMS_mid - dt;
                            SLAMS_window_stop = SLAMS_mid + dt;
                            t_3d = reshape(t, 1, 1, []);
                            logical_after_start = SLAMS_window_start < t_3d;
                            logical_before_stop = SLAMS_window_stop > t_3d;
                            logical_inside_window = logical_after_start & logical_before_stop;

                            n_inside_window = sum(logical_inside_window, 3);
                            probs3d = mean3d(probs.data, n_inside_window);
                            mahaldist3d = mean3d(mahaldist.data, n_inside_window);
                            logpdf3d = log(permute(sum(logical_inside_window.*reshape(exp(logpdf.data), 1, 1, []), 3)./n_inside_window, [2, 3, 1]));

                            probs3d_cell = mat2cell(probs3d, n_windows, 3, cell_converter_tool);
                            mahaldist3d_cell = mat2cell(mahaldist3d, n_windows, 3, cell_converter_tool);
                            logpdf3d_cell = mat2cell(logpdf3d, n_windows, 1, cell_converter_tool);
                            
                            [SLAMS.region_posterior_windows] = probs3d_cell{:};
                            [SLAMS.region_mahaldist_windows] = mahaldist3d_cell{:};
                            [SLAMS.region_logpdf_windows] = logpdf3d_cell{:};
                        end
                    end
                end
                
                if r.Include_B_stats
                    idx_first_last = sort([find(logical_manip(logic_SLAMS', 'firstOfOne')), find(logical_manip(logic_SLAMS', 'lastOfOne'))]);
                    for j = 1:n_SLAMS
                        idx_range = idx_first_last(j*2 - 1):idx_first_last(j*2);
                        B_abs = obj.ts.b_abs.data(idx_range);
                        B_bg = obj.ts.b_bg.data(idx_range);
                        SLAMS(j).B_bg_mean = mean(B_bg);
                        SLAMS(j).B_mean = mean(B_abs);
                        SLAMS(j).B_max = max(B_abs);
                        B_rel = B_abs./B_bg;
                        SLAMS(j).B_rel_mean = mean(B_rel);
                        SLAMS(j).B_rel_max = max(B_rel);
                    end
                end
            end
            % tints_SLAMS_SW = intersect_tints(tints_SW, tints_SLAMS);
            % tints_SLAMS = intersect_tints(tints_SW, tints_SLAMS);
            % tints_SLAMS = remove_edge_SLAMS(tints_SLAMS_SW, tints_SLAMS);

            % tints_SLAMS = intersect_tints(tints_SLAMS, r.tint);
            % tints_SLAMS = remove_empty_tints(tints_SLAMS);
            function out = mean3d(inp, n_inside_window)
                a = reshape(inp(:,1), 1, 1, []);
                b = reshape(inp(:,2), 1, 1, []);
                c = reshape(inp(:,3), 1, 1, []);
                a = sum(logical_inside_window.*a, 3);
                b = sum(logical_inside_window.*b, 3);
                c = sum(logical_inside_window.*c, 3);
                out = cat(3, a, b, c)./n_inside_window;
                out = permute(out, [2, 3, 1]);
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
                case 'E_ion_center'
                    if ~isfield(obj.ts, 'E_ion_center')
                        obj.ts.E_ion_center = load_tints_MMS('mms1_fpi_fast_l2_dis-moms', 'mms1_dis_energy_fast', obj.tint_load);
                    end
                case 'E_ion_delta'
                    if ~isfield(obj.ts, 'E_ion_delta')
                        obj.ts.E_ion_delta = load_tints_MMS('mms1_fpi_fast_l2_dis-moms', 'mms1_dis_energy_delta_fast', obj.tint_load);
                    end
                case 'pos'
                    if ~isfield(obj.ts, 'pos')
                        obj.ts.pos = load_tints_MMS('mms1_mec_srvy_l2_ephts04d', 'mms1_mec_r_gse', obj.tint_load);
                    end
                otherwise
                    error(['''' cell_string{i} ''' loader not implemented.'])
                end
            end
        end
    end
end

