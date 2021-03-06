classdef SLAMS_finder < handle
    %SLAMS_FINDER Finds SLAMS
    %   Instatiate this class and then run the method 'evaluate' to find SLAMS.
    %   Check the constructor and the 'evaluate' method for valid name value pairs.
    
    properties
        region_classifier
        classes
        n_classes
        priors
        ts
        last_run_results
        tint_load
        tint
        sc
        search_durations_classes
        search_durations
    end
    
    methods
        function obj = SLAMS_finder(varargin)

            p = inputParser;
            addParameter(p, 'Show_region_classifier_steps', false)
            addParameter(p, 'Spacecraft', 'MMS1')
            parse(p, varargin{:})
            r = p.Results;

            obj.ts = struct;
            obj.last_run_results = struct;
            obj.construct_region_classifier(r.Show_region_classifier_steps);
            obj.sc = lower(r.Spacecraft);
            obj.search_durations = containers.Map;
        end

        function data = create_region_data(~, E_ion_bin_center, E_ion_bin_delta, E_ion_omni, v_ion)
            E_bins_center = log(E_ion_bin_center);
            E_per_bin = E_ion_omni.*E_ion_bin_delta;
            E_tot = sum(E_per_bin, 2);
            E_mean = sum(E_per_bin.*E_bins_center, 2)./E_tot;
            E_MAD = sum(E_per_bin.*(abs(E_bins_center - E_mean)), 2)./E_tot;
            ang = acos(-v_ion(:,1)./sqrt(sum(v_ion.^2, 2)))/pi;
            ang = max(log(ang), -7);
            data = [E_MAD, E_mean, ang];
        end

        function construct_region_classifier(obj, plot_progress)
            disp('Constructing region classifier...')
            % X4 = readmatrix(['cache/', 'traindata4.txt']);
            % X5 = readmatrix(['cache/', 'traindata5.txt']);

            X = readmatrix(['cache/', 'randSampData.txt']);

            v_ion = X(:,2:4);
            E_ion_omni = X(:,5:36);
            E_ion_bin_delta = X(:,37:68);
            E_ion_bin_center = X(:,69:100);
            X = obj.create_region_data(E_ion_bin_center, E_ion_bin_delta, E_ion_omni, v_ion);

            [n, ~] = size(X);

            if plot_progress
                plot_points(X, 'Raw data')
            end

            % Normalize data
            m = mean(X);
            s = std(X);
            X_norm = (X - m)./s;

            % Hierarchical clustering
            tree = linkage(X_norm, 'ward', 'euclidean');

            % Lowest number of clusters that separates MSP, SW, MSH and FS
            n_clust = 9;
            y = cluster(tree, 'MaxClust', n_clust);
            
            if plot_progress
                plot_points(X, 'Hierarchical clusters', y)
            end
            
            % Merge clusters which belong to same class
            merge_map = {[2, 8], 9, [3, 4, 5, 6], [7, 8], 1};
            y = merge_clusters(y, merge_map);
            
            if plot_progress
                plot_points(X, 'Merged clusters', y)
            end

            % Discard last cluster
            n_clust = length(unique(y)) - 1;
            obj.n_classes = n_clust;
            obj.classes = {'MSP', 'SW', 'MSH', 'FS'};

            % Get mu, sigma and p of each cluster
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

            % Create GMM of each cluster
            model = gmdistribution(mu, sig, p);
            obj.priors = model.ComponentProportion;

            % Classify data using GMM model
            y = cluster(model, X);

            if plot_progress
                plot_points(X, 'GMM clusters', y)
            end

            obj.region_classifier = @classifierFunction;

            if plot_progress
                % Create boundary points
                n_bound = 100;
                [b_x, b_y, b_z] = meshgrid(linspace(0, 2, n_bound), linspace(5, 10, n_bound), linspace(-8, 1, n_bound));
                X_bound = [reshape(b_x, [], 1, 1), reshape(b_y, [], 1, 1), reshape(b_z, [], 1, 1)];
                [~, prob_test, ~, ~] = obj.region_classifier(X_bound);
                bound = any(prob_test > 0.45 & prob_test < 0.55, 2);
                X_bound = X_bound(bound,:);
                
                y = obj.region_classifier(X);

                plot_points(X, 'GMM clusters with boundary', y, X_bound)
            end

            function [y, posteriors, mahaldists, logpdfs] = classifierFunction(X)
                [y, ~, posteriors, logpdfs, d2] = cluster(model, X);
                % idx = find(isnan(logpdfs))
                % X(idx, :)
                mahaldists = sqrt(d2);
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
                % MAD log particle energy
                xlabel('F_3');
                % Mean log particle energy
                ylabel('F_2');
                % Log normalized cone angle
                zlabel('F_1');
            end
        end

        function [posteriors, mahaldists, logpdfs] = region_classify(obj)
            obj.load_ts({'E_ion_center', 'E_ion_delta', 'E_ion_omni', 'v_ion'});

            t = obj.ts.E_ion_omni.time;

            X = obj.create_region_data(obj.ts.E_ion_center.data, obj.ts.E_ion_delta.data, obj.ts.E_ion_omni.data, obj.ts.v_ion.data);

            [y, posteriors, mahaldists, logpdfs] = obj.region_classifier(X);

            posteriors = TSeries(t, posteriors);
            mahaldists = TSeries(t, mahaldists);
            logpdfs = TSeries(t, logpdfs);

            ts_label = TSeries(t, y);
            tints_cell_region = labels2tints(ts_label, 1:obj.n_classes);
            obj.last_run_results.tints_cell_region = tints_cell_region;

            obj.last_run_results.region_plotter = @plotter;
            function plotter(plt)
                plt.lineplot('class_probs', posteriors, 'ylabel', 'Posterior', 'legend', obj.classes, 'colorOrder', [1 0 0; 0 0.7 0; 0 0 1; 0 0.7 0.7])
            end

            function v = majorityVote(x)
                [~, v] = max(accumarray(x,1));
            end
        end

        function track_search_durations(obj, r)
            obj.load_ts({'pos'})

            pos = tlim(obj.ts.pos, obj.tint);
            t = pos.time;
            pos = pos.data/6378;
            BS = GSE2BS(pos);
            xBS_disc = ceil(BS(:,1));
            yBS_disc = ceil(BS(:,2));

            tmp = xBS_disc(2:end) ~= xBS_disc(1:end-1);
            starts_xBS = [1; tmp];
            stops_xBS = [tmp; 1];
            tmp = yBS_disc(2:end) ~= yBS_disc(1:end-1);
            starts_yBS = [1; tmp];
            stops_yBS = [tmp; 1];
            starts = find(starts_xBS | starts_yBS);
            stops = find(stops_xBS | stops_yBS);

            n_tints = length(starts);
            for i = 1:n_tints
                t_start = t(starts(i));
                t_stop = t(stops(i));
                dur = t_stop - t_start;
                key = [num2str(xBS_disc(starts(i))), ',', num2str(yBS_disc(starts(i)))];
                if r.Include_region_stats
                    arr_new = zeros(1, obj.n_classes + 1);
                    arr_new(end) = dur;
                    for j = 1:obj.n_classes 
                        class_tints = intersect_tints(obj.last_run_results.tints_cell_region{j}, [t_start, t_stop]);
                        dif = diff(class_tints.epoch);
                        dur_class = sum(dif(1:2:end)./1e9);
                        arr_new(j) = dur_class;
                    end
                    if isKey(obj.search_durations, key)
                        arr = obj.search_durations(key);
                        arr = arr + arr_new;
                        obj.search_durations(key) = arr;
                    else
                        obj.search_durations(key) = arr_new;
                    end
                else
                    if isKey(obj.search_durations, key)
                        obj.search_durations(key) = obj.search_durations(key) + dur;
                    else
                        obj.search_durations(key) = dur;
                    end
                end
            end
        end
        
        function SLAMS = SLAMS_find(obj, r)
            obj.load_ts({'b', 'b_abs'})

            t = obj.ts.b_abs.time;

            switch r.SLAMS_B0_method
            case 'mean'
                obj.ts.b0 = TSeries(t, movmean(obj.ts.b_abs.data, r.SLAMS_B0_window*1e9, 1, 'SamplePoints', obj.ts.b_abs.time.epoch));
            case 'median'
                obj.ts.b0 = TSeries(t, movmedian(obj.ts.b_abs.data, r.SLAMS_B0_window*1e9, 1, 'SamplePoints', obj.ts.b_abs.time.epoch));
            case 'harmmean'
                obj.ts.b0 = moving_window_func(obj.ts.b_abs, r.SLAMS_B0_window, @(x) harmmean(x.data), t);
            end

            b_upper_thresh = r.SLAMS_detection_constant*obj.ts.b0;
            b_lower_thresh = r.SLAMS_merging_constant*obj.ts.b0;

            logic_above = obj.ts.b_abs.data > b_upper_thresh.data;
            logic_under = obj.ts.b_abs.data < b_lower_thresh.data;
            logic_SLAMS = merge_SLAMS(logic_above, logic_under);

            SLAMS_label = TSeries(t, double(logic_SLAMS));
            
            tints_cell_NOTSLAMS_SLAMS = labels2tints(SLAMS_label, [0,1]);
            tints_SLAMS = tints_cell_NOTSLAMS_SLAMS{2};
            tints_SLAMS = remove_edge_SLAMS(tints_SLAMS, obj.tint);
            tints_SLAMS = remove_short_tints(tints_SLAMS, r.SLAMS_min_duration);
            obj.last_run_results.tints_SLAMS = tints_SLAMS;
            if ~isempty(tints_SLAMS)
                n_SLAMS = length(tints_SLAMS)/2;
                vec = repelem(1, n_SLAMS);
                starts = mat2cell(tints_SLAMS(1:2:end), vec, 1);
                stops = mat2cell(tints_SLAMS(2:2:end), vec, 1);
                SLAMS = struct('start', starts, 'stop', stops);
            else
                SLAMS = [];
            end
            
            obj.last_run_results.SLAMS_plotter = @plotter;
            function plotter(plt, name, arg)
                plt.lineplot(name, obj.ts.b0, 'color', 'c');
                plt.lineplot(name, b_lower_thresh, 'color', 'm');
                plt.lineplot(name, b_upper_thresh, 'color', 'r', arg{:});
            end

            function logic_SLAMS = merge_SLAMS(logic_above, logic_under)
                if all(logic_above == 0)
                    logic_SLAMS = logic_above;
                    return
                end
                index_pair_not_under = sort([find(logical_manip(logic_under, 'firstOfZero')); find(logical_manip(logic_under, 'lastOfZero'))]);
                index_above = find(logical_manip(logic_above, 'firstOfOne'));
                for i = 1:2:length(index_pair_not_under)
                    start_idx = index_pair_not_under(i);
                    stop_idx = index_pair_not_under(i + 1);
                    if all(~((start_idx <= index_above) & (stop_idx >= index_above)))
                        logic_under(start_idx:stop_idx) = true;
                    end
                end
                logic_SLAMS = ~logic_under;
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

            plt = modular_plot('title', [upper(obj.sc), ', automatic vs manual SLAMS detection']);
            if isfield(obj.last_run_results, 'region_plotter')
                obj.last_run_results.region_plotter(plt);
            end
            plt.lineplot('B_pred', obj.ts.b_abs);
            if isfield(obj.last_run_results, 'SLAMS_plotter')
                obj.last_run_results.SLAMS_plotter(plt, 'B_pred', {'ylabel', 'Automatic'})
            end
            plt.lineplot('B_true', obj.ts.b_abs, 'ylabel', 'Manual');
            plt.lineplot('B', {obj.ts.b, obj.ts.b_abs}, 'ylabel', 'B_{GSE} (nT)');
            plt.lineplot('n_ion', obj.ts.n_ion, 'ylabel', 'n_i (cm^{-3})')
            plt.spectrogram('E_ion_omni', {obj.ts.E_ion_omni, obj.ts.E_ion_center}, 'ylabel', 'W_i (eV)')
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

            plt = modular_plot('title', [upper(obj.sc), ', Predicted SLAMS']);
            if isfield(obj.last_run_results, 'region_plotter')
                obj.last_run_results.region_plotter(plt);
            end
            plt.lineplot('B_pred', obj.ts.b_abs);
            if isfield(obj.last_run_results, 'SLAMS_plotter')
                obj.last_run_results.SLAMS_plotter(plt, 'B_pred', {'ylabel', 'Pred'})
            end
            plt.lineplot('B', {obj.ts.b, obj.ts.b_abs}, 'ylabel', 'B_{GSE} (nT)');
            plt.lineplot('n_ion', obj.ts.n_ion, 'ylabel', 'n_i (cm^{-3})')
            plt.spectrogram('E_ion_omni', {obj.ts.E_ion_omni, obj.ts.E_ion_center}, 'ylabel', 'W_i (eV)')
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
            addParameter(p, 'Track_search_durations', true)
            addParameter(p, 'Extra_load_time', 30)
            addParameter(p, 'SLAMS_B0_method', 'median')
            addParameter(p, 'SLAMS_B0_window', 60)
            addParameter(p, 'SLAMS_detection_constant', 2)
            addParameter(p, 'SLAMS_merging_constant', 1.5)
            addParameter(p, 'SLAMS_min_duration', 0)
            parse(p, tint, varargin{:})
            r = p.Results;

            obj.newRun(r.tint, r.Extra_load_time)

            SLAMS = obj.SLAMS_find(r);
            if r.Include_region_stats
                [posteriors, mahaldists, logpdfs] = obj.region_classify();
                % any(isnan(posteriors.data))
                % any(isnan(mahaldists.data))
                % any(isnan(logpdfs.data))
            end

            if r.Track_search_durations
                obj.track_search_durations(r);
            end
            
            if ~isempty(SLAMS)
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
                        t = posteriors.time.epoch;
                        [~, idx_closest] = min(abs(SLAMS_mid - t'), [], 2);
                        SLAMS_posteriors = permute(posteriors.data(idx_closest, :), [3, 2, 1]);
                        SLAMS_mahaldists = permute(mahaldists.data(idx_closest, :), [3, 2, 1]);
                        SLAMS_logpdfs = permute(logpdfs.data(idx_closest, :), [3, 2, 1]);

                        n_sets = 1;

                        if ~isempty(r.Region_time_windows)
                            n_windows = length(r.Region_time_windows);
                            n_sets = n_sets + n_windows;
                            dt = int64((1e9*r.Region_time_windows)/2);
                            window_start = SLAMS_mid - dt;
                            window_stop = SLAMS_mid + dt;
                            % t_3d = reshape(t, 1, 1, []);
                            t_3d = permute(t, [3, 2, 1]);
                            logical_after_start = window_start < t_3d;
                            logical_before_stop = window_stop > t_3d;
                            logical_inside_window = logical_after_start & logical_before_stop;
                            n_inside_window = sum(logical_inside_window, 3);

                            SLAMS_posteriors = vertcat(SLAMS_posteriors, mean3d(posteriors.data, logical_inside_window, n_inside_window));
                            SLAMS_mahaldists = vertcat(SLAMS_mahaldists, mean3d(mahaldists.data, logical_inside_window, n_inside_window));
                            SLAMS_logpdfs = vertcat(SLAMS_logpdfs, log(mean3d(exp(logpdfs.data), logical_inside_window, n_inside_window)));
                        end

                        SLAMS_posteriors = mat2cell(SLAMS_posteriors, n_sets, obj.n_classes, cell_converter_tool);
                        SLAMS_mahaldists = mat2cell(SLAMS_mahaldists, n_sets, obj.n_classes, cell_converter_tool);
                        SLAMS_logpdfs = mat2cell(SLAMS_logpdfs, n_sets, 1, cell_converter_tool);

                        [SLAMS.region_posterior] = SLAMS_posteriors{:};
                        [SLAMS.region_mahaldist] = SLAMS_mahaldists{:};
                        [SLAMS.region_logpdf] = SLAMS_logpdfs{:};
                    end
                end
                
                if r.Include_B_stats
                    t = obj.ts.b_abs.time.epoch;
                    for j = 1:n_SLAMS
                        logical_idx = (SLAMS(j).start.epoch <= t) & (SLAMS(j).stop.epoch >= t);
                        B_abs = obj.ts.b_abs.data(logical_idx);
                        B0 = obj.ts.b0.data(logical_idx);
                        SLAMS(j).B0_mean = mean(B0);
                        SLAMS(j).B_mean = mean(B_abs);
                        SLAMS(j).B_max = max(B_abs);
                        B_rel = B_abs./B0;
                        SLAMS(j).B_rel_mean = mean(B_rel);
                        SLAMS(j).B_rel_max = max(B_rel);
                    end
                end
            end

            function out = mean3d(inp, logical_inside_window, n_inside_window)
                [~, n] = size(inp);
                c = cell(1, n);
                for i = 1:n
                    tmp = permute(inp(:, i), [3, 2, 1]);
                    c{i} = sum(logical_inside_window.*tmp, 3);
                end
                out = cat(3, c{:})./n_inside_window;
                out = permute(out, [2, 3, 1]);
            end
        end

        function load_ts(obj, cell_string)
            for i = 1:length(cell_string)
                switch cell_string{i}
                case 'b'
                    if ~isfield(obj.ts, 'b')
                        obj.ts.b = load_tints_MMS([obj.sc, '_fgm_srvy_l2'], [obj.sc, '_fgm_b_gse_srvy_l2'], obj.tint_load);
                    end
                case 'b_abs'
                    if ~isfield(obj.ts, 'b_abs')
                        obj.ts.b_abs = irf_abs(obj.ts.b);
                    end
                case 'v_ion'
                    if ~isfield(obj.ts, 'v_ion')
                        obj.ts.v_ion = load_tints_MMS([obj.sc, '_fpi_fast_l2_dis-moms'], [obj.sc, '_dis_bulkv_gse_fast'], obj.tint_load);
                    end
                case 'v_ion_abs'
                    if ~isfield(obj.ts, 'v_ion_abs')
                        obj.ts.v_ion_abs = irf_abs(obj.ts.v_ion);
                    end
                case 'n_ion'
                    if ~isfield(obj.ts, 'n_ion')
                        obj.ts.n_ion = load_tints_MMS([obj.sc, '_fpi_fast_l2_dis-moms'], [obj.sc, '_dis_numberdensity_fast'], obj.tint_load);
                    end
                case 'E_ion_omni'
                    if ~isfield(obj.ts, 'E_ion_omni')
                        obj.ts.E_ion_omni = load_tints_MMS([obj.sc, '_fpi_fast_l2_dis-moms'], [obj.sc, '_dis_energyspectr_omni_fast'], obj.tint_load);
                    end
                case 'E_ion_center'
                    if ~isfield(obj.ts, 'E_ion_center')
                        obj.ts.E_ion_center = load_tints_MMS([obj.sc, '_fpi_fast_l2_dis-moms'], [obj.sc, '_dis_energy_fast'], obj.tint_load);
                    end
                case 'E_ion_delta'
                    if ~isfield(obj.ts, 'E_ion_delta')
                        obj.ts.E_ion_delta = load_tints_MMS([obj.sc, '_fpi_fast_l2_dis-moms'], [obj.sc, '_dis_energy_delta_fast'], obj.tint_load);
                    end
                case 'pos'
                    if ~isfield(obj.ts, 'pos')
                        obj.ts.pos = load_tints_MMS([obj.sc, '_mec_srvy_l2_ephts04d'], [obj.sc, '_mec_r_gse'], obj.tint_load);
                    end
                otherwise
                    error(['''' cell_string{i} ''' loader not implemented.'])
                end
            end
        end
    end
end

