classdef SLAMS_finder < handle
    %SLAMS_FINDER Finds SLAMS
    %   Instatiate this class and then run the method evaluate to find SLAMS.
    %   Check the constructor and the evaluate method for valid name value pairs.
    
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
    end
    
    methods
        function obj = SLAMS_finder(varargin)

            p = inputParser;
            addParameter(p, 'Show_region_classifier_steps', false)
            addParameter(p, 'Show_region_classifier_SLAMS', false)
            addParameter(p, 'Spacecraft', 'MMS1')
            parse(p, varargin{:})
            r = p.Results;

            obj.ts = struct;
            obj.last_run_results = struct;
            obj.construct_region_classifier(r.Show_region_classifier_steps, r.Show_region_classifier_SLAMS);
            obj.sc = lower(r.Spacecraft);
        end

        function data = create_region_data(~, E_ion_bin_center, E_ion_bin_delta, E_ion_omni, v_ion)
            E_bins_center = log(E_ion_bin_center);
            E_per_bin = E_ion_omni.*E_ion_bin_delta;
            E_tot = sum(E_per_bin, 2);
            E_mean = sum(E_per_bin.*E_bins_center, 2)./E_tot;
            E_MAD = sum(E_per_bin.*(abs(E_bins_center - E_mean)), 2)./E_tot;
            ang = acos(-v_ion(:,1)./sqrt(sum(v_ion.^2, 2)))/pi;
            ang = log(ang);
            data = [E_MAD, E_mean, ang];
        end

        function construct_region_classifier(obj, plot_progress, include_SLAMS)
            disp('Constructing region classifier...')
            % X4 = readmatrix(['cache/', 'traindata4.txt']);
            % X5 = readmatrix(['cache/', 'traindata5.txt']);

            X = readmatrix(['cache/', 'randSampData.txt']);

            v_ion = X(:,2:4);
            E_ion_omni = X(:,5:36);
            E_ion_bin_delta = X(:,37:68);
            E_ion_bin_center = X(:,69:100);
            X = obj.create_region_data(E_ion_bin_center, E_ion_bin_delta, E_ion_omni, v_ion);

            if include_SLAMS
                X_SLAMS = readmatrix(['cache/', 'SLAMSSampData.txt']);

                v_ion = X_SLAMS(:,2:4);
                E_ion_omni = X_SLAMS(:,5:36);
                E_ion_bin_delta = X_SLAMS(:,37:68);
                E_ion_bin_center = X_SLAMS(:,69:100);
                X_SLAMS = obj.create_region_data(E_ion_bin_center, E_ion_bin_delta, E_ion_omni, v_ion);
            end
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
                if include_SLAMS
                    plot_points(X, 'GMM clusters', y, X_SLAMS)
                else
                    plot_points(X, 'GMM clusters', y)
                end
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

            function [y, probs, mahaldist, logpdf] = classifierFunction(X)
                [y, ~, probs, logpdf, d2] = cluster(model, X);
                mahaldist = sqrt(d2);
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

        function [probs, mahaldist, logpdf] = region_classify(obj)
            obj.load_ts({'E_ion_center', 'E_ion_delta', 'E_ion_omni', 'v_ion'});

            t = obj.ts.E_ion_omni.time;

            X = obj.create_region_data(obj.ts.E_ion_center.data, obj.ts.E_ion_delta.data, obj.ts.E_ion_omni.data, obj.ts.v_ion.data);

            [y, probs, mahaldist, logpdf] = obj.region_classifier(X);

            probs = TSeries(t, probs);
            mahaldist = TSeries(t, mahaldist);
            logpdf = TSeries(t, logpdf);

            ts_label = TSeries(t, y);
            tints_cell_region = labels2tints(ts_label, 1:obj.n_classes);
            obj.last_run_results.tints_cell_region = tints_cell_region;

            obj.last_run_results.region_plotter = @plotter;
            function plotter(plt)
                plt.lineplot('class_probs', probs, 'ylabel', 'Posterior', 'legend', obj.classes, 'colorOrder', [1 0 0; 0 0.7 0; 0 0 1; 0 0.7 0.7])
            end

            function v = majorityVote(x)
                [~, v] = max(accumarray(x,1));
            end
        end
        
        function SLAMS = SLAMS_find(obj, r)
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
                plt.lineplot(name, obj.ts.b_bg, 'color', 'c');
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

            plt = modular_plot('title', [upper(obj.sc), ', automatic vs manual SLAMS identification']);
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
            addParameter(p, 'Extra_load_time', 60)
            addParameter(p, 'SLAMS_B_bg_method', 'median')
            addParameter(p, 'SLAMS_B_bg_window', 60)
            addParameter(p, 'SLAMS_threshold', 2)
            addParameter(p, 'SLAMS_min_duration', 0)
            parse(p, tint, varargin{:})
            r = p.Results;

            obj.newRun(r.tint, r.Extra_load_time)

            SLAMS = obj.SLAMS_find(r);
            if r.Include_region_stats
                [probs, mahaldist, logpdf] = obj.region_classify();
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
                        t = probs.time.epoch;
                        [~, idx_closest] = min(abs(SLAMS_mid - t'), [], 2);
                        probs_cell = mat2cell(probs.data(idx_closest, :), cell_converter_tool, obj.n_classes);
                        mahaldist_cell = mat2cell(mahaldist.data(idx_closest, :), cell_converter_tool, obj.n_classes);
                        logpdf_cell = mat2cell(logpdf.data(idx_closest, :), cell_converter_tool, 1);

                        [SLAMS.region_posterior] = probs_cell{:};
                        [SLAMS.region_mahaldist] = mahaldist_cell{:};
                        [SLAMS.region_logpdf] = logpdf_cell{:};

                        if ~isempty(r.Region_time_windows)
                            n_windows = length(r.Region_time_windows);
                            dt = int64((1e9*r.Region_time_windows)/2);
                            SLAMS_window_start = SLAMS_mid - dt;
                            SLAMS_window_stop = SLAMS_mid + dt;
                            t_3d = reshape(t, 1, 1, []);
                            logical_after_start = SLAMS_window_start < t_3d;
                            logical_before_stop = SLAMS_window_stop > t_3d;
                            logical_inside_window = logical_after_start & logical_before_stop;

                            n_inside_window = sum(logical_inside_window, 3);
                            probs3d = mean3d(probs.data, logical_inside_window, n_inside_window);
                            mahaldist3d = mean3d(mahaldist.data, logical_inside_window, n_inside_window);
                            logpdf3d = log(mean3d(exp(logpdf.data), logical_inside_window, n_inside_window));

                            probs3d_cell = mat2cell(probs3d, n_windows, obj.n_classes, cell_converter_tool);
                            mahaldist3d_cell = mat2cell(mahaldist3d, n_windows, obj.n_classes, cell_converter_tool);
                            logpdf3d_cell = mat2cell(logpdf3d, n_windows, 1, cell_converter_tool);
                            
                            [SLAMS.region_posterior_windows] = probs3d_cell{:};
                            [SLAMS.region_mahaldist_windows] = mahaldist3d_cell{:};
                            [SLAMS.region_logpdf_windows] = logpdf3d_cell{:};
                        end
                    end
                end
                
                if r.Include_B_stats
                    t = obj.ts.b_abs.time.epoch;
                    for j = 1:n_SLAMS
                        logical_idx = (SLAMS(j).start.epoch <= t) & (SLAMS(j).stop.epoch >= t);
                        B_abs = obj.ts.b_abs.data(logical_idx);
                        B_bg = obj.ts.b_bg.data(logical_idx);
                        SLAMS(j).B_bg_mean = mean(B_bg);
                        SLAMS(j).B_mean = mean(B_abs);
                        SLAMS(j).B_max = max(B_abs);
                        B_rel = B_abs./B_bg;
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

