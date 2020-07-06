
function SLAMS_results(SLAMS_database)

    gridsize = 1;

    SLAMS_primary = SLAMS_database.SLAMS_primary;
    % SLAMS_primary = filter_short(SLAMS_primary, 1);
    % SLAMS_primary = filter_weak(SLAMS_primary, 2.5);
    search_durations_primary = SLAMS_database.search_durations_primary;
    
    SLAMS_secondary = SLAMS_database.SLAMS_secondary;
    % SLAMS_secondary = filter_short(SLAMS_secondary, 1);
    % SLAMS_secondary = filter_weak(SLAMS_secondary, 2.5);
    search_durations_secondary = SLAMS_database.search_durations_secondary;

    SLAMS_combined = combine_SLAMS(SLAMS_primary, SLAMS_secondary);
    search_durations_combined = combine_map(search_durations_primary, 5, search_durations_secondary, 1);

    n_classes = SLAMS_database.n_region_classes;
    color_scheme = {[1 0 0], [0 0.7 0], [0 0 1], [0 1 1]};
    ordering = [3, 4, 2, 1];
    
    % windows: [instant, 15, 30, 60, 120, 240, 480]
    window_idx = 3;
    posterior_threshold = 0.5;
    SLAMS_classes = split_SLAMS(SLAMS_primary, posterior_threshold, window_idx);

    disp(['n_SLAMS_primary = ', num2str(length(SLAMS_primary))])
    disp(['n_SLAMS_secondary = ', num2str(length(SLAMS_secondary))])
    disp(['n_SLAMS_tot = ', num2str(length(SLAMS_primary) + length(SLAMS_secondary))])
    disp(['n_SLAMS_MSP = ', num2str(length(SLAMS_classes{1}))])
    disp(['n_SLAMS_SW = ', num2str(length(SLAMS_classes{2}))])
    disp(['n_SLAMS_MSH = ', num2str(length(SLAMS_classes{3}))])
    disp(['n_SLAMS_FS = ', num2str(length(SLAMS_classes{4}))])

    % Plot position of classified and secondary SLAMS
    figure;
    hold on
    grid on
    axis equal
    plot_pos(SLAMS_primary, [0, 0.4470, 0.7410])
    plot_pos(SLAMS_secondary, [0.8500, 0.3250, 0.0980])
    plot_reference(true, 'k')
    legend('fpi active', 'fpi inactive')
    title('Detected SLAMS')
    
    % Plot position of classified SLAMS
    figure;
    hold on
    grid on
    axis equal
    plot_pos_classes(SLAMS_classes);
    plot_reference(true, 'k')
    legend(SLAMS_database.region_classes{ordering})
    title('Detected SLAMS')
    
    % Plot search durations when fpi is active
    gridplot_settings = {'use_GSE', true, 'colormap', 'hot', 'color_label', 'Hours [h]', 'ref_color', 'c', 'gridsize', gridsize};
    plot_search_durations('Search time: fpi active', search_durations_primary, 5, gridplot_settings, false)
    
    % % Plot search durations when fpi is inactive
    % gridplot_settings = {'use_GSE', true, 'colormap', 'hot', 'color_label', 'Hours [h]', 'ref_color', 'c', 'gridsize', gridsize};
    % plot_search_durations('Search time: fpi inactive', search_durations_secondary, 1, gridplot_settings, false)
    
    % Plot search durations when fgm is active
    gridplot_settings = {'use_GSE', true, 'colormap', 'hot', 'color_label', 'Hours [h]', 'ref_color', 'c', 'gridsize', gridsize};
    plot_search_durations('Search time: fgm active', search_durations_combined, 1, gridplot_settings, false)

    % Plot portion of time spent in each region
    for j = 1:n_classes
        gridplot_settings = {'use_GSE', true, 'colormap', 'hot', 'max_value', 100, 'color_label', 'Percent [%]', 'ref_color', 'c', 'gridsize', gridsize};
        plot_search_durations(['Portion of time in ', SLAMS_database.region_classes{j}], search_durations_primary, [j, 5], gridplot_settings, true)
    end




    % Plot SLAMS rate
    gridplot_settings = {'use_GSE', false, 'colormap', 'hot', 'max_value', 50, 'color_label', 'Count per hour [n/h]', 'ref_color', 'c', 'gridsize', gridsize};
    plot_SLAMS_rate('SLAMS rate', SLAMS_combined, search_durations_combined, 1, gridplot_settings)
    
    % Plot SLAMS rel strength
    gridplot_settings = {'use_GSE', false, 'colormap', 'jet', 'min_value', 2, 'max_value', 3.5, 'color_label', 'Mean B_{max}/B_0 [-]', 'ref_color', 'r', 'gridsize', gridsize};
    plot_SLAMS_strength('SLAMS relative strength', SLAMS_combined, gridplot_settings, true)

    % Plot SLAMS strength
    gridplot_settings = {'use_GSE', false, 'colormap', 'jet', 'min_value', 0, 'max_value', 50, 'color_label', 'Mean B_{max} [nT]', 'ref_color', 'r', 'gridsize', gridsize};
    plot_SLAMS_strength('SLAMS strength', SLAMS_combined, gridplot_settings, false)

    rate_specific = [nan, nan, 100, 50];
    rate_global = [nan, nan, 50, 10];

    for j = [3,4]
        % Plot SLAMS classes rate
        gridplot_settings = {'use_GSE', false, 'colormap', 'hot', 'max_value', rate_global(j), 'color_label', 'Count per hour [n/h]', 'ref_color', 'c', 'gridsize', gridsize};
        plot_SLAMS_rate([SLAMS_database.region_classes{j}, ' SLAMS rate'], SLAMS_classes{j}, search_durations_primary, 5, gridplot_settings)

        gridplot_settings = {'use_GSE', false, 'colormap', 'hot', 'max_value', rate_specific(j), 'color_label', ['Count per ', SLAMS_database.region_classes{j}, ' hour [n/h]'], 'ref_color', 'c', 'gridsize', gridsize};
        plot_SLAMS_rate([SLAMS_database.region_classes{j}, ' SLAMS rate in ' SLAMS_database.region_classes{j}], SLAMS_classes{j}, search_durations_primary, j, gridplot_settings)
        
        % Plot SLAMS classes rel strength
        gridplot_settings = {'use_GSE', false, 'colormap', 'jet', 'min_value', 2, 'max_value', 3.5, 'color_label', 'Mean B_{max}/B_0 [-]', 'ref_color', 'r', 'gridsize', gridsize};
        plot_SLAMS_strength([SLAMS_database.region_classes{j}, ' SLAMS relative strength'], SLAMS_classes{j}, gridplot_settings, true)

        % Plot SLAMS classes strength
        gridplot_settings = {'use_GSE', false, 'colormap', 'jet', 'min_value', 0, 'max_value', 50, 'color_label', 'Mean B_{max} [-]', 'ref_color', 'r', 'gridsize', gridsize};
        plot_SLAMS_strength([SLAMS_database.region_classes{j}, ' SLAMS strength'], SLAMS_classes{j}, gridplot_settings, false)
    end
    
    % plot_strength(SLAMS_classes);
    
    function plot_pos_classes(SLAMS_classes)
        for i = ordering
            plot_pos(SLAMS_classes{i}, color_scheme{i})
        end
    end

    function plot_reference(use_GSE, c)
        if use_GSE
            K = 25;
            epsilon = 0.8;
            a = linspace(0, 2*pi);
            r = K./(1 + epsilon*cos(a));
            x = cos(a);
            y = sin(a);
            plot(x, y, 'color', c, 'LineWidth', 1.5)
            plot(x.*r, y.*r, 'color', c, 'LineWidth', 1.5)
            xlim([-22, 30])
            ylim([0, 35])
            xlabel('x R_E')
        else
            plot([0, 0], [0,35], 'color', c, 'LineWidth', 1.5)
            xlim([-30, 18])
            ylim([0, 35])
            xlabel('x - x_{BS} R_E')
        end
        ylabel('\rho R_E')
        % plot(x.*r + 11.5, y.*r, 'k')
        % plot(x.*r + 15, y.*r, 'k')

        % pos = [repelem(30, 100)', zeros(100,1), (y.*r)'];
        % d = dist_to_bow_shock(pos);
        % plot(x.*r + d', y.*r, 'k')
    end
    
    function plot_pos(SLAMS, c)
        pos = vertcat(SLAMS.pos_GSE)/6378;
        pos_yz = irf_abs(pos(:, 2:3), 1);
        scatter(pos(:, 1), pos_yz, 36, c, '.')
    end
    
    % function histogram_angles(SLAMS, SLAMS_classes)
    %     edges = linspace(0, 180, 19);

    %     figure;
    %     hold on
    %     histogram(calc_angle(SLAMS), edges, 'FaceColor', [0.5 0.5 0.5]);
    %     for i = ordering
    %         histogram(calc_angle(SLAMS_classes{i}), edges, 'FaceColor', color_scheme{i});
    %     end
    %     ylabel('Count SLAMS')
    %     xlabel('Angle from [1 0 0] GSE')
    %     grid on
    %     legend('ALL', SLAMS_database.region_classes{ordering})
        
    %     function ang = calc_angle(SLAMS)
    %         pos_GSE = vertcat(SLAMS.pos_GSE);
    %         ang = acosd(pos_GSE(:,1)./(sqrt(sum(pos_GSE.^2, 2))));
    %     end
    % end
    
    function plot_strength(SLAMS_classes, plot_3d)
        if nargin == 1
            plot_3d = false;
        end
        figure;
        hold on
        for i = ordering
            pos = vertcat(SLAMS_classes{i}.pos_GSE)/6378;
            yz = irf_abs(pos(:,2:3), 1);
            xBS = bowshock_pos(yz);
            x = pos(:,1) - xBS;
            s = vertcat(SLAMS_classes{i}.B_max);
            if plot_3d
                scatter3(x, yz, s, 36, color_scheme{i}, '.')
            else
                scatter(x, s, 36, color_scheme{i}, '.')
            end
        end
    end
    
    function SLAMS_classes = split_SLAMS(SLAMS, posterior_threshold, window_idx)
        
        SLAMS_classes = cell(1, n_classes);

        posteriors = cat(3, SLAMS.region_posterior);
        posteriors = permute(posteriors(window_idx,:,:), [3, 2, 1]);
        [~, idx] = max(posteriors, [], 2);

        for i = 1:n_classes
            SLAMS_classes{i} = SLAMS(idx == i);
            % SLAMS_classes{i} = filter_SLAMS(SLAMS, 'Region_filter', @(p,m,l) probable_class(p, i));
            % % SLAMS_classes{i} = filter_SLAMS(SLAMS, 'Region_filter', @(p,m,l) p(:,i,window_idx) > posterior_threshold);
        end

        % function bool = probable_class(p, clas)
        %     [~, idx] = max(p(:,:,window_idx), [], 2);
        %     bool = idx == clas;
        % end
    end

    function grid_plot(gridd, varargin)

        p = inputParser;
        addRequired(p, 'gridd');
        addParameter(p, 'use_GSE', false)
        addParameter(p, 'colormap', 'hot')
        addParameter(p, 'color_label', [])
        addParameter(p, 'max_value', [])
        addParameter(p, 'min_value', 0)
        addParameter(p, 'ref_color', 'k')
        addParameter(p, 'gridsize', 1)
        parse(p, gridd, varargin{:})
        r = p.Results;

        [n_gridd, col] = size(gridd);
        if isempty(r.max_value)
            r.max_value = max(gridd(:,3));
        end
        if isempty(r.min_value)
            r.min_value = min(gridd(:,3));
        end
        figure;
        hold on
        if ischar(r.colormap)
            cm = colormap(r.colormap);
        else
            cm = r.colormap;
        end

        for i = 1:n_gridd
            if ~isnan(gridd(i,3))
                c = interpolate_color(gridd(i,3));
                plot_patch(gridd(i,1), gridd(i,2), c);
            end
        end

        cb = colorbar();
        caxis([r.min_value, r.max_value]);
        if ~isempty(r.color_label)
            cb.Label.String = r.color_label;
        end
        axis equal
        grid on

        plot_reference(r.use_GSE, r.ref_color)
        
        function c = interpolate_color(value)
            row_idx = ceil(256*min((value - r.min_value)/(r.max_value - r.min_value), 1));
            if row_idx <= 0
                row_idx = 1;
            end
            c = cm(row_idx, :);
        end

        function plot_patch(x, y, c)
            down = y - r.gridsize;
            up = y;
            if r.use_GSE
                y_int = linspace(down, up, 10*r.gridsize)';
                xBS = bowshock_pos(y_int);
                X = [xBS + x - r.gridsize; flip(xBS) + x];
                Y = [y_int; flip(y_int)];
            else
                left = x - r.gridsize;
                right = x;
                X = [left left right right];
                Y = [down up up down];
            end
            patch(X, Y, c, 'LineStyle', '-');
        end
    end

    function map = mapify_SLAMS(pos_GSE, property, keys)
        map = containers.Map;
        if nargin == 3
            for i = keys
                map(i{:}) = 0;
            end
        end

        pos = pos_GSE/6378;
        BS = GSE2BS(pos);

        xBS_disc = ceil(BS(:,1));
        yBS_disc = ceil(BS(:,2));

        for i = 1:length(xBS_disc)
            if ~isnan(xBS_disc(i)) && ~isnan(yBS_disc(i))
                key = [num2str(xBS_disc(i)), ',', num2str(yBS_disc(i))];
                if isKey(map, key)
                    map(key) = map(key) + property(i);
                else
                    map(key) = property(i);
                end
            end
        end
    end

    function rescaled_map = rescale_map(map, gridsize)
        if gridsize == 1
            rescaled_map = map;
            return
        end
        rescaled_map = containers.Map;
        k = keys(map);
        for i = k
            coords = sscanf(i{:}, '%d,%d');
            rescaled_coords = ceil(coords/gridsize)*gridsize;
            key = sprintf('%d,%d', rescaled_coords(:));
            if isKey(rescaled_map, key)
                rescaled_map(key) = rescaled_map(key) + map(i{:});
            else
                rescaled_map(key) = map(i{:});
            end
        end
    end

    function array = arrayify_map(map)
        k = keys(map);
        n_keys = length(k);
        array = zeros(n_keys, 2);
        for i = 1:n_keys
            coords = sscanf(k{i}, '%d,%d');
            array(i,1:2) = coords;
            value = map(k{i});
            n_value = length(value);
            array(i,3:(3 + n_value - 1)) = value;
        end
    end

    function map = combine_map(map1, idx1, map2, idx2)
        map = containers.Map;
        k = keys(map1);
        for i = k
            v = map1(i{:});
            map(i{:}) = v(idx1);
        end
        k = keys(map2);
        for i = k
            v = map2(i{:});
            if isKey(map, i{:})
                map(i{:}) = map(i{:}) + v(idx2);
            else
                map(i{:}) = v(idx2);
            end
        end
    end

    function flat_arr = flatten_array(arr)
        uniq = unique(arr(:,1));
        l = length(uniq);
        flat_arr = [uniq, zeros(l,1)];
        for i = 1:l
            % row_idx = find(arr(:,1) == uniq(i));
            flat_arr(i, 2) = sum(arr(arr(:,1) == uniq(i), 3));
        end
    end

    function map = map_zero_to_nan(map)
        k = keys(map);
        for i = k
            if map(i{:}) == 0
                map(i{:}) = nan;
            end
        end
    end

    function SLAMS_combined = combine_SLAMS(SLAMS, SLAMS_secondary)
        SLAMS = rmfield(SLAMS, 'region_posterior');
        SLAMS = rmfield(SLAMS, 'region_mahaldist');
        SLAMS = rmfield(SLAMS, 'region_logpdf');
        SLAMS_combined = vertcat(SLAMS, SLAMS_secondary);
    end

    function plot_search_durations(name, search_dur_map, k, gridplot_settings, relative)
        search_dur = arrayify_map(rescale_map(search_dur_map, setting_get(gridplot_settings, 'gridsize')));
        search_dur = search_dur(:, [1, 2, 2 + k]);
        if relative
            search_dur = [search_dur(:, 1:2), search_dur(:, 3)./search_dur(:, 4)*100];
        else
            search_dur = [search_dur(:, 1:2), search_dur(:, 3)/3600];
        end
        grid_plot(search_dur, gridplot_settings{:})
        title(name)
    end

    function plot_SLAMS_rate(name, SLAMS, search_dur, k, gridplot_settings)
        map_SLAMS_count = mapify_SLAMS(vertcat(SLAMS.pos_GSE), repelem(1, length(SLAMS))', keys(search_dur));
        gridsize = setting_get(gridplot_settings, 'gridsize');
        arr_SLAMS_count = arrayify_map(map_zero_to_nan(rescale_map(map_SLAMS_count, gridsize)));
        arr_search_dur = arrayify_map(rescale_map(search_dur, gridsize));
        assert(all(arr_SLAMS_count(:,1) == arr_search_dur(:,1)));
        SLAMS_per_h = [arr_search_dur(:, 1:2), 3600*arr_SLAMS_count(:, 3)./arr_search_dur(:, 2 + k)];
        grid_plot(SLAMS_per_h, gridplot_settings{:})
        title(name)
    end

    function plot_SLAMS_strength(name, SLAMS, gridplot_settings, relative)
        pos = vertcat(SLAMS.pos_GSE);
        gridsize = setting_get(gridplot_settings, 'gridsize');
        map_SLAMS_count = mapify_SLAMS(pos, repelem(1, length(SLAMS))');
        arr_SLAMS_count = arrayify_map(rescale_map(map_SLAMS_count, gridsize));
        if relative
            map_SLAMS_strength = mapify_SLAMS(pos, vertcat(SLAMS.B_rel_max));
        else
            map_SLAMS_strength = mapify_SLAMS(pos, vertcat(SLAMS.B_max));
        end
        arr_SLAMS_strength = arrayify_map(rescale_map(map_SLAMS_strength, gridsize));
        arr_SLAMS_mean_strength = [arr_SLAMS_strength(:, 1:2), arr_SLAMS_strength(:, 3)./arr_SLAMS_count(:, 3)];
        grid_plot(arr_SLAMS_mean_strength, gridplot_settings{:});
        title(name);
    end

    % function plot_SLAMS_monolithness(name, SLAMS, gridplot_settings)
    %     pos = vertcat(SLAMS.pos_GSE);
    %     gridsize = setting_get(gridplot_settings, 'gridsize');
    %     map_SLAMS_monolithness = mapify_SLAMS(pos, vertcat(SLAMS.stop) - vertcat(SLAMS.start));
    %     arr_SLAMS_monolithness = arrayify_map(rescale_map(map_SLAMS_monolithness, gridsize));
    %     map_SLAMS_count = mapify_SLAMS(pos, repelem(1, length(SLAMS))');
    %     arr_SLAMS_count = arrayify_map(rescale_map(map_SLAMS_count, gridsize));
    %     arr_SLAMS_mean_monolithness = [arr_SLAMS_monolithness(:, 1:2), arr_SLAMS_monolithness(:, 3)./arr_SLAMS_count(:, 3)];
    %     grid_plot(arr_SLAMS_mean_monolithness, gridplot_settings{:});
    %     title(name);
    % end

    function SLAMS = filter_short(SLAMS, lim)
        dur = [SLAMS.stop] - [SLAMS.start];
        SLAMS = SLAMS(dur > lim);
    end

    function SLAMS = filter_weak(SLAMS, lim)
        strength = [SLAMS.B_max];
        SLAMS = SLAMS(strength > lim);
    end    
end