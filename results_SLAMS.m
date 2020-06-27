
function results_SLAMS(SLAMS_database)

    % SLAMS_filter = {
    %     'ID', [], ...
    %     'After_date', [], ...
    %     'Before_date', [], ...
    %     'Min_duration', [], ...
    %     'Max_duration', [], ...
    %     'GSE_filter', [], ...
    %     'B_filter', [], ...
    %     'Region_filter', []
    % };

    % % SLAMS_filter = setting_set(SLAMS_filter, 'B_filter', @(~,~,~,~,B_rel_max) B_rel_max > 3);
    % SLAMS = filter_SLAMS(SLAMS_database.SLAMS, SLAMS_filter{:});

    gridsize = 1;

    SLAMS = SLAMS_database.SLAMS;
    search_durations = SLAMS_database.search_durations;
    SLAMS_unclassified = SLAMS_database.SLAMS_unclassified;
    search_durations_unclassified = SLAMS_database.search_durations_unclassified;
    search_durations_combined = combine_map(search_durations, 5, search_durations_unclassified, 1);

    n_classes = SLAMS_database.n_region_classes;
    color_scheme = {[1 0 0], [0 0.7 0], [0 0 1], [0 1 1]};
    ordering = [3, 4, 2, 1];
    
    % windows: [instant, 15, 30, 60, 120, 240, 480]
    window_idx = 3;
    posterior_threshold = 0.5;
    SLAMS_classes = split_SLAMS(SLAMS, posterior_threshold, window_idx);

    disp(['n_SLAMS = ', num2str(length(SLAMS))])
    disp(['n_SLAMS_unclassified = ', num2str(length(SLAMS_unclassified))])
    disp(['n_SLAMS_tot = ', num2str(length(SLAMS) + length(SLAMS_unclassified))])
    disp(['n_SLAMS_MSP = ', num2str(length(SLAMS_classes{1}))])
    disp(['n_SLAMS_SW = ', num2str(length(SLAMS_classes{2}))])
    disp(['n_SLAMS_MSH = ', num2str(length(SLAMS_classes{3}))])
    disp(['n_SLAMS_FS = ', num2str(length(SLAMS_classes{4}))])

    % Plot position of classified and unclassified SLAMS
    figure;
    hold on
    grid on
    axis equal
    plot_pos(SLAMS, [0, 0.4470, 0.7410])
    plot_pos(SLAMS_unclassified, [0.8500, 0.3250, 0.0980])
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
    search_durations = arrayify_map(rescale_map(SLAMS_database.search_durations, gridsize));
    search_dur_tot = search_durations(:,[1, 2, n_classes + 3]);
    search_dur_tot_h = [search_durations(:, 1:2), search_dur_tot(:, 3)/3600];
    grid_plot(search_dur_tot_h, 'use_GSE', true, 'colormap', 'hot', 'color_label', 'Hours [h]', 'ref_color', 'c', 'gridsize', gridsize)
    title('Search time: fpi active')
    
    % % Plot search durations when fpi is inactive
    % search_durations_unclassified = arrayify_map(SLAMS_database.search_durations_unclassified);
    % % search_dur_tot = search_durations_unclassified(:,[1, 2, 3]);
    % search_dur_tot_h = [search_durations_unclassified(:, 1:2), search_durations_unclassified(:, 3)/3600];
    % grid_plot(search_dur_tot_h, 'use_GSE', true, 'colormap', 'hot', 'color_label', 'Hours [h]')
    % title('Search time: fpi inactive')
    
    % Plot search durations when fgm is active
    search_durations_combined = arrayify_map(rescale_map(search_durations_combined, gridsize));
    search_durations_combined_h = [search_durations_combined(:,1:2), search_durations_combined(:,3)/3600];
    grid_plot(search_durations_combined_h, 'use_GSE', true, 'colormap', 'hot', 'color_label', 'Hours [h]', 'ref_color', 'c', 'gridsize', gridsize)
    title('Search time: fgm active')
    
    % Plot portion of time spent in each region
    search_dur_classes = cell(1, n_classes);
    class_ratio = cell(1, n_classes);
    label = {
        'Portion of time in MSP';
        'Portion of time in SW';
        'Portion of time in MSH';
        'Portion of time in FS'
    };
    for j = 1:n_classes
        search_dur_classes{j} = search_durations(:, [1, 2, j + 2]);

        class_ratio{j} = [search_durations(:, 1:2), search_dur_classes{j}(:, 3)./search_dur_tot(:, 3)*100];

        grid_plot(class_ratio{j}, 'use_GSE', true, 'colormap', 'hot', 'max_value', 100, 'color_label', 'Percent [%]', 'ref_color', 'c', 'gridsize', gridsize)
        title(label{j})
    end

    % Plot SLAMS rate
    map_SLAMS_count = mapify_SLAMS(vertcat(SLAMS.pos_GSE), repelem(1, length(SLAMS))', keys(SLAMS_database.search_durations));
    map_SLAMS_unclassified_count = mapify_SLAMS(vertcat(SLAMS_unclassified.pos_GSE), repelem(1, length(SLAMS_unclassified))', keys(SLAMS_database.search_durations_unclassified));
    map_SLAMS_combined_count = combine_map(map_SLAMS_count, 1, map_SLAMS_unclassified_count, 1);
    arr_SLAMS_combined_count = arrayify_map(map_zero_to_nan(rescale_map(map_SLAMS_combined_count, gridsize)));
    SLAMS_per_h = [search_durations_combined(:, 1:2), arr_SLAMS_combined_count(:, 3)./search_durations_combined_h(:, 3)];
    grid_plot(SLAMS_per_h, 'use_GSE', false, 'colormap', 'hot', 'color_label', 'Count per hour [n/h]', 'ref_color', 'r', 'gridsize', gridsize)
    title('SLAMS rate')
    
    % Plot SLAMS rel strength
    map_SLAMS_strength = mapify_SLAMS(vertcat(SLAMS.pos_GSE), vertcat(SLAMS.B_rel_max), keys(SLAMS_database.search_durations));
    map_SLAMS_unclassified_strength = mapify_SLAMS(vertcat(SLAMS_unclassified.pos_GSE), vertcat(SLAMS_unclassified.B_rel_max), keys(SLAMS_database.search_durations_unclassified));
    map_SLAMS_combined_strength = combine_map(map_SLAMS_strength, 1, map_SLAMS_unclassified_strength, 1);
    arr_SLAMS_combined_strength = arrayify_map(rescale_map(map_SLAMS_combined_strength, gridsize));
    SLAMS_mean_B_max = [search_durations_combined(:, 1:2), arr_SLAMS_combined_strength(:, 3)./arr_SLAMS_combined_count(:, 3)];
    grid_plot(SLAMS_mean_B_max, 'use_GSE', false, 'colormap', 'jet', 'min_value', 2, 'max_value', 3, 'color_label', 'Mean B_{max}/B_0 [-]', 'ref_color', 'r', 'gridsize', gridsize)
    title('SLAMS relative strength')

    % Plot SLAMS strength
    map_SLAMS_strength = mapify_SLAMS(vertcat(SLAMS.pos_GSE), vertcat(SLAMS.B_max), keys(SLAMS_database.search_durations));
    map_SLAMS_unclassified_strength = mapify_SLAMS(vertcat(SLAMS_unclassified.pos_GSE), vertcat(SLAMS_unclassified.B_max), keys(SLAMS_database.search_durations_unclassified));
    map_SLAMS_combined_strength = combine_map(map_SLAMS_strength, 1, map_SLAMS_unclassified_strength, 1);
    arr_SLAMS_combined_strength = arrayify_map(rescale_map(map_SLAMS_combined_strength, gridsize));
    SLAMS_mean_B_max = [search_durations_combined(:, 1:2), arr_SLAMS_combined_strength(:, 3)./arr_SLAMS_combined_count(:, 3)];
    grid_plot(SLAMS_mean_B_max, 'use_GSE', false, 'colormap', 'jet', 'min_value', 0, 'max_value', 50, 'color_label', 'Mean B_{max} [nT]', 'ref_color', 'r', 'gridsize', gridsize)
    title('SLAMS strength')
    
    label_count = {
        'MSP SLAMS rate';
        'SW SLAMS rate';
        'MSH SLAMS rate';
        'FS SLAMS rate'
    };

    label_rel_strength = {
        'MSP SLAMS relative strength';
        'SW SLAMS relative strength';
        'MSH SLAMS relative strength';
        'FS SLAMS relative strength'
    };
        
    label_strength = {
        'MSP SLAMS strength';
        'SW SLAMS strength';
        'MSH SLAMS strength';
        'FS SLAMS strength'
    };
    for j = [3,4]
        % Plot SLAMS classes rate
        map_SLAMS_count = mapify_SLAMS(vertcat(SLAMS_classes{j}.pos_GSE), repelem(1, length(SLAMS_classes{j}))', keys(SLAMS_database.search_durations));
        arr_SLAMS_count = arrayify_map(map_zero_to_nan(rescale_map(map_SLAMS_count, gridsize)));
        SLAMS_per_h = [search_durations(:, 1:2), arr_SLAMS_count(:, 3)./search_dur_tot_h(:, 3)];
        grid_plot(SLAMS_per_h, 'use_GSE', false, 'colormap', 'hot', 'color_label', 'Count per hour [n/h]', 'ref_color', 'r', 'gridsize', gridsize)
        title(label_count{j});
        
        % Plot SLAMS classes rel strength
        map_SLAMS_strength = mapify_SLAMS(vertcat(SLAMS_classes{j}.pos_GSE), vertcat(SLAMS_classes{j}.B_rel_max), keys(SLAMS_database.search_durations));
        arr_SLAMS_strength = arrayify_map(rescale_map(map_SLAMS_strength, gridsize));
        SLAMS_mean_B_max = [search_durations(:, 1:2), arr_SLAMS_strength(:, 3)./arr_SLAMS_count(:, 3)];
        grid_plot(SLAMS_mean_B_max, 'use_GSE', false, 'colormap', 'jet', 'min_value', 2, 'max_value', 3, 'color_label', 'Mean B_{max}/B_0 [-]', 'ref_color', 'r', 'gridsize', gridsize)
        title(label_rel_strength{j});

        % Plot SLAMS classes strength
        map_SLAMS_strength = mapify_SLAMS(vertcat(SLAMS_classes{j}.pos_GSE), vertcat(SLAMS_classes{j}.B_max), keys(SLAMS_database.search_durations));
        arr_SLAMS_strength = arrayify_map(rescale_map(map_SLAMS_strength, gridsize));
        SLAMS_mean_B_max = [search_durations(:, 1:2), arr_SLAMS_strength(:, 3)./arr_SLAMS_count(:, 3)];
        grid_plot(SLAMS_mean_B_max, 'use_GSE', false, 'colormap', 'jet', 'min_value', 0, 'max_value', 50, 'color_label', 'Mean B_{max} [-]', 'ref_color', 'r', 'gridsize', gridsize)
        title(label_strength{j});
    end
    
    % histogram_angles(SLAMS, SLAMS_classes);
    
    plot_strength(SLAMS_classes);
    % SLAMS_filter = setting_set(SLAMS_filter, 'B_filter', @(B_bg_mean, B_mean, B_max, B_rel_mean, B_rel_max) B_rel_max > 4);
    % SLAMS_strong = filter_SLAMS(SLAMS, SLAMS_filter{:});
    % figure;
    % plot_pos(SLAMS_strong, 'k.')
    
    % figure;
    % hold on
    % plot([SLAMS_MSH.B_mean], [SLAMS_MSH.B_max], 'b.')
    % plot([SLAMS_FS.B_mean], [SLAMS_FS.B_max], 'c.')

    % function cm = create_colormap(minColor, maxColor)
    %     cm = interp1([0; 1], [minColor; maxColor], linspace(0,1,256));
    % end
    
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
            xlabel('x R_{GSE}')
        else
            plot([0, 0], [0,35], 'color', c, 'LineWidth', 1.5)
            xlim([-30, 18])
            ylim([0, 35])
            xlabel('x - x_{BS} R_{GSE}')
        end
        ylabel('\rho R_{GSE}')
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
    
    function histogram_angles(SLAMS, SLAMS_classes)
        edges = linspace(0, 180, 19);

        figure;
        hold on
        histogram(calc_angle(SLAMS), edges, 'FaceColor', [0.5 0.5 0.5]);
        for i = ordering
            histogram(calc_angle(SLAMS_classes{i}), edges, 'FaceColor', color_scheme{i});
        end
        ylabel('Count SLAMS')
        xlabel('Angle from [1 0 0] GSE')
        grid on
        legend('ALL', SLAMS_database.region_classes{ordering})
        
        function ang = calc_angle(SLAMS)
            pos_GSE = vertcat(SLAMS.pos_GSE);
            ang = acosd(pos_GSE(:,1)./(sqrt(sum(pos_GSE.^2, 2))));
        end
    end
    
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
            s = vertcat(SLAMS_classes{i}.B_rel_max);
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
end