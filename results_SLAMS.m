
function results_SLAMS(SLAMS, dataset_info)

    SLAMS_filter = {
        'ID', [], ...
        'After_date', [], ...
        'Before_date', [], ...
        'Min_duration', [], ...
        'Max_duration', [], ...
        'GSE_filter', [], ...
        'B_filter', [], ...
        'Region_filter', []
    };

    % SLAMS_filter = setting_set(SLAMS_filter, 'B_filter', @(~,~,~,~,B_rel_max) B_rel_max > 3);
    SLAMS = filter_SLAMS(SLAMS, SLAMS_filter{:});

    n_classes = dataset_info.n_region_classes;
    color_scheme = {[1 0 0], [0 0.7 0], [0 0 1], [0 1 1]};
    ordering = [3, 4, 2, 1];
    
    % windows: [instant, 15, 30, 60, 120, 240, 480]
    window_idx = 1;
    posterior_threshold = 0.5;
    SLAMS_classes = split_SLAMS(SLAMS, posterior_threshold, window_idx);

    plot_search_time();


    plot_pos_classes(SLAMS_classes);
    
    histogram_angles(SLAMS, SLAMS_classes);
    
    plot_strength(SLAMS_classes);
    % SLAMS_filter = setting_set(SLAMS_filter, 'B_filter', @(B_bg_mean, B_mean, B_max, B_rel_mean, B_rel_max) B_rel_max > 4);
    % SLAMS_strong = filter_SLAMS(SLAMS, SLAMS_filter{:});
    % figure;
    % plot_pos(SLAMS_strong, 'k.')
    
    % figure;
    % hold on
    % plot([SLAMS_MSH.B_mean], [SLAMS_MSH.B_max], 'b.')
    % plot([SLAMS_FS.B_mean], [SLAMS_FS.B_max], 'c.')

    function plot_search_time()
        search_durations = arrayify_map(dataset_info.search_durations);
        search_dur_tot = search_durations(:,[1, 2, n_classes + 3]);
    
        search_dur_tot_h = [search_durations(:, 1:2), search_dur_tot(:, 3)/3600];
        grid_plot(search_dur_tot_h, 'use_GSE', true, 'colormap', 'hot', 'color_label', 'Search time [h]')
    
        search_dur_classes = cell(1, n_classes);
        class_ratio = cell(1, n_classes);
        color_label = {
            'Portion of time in MSP';
            'Portion of time in SW';
            'Portion of time in MSH';
            'Portion of time in FS'
        };
        for j = 1:n_classes
            search_dur_classes{j} = search_durations(:, [1, 2, j + 2]);
    
            class_ratio{j} = [search_durations(:, 1:2), search_dur_classes{j}(:, 3)./search_dur_tot(:, 3)];
    
            grid_plot(class_ratio{j}, 'use_GSE', true, 'colormap', 'hot', 'max_value', 1, 'color_label', color_label{j})
        end

        map_SLAMS_count = mapify_SLAMS(vertcat(SLAMS.pos_GSE), repelem(1, length(SLAMS))', keys(dataset_info.search_durations));
        arr_SLAMS_count = arrayify_map(map_SLAMS_count);
        % size(arr_SLAMS_count)
        % size(search_dur_tot_h)
        SLAMS_per_h = [search_durations(:, 1:2), arr_SLAMS_count(:, 3)./search_dur_tot_h(:, 3)];
        grid_plot(SLAMS_per_h, 'use_GSE', true, 'colormap', 'hot', 'color_label', 'SLAMS per hour [n/h]')
        
        map_SLAMS_strength = mapify_SLAMS(vertcat(SLAMS.pos_GSE), vertcat(SLAMS.B_rel_max), keys(dataset_info.search_durations));
        arr_SLAMS_strength = arrayify_map(map_SLAMS_strength);
        SLAMS_mean_B_max = [search_durations(:, 1:2), arr_SLAMS_strength(:, 3)./arr_SLAMS_count(:, 3)];
        grid_plot(SLAMS_mean_B_max, 'use_GSE', true, 'colormap', 'hot', 'color_label', 'SLAMS mean B_{rel,max} [nT]')

        color_label_count = {
            'MSP SLAMS per hour [n/h]';
            'SW SLAMS per hour [n/h]';
            'MSH SLAMS per hour [n/h]';
            'FS SLAMS per hour [n/h]'
        };

        color_label_strength = {
            'MSP SLAMS mean B_{rel,max} [nT]';
            'SW SLAMS mean B_{rel,max} [nT]';
            'MSH SLAMS mean B_{rel,max} [nT]';
            'FS SLAMS mean B_{rel,max} [nT]'
        };
        for j = [3,4]
            map_SLAMS_count = mapify_SLAMS(vertcat(SLAMS_classes{j}.pos_GSE), repelem(1, length(SLAMS_classes{j}))', keys(dataset_info.search_durations));
            arr_SLAMS_count = arrayify_map(map_SLAMS_count);
            SLAMS_per_h = [search_durations(:, 1:2), arr_SLAMS_count(:, 3)./search_dur_tot_h(:, 3)];
            grid_plot(SLAMS_per_h, 'use_GSE', false, 'colormap', 'hot', 'color_label', color_label_count{j})

            map_SLAMS_strength = mapify_SLAMS(vertcat(SLAMS_classes{j}.pos_GSE), vertcat(SLAMS_classes{j}.B_rel_max), keys(dataset_info.search_durations));
            arr_SLAMS_strength = arrayify_map(map_SLAMS_strength);
            SLAMS_mean_B_max = [search_durations(:, 1:2), arr_SLAMS_strength(:, 3)./arr_SLAMS_count(:, 3)];
            grid_plot(SLAMS_mean_B_max, 'use_GSE', false, 'colormap', 'hot', 'color_label', color_label_strength{j})
        end
    end

    function cm = create_colormap(minColor, maxColor)
        cm = interp1([0; 1], [minColor; maxColor], linspace(0,1,256));
    end
    
    function plot_pos_classes(SLAMS_classes)
        figure;
        hold on
        % plot_pos(SLAMS, [0.5 0.5 0.5])
        for i = ordering
            plot_pos(SLAMS_classes{i}, color_scheme{i})
        end
        plot_reference(true, 'k')
        ylabel('\rho R_{GSE}')
        xlabel('x R_{GSE}')
        axis equal
        grid on
        legend(dataset_info.region_classes{ordering})

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
        legend('ALL', dataset_info.region_classes{ordering})
        
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

        for i = 1:n_classes
            SLAMS_classes{i} = filter_SLAMS(SLAMS, 'Region_filter', @(p,m,l) probable_class(p, i));
            % SLAMS_classes{i} = filter_SLAMS(SLAMS, 'Region_filter', @(p,m,l) p(:,i,window_idx) > posterior_threshold);
        end

        function bool = probable_class(p, clas)
            [~, idx] = max(p(:,:,window_idx), [], 2);
            bool = idx == clas;
        end
    end

    function grid_plot(gridd, varargin)

        p = inputParser;
        addRequired(p, 'gridd');
        addParameter(p, 'use_GSE', false)
        addParameter(p, 'colormap', 'hot')
        addParameter(p, 'color_label', [])
        addParameter(p, 'max_value', [])
        parse(p, gridd, varargin{:})
        r = p.Results;

        [n_gridd, col] = size(gridd);
        if isempty(r.max_value)
            r.max_value = max(gridd(:,3));
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
        caxis([0, r.max_value]);
        if ~isempty(r.color_label)
            cb.Label.String = r.color_label;
        end
        axis equal

        plot_reference(r.use_GSE, [0.2 0.2 1])
        
        function c = interpolate_color(value)
            row_idx = ceil(256*min(value/r.max_value, 1));
            if row_idx == 0
                row_idx = 1;
            end
            c = cm(row_idx, :);
        end

        function plot_patch(x, y, c)
            down = y - 1;
            up = y;
            if r.use_GSE
                y_int = linspace(down, up, 10)';
                xBS = bowshock_pos(y_int);
                X = [xBS + x - 1; flip(xBS) + x];
                Y = [y_int; flip(y_int)];
            else
                left = x - 1;
                right = x;
                X = [left left right right];
                Y = [down up up down];
            end
            patch(X, Y, c);
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
end