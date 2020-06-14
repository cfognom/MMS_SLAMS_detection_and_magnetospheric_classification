
function results_SLAMS(SLAMS, dataset_info)

    SLAMS_filter = {
        'ID', [], ...
        'After_date', '2017-09-08T00:00:00.000000000Z', ...
        'Before_date', [], ...
        'Min_duration', [], ...
        'Max_duration', [], ...
        'GSE_filter', [], ...
        'B_filter', [], ...
        'Region_filter', []
    };

    SLAMS = filter_SLAMS(SLAMS, SLAMS_filter{:});

    n_classes = dataset_info.n_region_classes;
    color_scheme = {[1 0 0], [0 0.7 0], [0 0 1], [0 1 1]};
    ordering = [3, 4, 2, 1];

    % windows: [instant, 15, 30, 60, 120, 240, 480]
    window_idx = 4;
    posterior_threshold = 0.5;
    SLAMS_classes = split_SLAMS(SLAMS, posterior_threshold, window_idx);
    
    plot_search_durations(dataset_info.search_durations, true)
    df

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
    
    
    function plot_pos_classes(SLAMS_classes)
        figure;
        hold on
        % plot_pos(SLAMS, [0.5 0.5 0.5])
        for i = ordering
            plot_pos(SLAMS_classes{i}, color_scheme{i})
        end
        plot_reference('k')
        ylabel('\rho R_{GSE}')
        xlabel('x R_{GSE}')
        axis equal
        grid on
        legend(dataset_info.region_classes{ordering})

    end

    function plot_reference(c)
        K = 25;
        epsilon = 0.8;
        a = linspace(0, 2*pi);
        r = K./(1 + epsilon*cos(a));
        x = cos(a);
        y = sin(a);
        plot(x, y, c)
        plot(x.*r, y.*r, c)
        xlim([-22, 30])
        ylim([0, 35])
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
            % SLAMS_classes{i} = filter_SLAMS(SLAMS, 'Region_filter', @(p,m,l) probable_class(p, i));
            SLAMS_classes{i} = filter_SLAMS(SLAMS, 'Region_filter', @(p,m,l) p(:,i,window_idx) > posterior_threshold);
        end

        function bool = probable_class(p, clas)
            [~, idx] = max(p(:,:,window_idx), [], 2);
            bool = idx == clas;
        end
    end

    function plot_search_durations(map, use_BS_coords)
        if nargin == 1
            use_BS_coords = false;
        end
        figure;
        hold on
        k = keys(map);
        v = values(map);
        v = vertcat(v{:});
        m = max(v(:,end));
        max_h = m/(60*60);
        sum(v(:,end))/(60*60)
        cm = jet;
        for i = k
            arr_key = sscanf(i{:}, '%d,%d');
            v = map(i{:});
            c = cm(round(1 + 255*v(end)/m), :);
            plot_rec(arr_key(1), arr_key(2), c);
        end
        colorbar
        set(gca, 'ColorScale', 'log')
        
        plot_reference('r')

        function plot_rec(x, y, c)
            down = y - 1;
            up = y;
            if ~use_BS_coords
                left = x - 1;
                right = x;
                X = [left left right right];
                Y = [down up up down];
            else
                y_int = linspace(down, up, 10)';
                xBS = bowshock_pos(y_int);
                X = [xBS + x - 1; flip(xBS) + x];
                Y = [y_int; flip(y_int)];
            end
            patch(X, Y, c);
            axis equal
        end
    end
end