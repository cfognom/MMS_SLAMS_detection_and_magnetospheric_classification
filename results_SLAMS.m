
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

    SLAMS = filter_SLAMS(SLAMS, SLAMS_filter{:});

    n_classes = dataset_info.n_region_classes;
    color_scheme = {[1 0 0], [0 0.8 0], [0 0 1], [0 1 1]};
    ordering = [3, 4, 2, 1];

    % windows: [instant, 15, 30, 60, 120, 240, 480]
    window_idx = 4;
    posterior_threshold = 0.8;
    SLAMS_classes = split_SLAMS(SLAMS, posterior_threshold, window_idx);
    
    
    plot_pos_classes(SLAMS_classes);
    
    histogram_angles(SLAMS, SLAMS_MSP, SLAMS_SW, SLAMS_MSH, SLAMS_FS);
    
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
        for i = ordering
            plot_pos(SLAMS_classes{i}, color_scheme{i})
        end
        plot_reference()
        ylabel('\rho R_{GSE}')
        xlabel('x R_{GSE}')
        axis equal
        grid on

        function plot_reference()
            K = 25;
            epsilon = 0.8;
            a = linspace(0, 2*pi);
            r = K./(1 + epsilon*cos(a));
            x = cos(a);
            y = sin(a);
            plot(x, y, 'k')
            plot(x.*r, y.*r, 'k')
            % pos = [repelem(30, 100)', zeros(100,1), (y.*r)'];
            % d = dist_to_bow_shock(pos);
            % plot(x.*r + d', y.*r, 'k')
        end
    
    end
    
    function plot_pos(SLAMS, c)
        pos = vertcat(SLAMS.pos_GSE)/6378;
        pos_yz = irf_abs(pos(:, 2:3), 1);
        scatter(pos(:, 1), pos_yz, 36, c, '.')
    end
    
    function histogram_angles(SLAMS, SLAMS_MSP, SLAMS_SW, SLAMS_MSH, SLAMS_FS)
        edges = linspace(0,180, 19);
        all_ang = calc_angle(SLAMS);
        MSP_ang = calc_angle(SLAMS_MSP);
        SW_ang = calc_angle(SLAMS_SW);
        MSH_ang = calc_angle(SLAMS_MSH);
        FS_ang = calc_angle(SLAMS_FS);
    
        figure;
        histogram(all_ang, edges, 'FaceColor', [0.5 0.5 0.5]);
        hold on
        histogram(MSH_ang, edges, 'FaceColor', [0 0 1]);
        histogram(FS_ang, edges, 'FaceColor', [0 1 1]);
        histogram(SW_ang, edges, 'FaceColor', [0 1 0]);
        histogram(MSP_ang, edges, 'FaceColor', [1 0 0]);
        ylabel('Count SLAMS')
        xlabel('Angle from [1 0 0] GSE')
        legend('ALL', 'MSH', 'FS', 'SW', 'MSP')
    
        function bool = angle_filt(pos_GSE)
        end
        
        function ang = calc_angle(SLAMS)
            pos_GSE = vertcat(SLAMS.pos_GSE);
            ang = acosd(pos_GSE(:,1)./(sqrt(sum(pos_GSE.^2, 2))));
        end
    end
    
    function plot_strength(SLAMS, SLAMS_MSP, SLAMS_SW, SLAMS_MSH, SLAMS_FS)
        
    end
    
    function d = dist_to_bow_shock(pos)
        K = 25;
        epsilon = 0.8;
        pos_yz = irf_abs(pos(:, 2:3), 1);
        bowshock = 2*(epsilon*K - sqrt(K^2 + 0.25*(4*epsilon^2 - 4)*pos_yz.^2))/(2*epsilon^2 - 2);
        d = pos(:,1) - bowshock;
    end
    
    function SLAMS_classes = split_SLAMS(SLAMS, posterior_threshold, window_idx)
        
        SLAMS_classes = cell(1, n_classes);

        SLAMS.region_posterior

        for i = 1:n_classes
            SLAMS_classes{i} = filter_SLAMS(SLAMS, 'Region_filter', @(p,m,l) p(:,i,window_idx) > posterior_threshold);
        end

        function bool = probable_class(p, m, l)
            p(:,i,window_idx)
        end
    end
end