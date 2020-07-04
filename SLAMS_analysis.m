function SLAMS_analysis(SLAMS_database)
    
    sc = 'mms1';

    time_intervals = {
        '2015-10-01T00:00:00.000000000Z', '2016-01-30T00:00:00.000000000Z';
        '2016-11-01T00:00:00.000000000Z', '2017-05-09T00:00:00.000000000Z';
        '2017-09-08T00:00:00.000000000Z', '2018-05-23T00:00:00.000000000Z';
        '2018-09-23T00:00:00.000000000Z', '2019-06-16T00:00:00.000000000Z';
        '2019-09-26T00:00:00.000000000Z', '2020-06-08T00:00:00.000000000Z';
    };

    n = 1000;

    SLAMS = SLAMS_database.SLAMS;

    window_idx = 3;

    posteriors = cat(3, SLAMS.region_posterior);
    posteriors = permute(posteriors(window_idx,:,:), [3, 2, 1]);
    [~, idx] = max(posteriors, [], 2);

    SLAMS_MSH = pick_n_SLAMS(SLAMS(idx == 3), n);
    SLAMS_FS = pick_n_SLAMS(SLAMS(idx == 4), n);
    % SLAMS_SW = pick_n_SLAMS(SLAMS(idx == 2), 472);

    % n_SLAMS_MSH = length(SLAMS_MSH);
    % n_SLAMS_FS = length(SLAMS_FS);

    starts_MSH = [SLAMS_MSH.start];
    stops_MSH = [SLAMS_MSH.stop];
    starts_FS = [SLAMS_FS.start];
    stops_FS = [SLAMS_FS.stop];
    % starts_SW = [SLAMS_SW.start];
    % stops_SW = [SLAMS_SW.stop];

    SLAMS_MSH_tints = [starts_MSH; stops_MSH];
    SLAMS_MSH_tints = reorder_tints(SLAMS_MSH_tints);
    SLAMS_FS_tints = [starts_FS; stops_FS];
    SLAMS_FS_tints = reorder_tints(SLAMS_FS_tints);

    SLAMS_MSH_mids = EpochTT(int64(([starts_MSH.epoch] + [stops_MSH.epoch])/2));
    SLAMS_FS_mids = EpochTT(int64(([starts_FS.epoch] + [stops_FS.epoch])/2));
    % SLAMS_SW_mids = EpochTT(int64(([starts_SW.epoch] + [stops_SW.epoch])/2));
    tints_active = get_tints_active(UTC2tints(time_intervals), 'search');
    control_points = sample_tints(tints_active, n);

    dt = 60;

    tints_SLAMS_MSH = point_to_tint(SLAMS_MSH_mids, dt);
    tints_SLAMS_FS = point_to_tint(SLAMS_FS_mids, dt);
    % tints_SLAMS_SW = point_to_tint(SLAMS_SW_mids, dt);
    tints_control = point_to_tint(control_points, dt);

    Ma_v = get_cached_matrix('Ma_SLAMS_MSH', @() process_tints(tints_SLAMS_MSH, SLAMS_MSH_tints));
    Ma_SLAMS_MSH = Ma_v(:,1);
    v_MSH = Ma_v(:,2);
    Ma_v = get_cached_matrix('Ma_SLAMS_FS', @() process_tints(tints_SLAMS_FS, SLAMS_FS_tints));
    Ma_SLAMS_FS = Ma_v(:,1);
    v_FS = Ma_v(:,2);
    % Ma_SLAMS_SW = get_cached_matrix('Ma_SLAMS_SW', @() process_tints(tints_SLAMS_SW));
    Ma_v = get_cached_matrix('Ma_control', @() process_tints(tints_control));
    Ma_control = Ma_v(:,1);

    edges = 0:1:20;

    figure;
    hold on
    grid on
    histogram(Ma_control, edges, 'FaceColor', 'k')
    histogram(Ma_SLAMS_MSH, edges, 'FaceColor', 'b')
    histogram(Ma_SLAMS_FS, edges, 'FaceColor', 'c')
    % histogram(Ma_SLAMS_SW, edges, 'FaceColor', 'g')
    ylabel('SLAMS count [n]')
    xlabel('Mach number of surrounding flow [-]')
    title('SLAMS mach number')
    xticks(edges(1:2:end))
    legend('Control', 'MSH SLAMS', 'FS SLAMS')

    delta_t_MSH = [SLAMS_MSH.stop] - [SLAMS_MSH.start];
    delta_t_FS = [SLAMS_FS.stop] - [SLAMS_FS.start];
    size_SLAMS_MSH_nans = (delta_t_MSH).*v_MSH;
    size_SLAMS_MSH = size_SLAMS_MSH_nans(~isnan(size_SLAMS_MSH_nans));
    size_SLAMS_FS_nans = (delta_t_FS).*v_FS;
    size_SLAMS_FS = size_SLAMS_FS_nans(~isnan(size_SLAMS_FS_nans));
    
    l_MSH = length(size_SLAMS_MSH);
    l_FS = length(size_SLAMS_FS);
    
    l_min = min(l_MSH, l_FS);
    size_SLAMS_MSH = size_SLAMS_MSH(1:l_min);
    size_SLAMS_FS = size_SLAMS_FS(1:l_min);
    
    disp(['MSH SLAMS count for size plot = ', num2str(length(size_SLAMS_MSH))])
    disp(['FS SLAMS count for size plot = ', num2str(length(size_SLAMS_FS))])

    edges = 0:500:8000;
    
    figure;
    hold on
    grid on
    histogram(size_SLAMS_MSH, edges, 'FaceColor', 'b')
    histogram(size_SLAMS_FS, edges, 'FaceColor', 'c')
    ylabel('Count [n]')
    xlabel('Size [km]')
    title('SLAMS size')
    xticks(edges(1:2:end))
    legend('MSH SLAMS', 'FS SLAMS')
    
    edges = 0:1:20;

    figure;
    hold on
    grid on
    histogram(delta_t_MSH, edges, 'FaceColor', 'b')
    histogram(delta_t_FS, edges, 'FaceColor', 'c')
    ylabel('Count [n]')
    xlabel('Duration [s]')
    title('SLAMS duration')
    xticks(edges(1:2:end))
    legend('MSH SLAMS', 'FS SLAMS')

    edges = 2:0.1:4;

    figure;
    hold on
    grid on
    histogram([SLAMS_MSH.B_rel_max], edges, 'FaceColor', 'b')
    histogram([SLAMS_FS.B_rel_max], edges, 'FaceColor', 'c')
    ylabel('Count [n]')
    xlabel('B_{max}/B_0 [-]')
    title('SLAMS relative strength')
    xticks(edges(1:2:end))
    legend('MSH SLAMS', 'FS SLAMS')

    edges = 0:1:20;

    figure;
    hold on
    grid on
    counts = histogram_mean(Ma_SLAMS_MSH, [SLAMS_MSH.B_rel_max], edges);
    histogram('BinEdges', edges, 'BinCounts', counts, 'FaceColor', 'b')
    counts = histogram_mean(Ma_SLAMS_FS, [SLAMS_FS.B_rel_max], edges);
    histogram('BinEdges', edges, 'BinCounts', counts, 'FaceColor', 'c')
    xlabel('Mach number of surrounding flow [-]')
    ylabel('Mean B_{max}/B_0 [-]')
    title('SLAMS relative strength vs mach number')
    legend('MSH SLAMS', 'FS SLAMS', 'Location', 'NorthWest')

    figure;
    hold on
    grid on
    counts = histogram_mean(Ma_SLAMS_MSH, [SLAMS_MSH.B_max], edges);
    histogram('BinEdges', edges, 'BinCounts', counts, 'FaceColor', 'b')
    counts = histogram_mean(Ma_SLAMS_FS, [SLAMS_FS.B_max], edges);
    histogram('BinEdges', edges, 'BinCounts', counts, 'FaceColor', 'c')
    xlabel('Mach number of surrounding flow [-]')
    ylabel('Mean B_{max} [nT]')
    title('SLAMS strength vs mach number')
    legend('MSH SLAMS', 'FS SLAMS')

    figure;
    hold on
    grid on
    counts = histogram_mean(Ma_SLAMS_MSH, size_SLAMS_MSH_nans, edges);
    histogram('BinEdges', edges, 'BinCounts', counts, 'FaceColor', 'b')
    counts = histogram_mean(Ma_SLAMS_FS, size_SLAMS_FS_nans, edges);
    histogram('BinEdges', edges, 'BinCounts', counts, 'FaceColor', 'c')
    xlabel('Mach number of surrounding flow [-]')
    ylabel('Mean size [km]')
    title('SLAMS size vs mach number')
    legend('MSH SLAMS', 'FS SLAMS')

    figure;
    hold on
    grid on
    counts = histogram_mean(Ma_SLAMS_MSH, delta_t_MSH, edges);
    histogram('BinEdges', edges, 'BinCounts', counts, 'FaceColor', 'b')
    counts = histogram_mean(Ma_SLAMS_FS, delta_t_FS, edges);
    histogram('BinEdges', edges, 'BinCounts', counts, 'FaceColor', 'c')
    xlabel('Mach number of surrounding flow [-]')
    ylabel('Mean duration [s]')
    title('SLAMS duration vs mach number')
    legend('MSH SLAMS', 'FS SLAMS')

    function SLAMS = pick_n_SLAMS(SLAMS, n)
        rng(100)
        rdm_idx = randperm(length(SLAMS), n);
        SLAMS = SLAMS(rdm_idx);
    end

    function tints = point_to_tint(t, dt)
        tints = [t + -dt/2; t + dt/2];
        ordering = reshape(reshape(1:length(tints), [], 2)', 1, []);
        tints = tints(ordering);
    end

    function Ma_v = process_tints(tints, tints_SLAMS)
        m_ion = 1.6726219*1e-27;
        mu_0 = 1.25663706212*1e-6;

        n = length(tints)/2;

        Ma_v = zeros(n, 2);
        for i = 1:n
            fprintf('Processing tint %u/%u\n', i, n)
            tint = select_tint(tints, i);
            if nargin == 1
                tint_SLAMS = [];
            else
                tint_SLAMS = select_tint(tints_SLAMS, i);
            end
            [numberdensity, bulkv, b, Ma_v(i,2)] = load_data(tint, tint_SLAMS);

            alfven_speed = (b*1e-9)/sqrt(mu_0*numberdensity*1e6*m_ion);
            Ma_v(i,1) = (bulkv*1e3)/alfven_speed;
        end

        function [numberdensity, bulkv, b, v] = load_data(tint, tint_SLAMS)
            % tint
            numberdensity = load_tints_MMS([sc, '_fpi_fast_l2_dis-moms'], [sc, '_dis_numberdensity_fast'], tint);
            numberdensity = median(numberdensity.data);
            bulkv = load_tints_MMS([sc, '_fpi_fast_l2_dis-moms'], [sc, '_dis_bulkv_gse_fast'], tint);
            if isempty(tint_SLAMS)
                bulkv_lim = [];
            else
                bulkv_lim = bulkv.tlim(tint_SLAMS);
            end
            % bulkv
            % length(bulkv)
            % ceil(length(bulkv)/2)
            % idx = ceil(length(bulkv)/2);
            if isempty(bulkv_lim)
                v = NaN;
            else
                v = mean(irf_abs(bulkv_lim.data, 1));
            end
            bulkv = median(irf_abs(bulkv.data, 1));
            b = load_tints_MMS([sc, '_fgm_srvy_l2'], [sc, '_fgm_b_gse_srvy_l2'], tint);
            b = median(irf_abs(b.data, 1));
        end
    end

    function tints = reorder_tints(tints)
        l = length(tints);
        ord = reshape(reshape(1:l, [], 2)', 1, [])';
        tints = tints(ord);
    end

    function counts = histogram_mean(x, y, edges)
        idx_nan = isnan(y);
        x(idx_nan) = [];
        y(idx_nan) = [];
        [~, ~, bins] = histcounts(x, edges);
        n_bins = length(edges) - 1;
        counts = zeros(1, n_bins);
        for j = 1:n_bins
            counts(j) = mean(y(bins == j));
        end
        counts(isnan(counts)) = 0;
    end
end