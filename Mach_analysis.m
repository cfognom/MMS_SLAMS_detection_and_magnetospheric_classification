function example_analysis(SLAMS_database)
    
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

    M_a_v = get_cached_matrix('M_a_SLAMS_MSH', @() process_tints(tints_SLAMS_MSH));
    M_a_SLAMS_MSH = M_a_v(:,1);
    v_MSH = M_a_v(:,1);
    M_a_v = get_cached_matrix('M_a_SLAMS_FS', @() process_tints(tints_SLAMS_FS));
    M_a_SLAMS_FS = M_a_v(:,1);
    v_FS = M_a_v(:,1);
    % M_a_SLAMS_SW = get_cached_matrix('M_a_SLAMS_SW', @() process_tints(tints_SLAMS_SW));
    M_a_v = get_cached_matrix('M_a_control', @() process_tints(tints_control));
    M_a_control = M_a_v(:,1);

    edges = 0:1:15;

    figure;
    hold on
    grid on
    histogram(M_a_control, edges, 'FaceColor', 'k')
    histogram(M_a_SLAMS_MSH, edges, 'FaceColor', 'b')
    histogram(M_a_SLAMS_FS, edges, 'FaceColor', 'c')
    % histogram(M_a_SLAMS_SW, edges, 'FaceColor', 'g')
    ylabel('Count')
    xlabel('Mach number')
    title(['Mach number'])
    xticks(edges)
    legend('Control', 'MSH SLAMS', 'FS SLAMS')

    figure;
    hold on
    % scatter(M_a_SLAMS_MSH, [SLAMS_MSH.B_rel_max], 36, 'b', '.')
    % scatter(M_a_SLAMS_FS, [SLAMS_FS.B_rel_max], 36, 'c', '.')
    scatter(M_a_SLAMS_MSH, ([SLAMS_MSH.stop] - [SLAMS_MSH.start]).*v_MSH, 36, 'b', '.')
    scatter(M_a_SLAMS_FS, ([SLAMS_FS.stop] - [SLAMS_FS.start]).*v_FS, 36, 'c', '.')

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

    function M_a_v = process_tints(tints)
        m_ion = 1.6726219*1e-27;
        mu_0 = 1.25663706212*1e-6;

        n = length(tints)/2;

        M_a_v = zeros(n, 2);
        for i = 1:n
            fprintf('Processing tint %u/%u\n', i, n)
            tint = select_tint(tints, i);
            [numberdensity, bulkv, b, M_a_v(i,2)] = load_data(tint);

            alfven_speed = (b*1e-9)/sqrt(mu_0*numberdensity*1e6*m_ion);
            M_a_v(i,1) = (bulkv*1e3)/alfven_speed;
        end

        function [numberdensity, bulkv, b, v] = load_data(tint)
            % tint
            numberdensity = load_tints_MMS([sc, '_fpi_fast_l2_dis-moms'], [sc, '_dis_numberdensity_fast'], tint);
            numberdensity = median(numberdensity.data);
            bulkv = load_tints_MMS([sc, '_fpi_fast_l2_dis-moms'], [sc, '_dis_bulkv_gse_fast'], tint);
            bulkv = irf_abs(bulkv.data, 1);
            % bulkv
            % length(bulkv)
            % ceil(length(bulkv)/2)
            idx = ceil(length(bulkv)/2);
            if idx == 0
                v = NaN;
            else
                v = bulkv(idx);
            end
            bulkv = median(bulkv);
            b = load_tints_MMS([sc, '_fgm_srvy_l2'], [sc, '_fgm_b_gse_srvy_l2'], tint);
            b = median(irf_abs(b.data, 1));
        end
    end
end