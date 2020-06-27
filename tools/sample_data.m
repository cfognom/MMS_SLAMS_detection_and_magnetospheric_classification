function sample_data()

    % Time intervals when MMS apogee is touching the magnetopause
    % Taken from MMS science data center orbit plots
    time_intervals = {
        '2015-10-01T00:00:00.000000000Z', '2016-01-30T00:00:00.000000000Z';
        '2016-11-01T00:00:00.000000000Z', '2017-05-09T00:00:00.000000000Z';
        '2017-09-08T00:00:00.000000000Z', '2018-05-23T00:00:00.000000000Z';
        '2018-09-23T00:00:00.000000000Z', '2019-06-16T00:00:00.000000000Z';
        '2019-09-26T00:00:00.000000000Z', '2020-06-08T00:00:00.000000000Z';
    };

    tints_search = UTC2tints(time_intervals);

    tints_active = get_tints_active(tints_search, 'sample');

    rng(100)
    n = 10000;
    t_train = sample_tints(tints_active, n);
    % starts = [SLAMS.start];
    % stops = [SLAMS.stops];
    % t_train = EpochTT(int64((starts.epoch + stops.epoch)/2));
    X_t = double(sort(t_train).epoch);

    % randSampData
    fpi_time_window = 4.5*4;
    mec_time_window = 2*30;
    save_path = 'cache/randSampData.txt';
    % save_path = 'cache/SLAMSSampData.txt';
    check();
    X_v_ion = get_data_samples('mms1_fpi_fast_l2_dis-moms', 'mms1_dis_bulkv_gse_fast', t_train, fpi_time_window);
    X_E_ion_spectr = get_data_samples('mms1_fpi_fast_l2_dis-moms', 'mms1_dis_energyspectr_omni_fast', t_train, fpi_time_window);
    X_E_ion_delta = get_data_samples('mms1_fpi_fast_l2_dis-moms', 'mms1_dis_energy_delta_fast', t_train, fpi_time_window);
    X_E_ion = get_data_samples('mms1_fpi_fast_l2_dis-moms', 'mms1_dis_energy_fast', t_train, fpi_time_window);
    X_gse = get_data_samples('mms1_mec_srvy_l2_ephts04d', 'mms1_mec_r_gse', t_train, mec_time_window);
    X = [X_t, X_v_ion, X_E_ion_spectr, X_E_ion_delta, X_E_ion, X_gse];

    writematrix(X, save_path);

    function check()
        if isfile(save_path)
            inp = input(['''', save_path, ''' already exists, do you want to replace it? y/n [n]: '], 's');
            if isempty(inp) || inp ~= 'y'
                return
            end
        end
    end
end

