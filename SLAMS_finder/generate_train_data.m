function generate_train_data()
%GENERATE_TRAIN_DATA Summary of this function goes here
%   Detailed explanation goes here

    tints_search_user = get_tints_user();

    tints_active = get_tints_active(tints_search_user);

    rng(100)
    n = 10000;
    t_train = sample_tints(tints_active, n);
    X_t = double(sort(t_train).epoch);

    % % traindata1
    % fgm_time_window = 0.125*4;
    % fpi_time_window = 4.5*4;
    % save_path = 'cache/traindata1.txt';
    % check();
    % X_b = get_data_samples('mms1_fgm_srvy_l2',          'mms1_fgm_b_gse_srvy_l2',      t_train, fgm_time_window);
    % X_v = get_data_samples('mms1_fpi_fast_l2_dis-moms', 'mms1_dis_bulkv_gse_fast',     t_train, fpi_time_window);
    % X_n = get_data_samples('mms1_fpi_fast_l2_dis-moms', 'mms1_dis_numberdensity_fast', t_train, fpi_time_window);
    % X_Tpara = get_data_samples('mms1_fpi_fast_l2_dis-moms', 'mms1_dis_temppara_fast', t_train, fpi_time_window);
    % X_Tperp = get_data_samples('mms1_fpi_fast_l2_dis-moms', 'mms1_dis_tempperp_fast', t_train, fpi_time_window);
    % X = [X_t, X_b, X_v, X_n, X_Tpara, X_Tperp];

    % % traindata2
    % fgm_time_window = 120;
    % fpi_time_window = 120;
    % save_path = 'cache/traindata2.txt';
    % check();
    % X_b_vari_mean = get_data_samples('mms1_fgm_srvy_l2',          'mms1_fgm_b_gse_srvy_l2',      t_train, fgm_time_window, 'modfunc', @vari_mean);
    % X_v_vari_mean = get_data_samples('mms1_fpi_fast_l2_dis-moms', 'mms1_dis_bulkv_gse_fast',     t_train, fpi_time_window, 'modfunc', @vari_mean);
    % X_n_vari_mean = get_data_samples('mms1_fpi_fast_l2_dis-moms', 'mms1_dis_numberdensity_fast', t_train, fpi_time_window, 'modfunc', @vari_mean);
    % X = [X_t, X_b_vari_mean, X_v_vari_mean, X_n_vari_mean];

    % % traindata3
    % fpi_time_window = 4.5*4;
    % save_path = 'cache/traindata3.txt';
    % check();
    % X_E_ion_spectr = get_data_samples('mms1_fpi_fast_l2_dis-moms', 'mms1_dis_energyspectr_omni_fast', t_train, fpi_time_window);
    % X_E_ion = get_data_samples('mms1_fpi_fast_l2_dis-moms', 'mms1_dis_energy_fast', t_train, fpi_time_window);
    % X = [X_t, X_E_ion, X_E_ion_spectr];

    % % traindata4
    % fpi_time_window = 4.5*4;
    % mec_time_window = 2*30;
    % save_path = 'cache/traindata4.txt';
    % check();
    % X_v_ion = get_data_samples('mms1_fpi_fast_l2_dis-moms', 'mms1_dis_bulkv_gse_fast', t_train, fpi_time_window);
    % X_E_ion_spectr = get_data_samples('mms1_fpi_fast_l2_dis-moms', 'mms1_dis_energyspectr_omni_fast', t_train, fpi_time_window);
    % X_E_ion = get_data_samples('mms1_fpi_fast_l2_dis-moms', 'mms1_dis_energy_fast', t_train, fpi_time_window);
    % X_gse = get_data_samples('mms1_mec_srvy_l2_ephts04d', 'mms1_mec_r_gse', t_train, mec_time_window);
    % X = [X_t, X_v_ion, X_E_ion_spectr, X_E_ion, X_gse];

    % traindata5
    fpi_time_window = 4.5*4;
    save_path = 'cache/traindata5.txt';
    check();
    X_E_ion_delta = get_data_samples('mms1_fpi_fast_l2_dis-moms', 'mms1_dis_energy_delta_fast', t_train, fpi_time_window);
    X_heat_ion = get_data_samples('mms1_fpi_fast_l2_dis-moms', 'mms1_dis_heatq_gse_fast', t_train, fpi_time_window);
    X = [X_t, X_E_ion_delta, X_heat_ion];

    writematrix(X, save_path);

    function check()
        if isfile(save_path)
            inp = input(['''', save_path, ''' already exists, do you want to replace it? y/n [n]: '], 's');
            if isempty(inp) || inp ~= 'y'
                return
            end
        end
    end

    function out = vari_mean(ts, t, dt)
        n = length(t);
        [~, d] = size(ts.data);
        out = zeros(n, d*2);
        for i = 1:n
            tint = [t(i) + -dt, t(i) + dt];
            ts_lim = ts.tlim(tint);
            out(i, 1:d) = double(var(ts_lim.data, 0, 1));
            out(i, d+1:d*2) = double(mean(ts_lim.data, 1));
        end
    end
end

