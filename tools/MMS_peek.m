function MMS_peek()
%MMS_PEEK Summary of this function goes here
%   Detailed explanation goes here

    % tint_user = [EpochTT('2018-11-13T10:35:53.0Z'), EpochTT('2018-11-13T10:40:13.0Z')]; % SLAMS 1
    % tint_user = [EpochTT('2017-10-20T23:09:00.0Z'), EpochTT('2017-10-20T23:15:00.0Z')]; % SLAMS 2 % Den
    % tint_user = [EpochTT('2017-11-23T04:54:23Z'), EpochTT('2017-11-23T04:56:43Z')]; % SLAMS 3
    % tint_user = [EpochTT('2017-11-23T01:54:23Z'), EpochTT('2017-11-23T05:56:43Z')]; % SLAMS 4 % och Den
    % tint_user = [EpochTT('2018-01-09T05:00:00Z'), EpochTT('2018-01-09T12:00:00Z')]; % Des
    % tint_user = [EpochTT('2019-04-13T16:50:00Z'), EpochTT('2019-04-13T17:30:00Z')]; % HFA?
    % tint_user = [EpochTT('2017-12-21T08:00:00Z'), EpochTT('2017-12-21T09:00:00Z')];
    % tint_user = [EpochTT('2017-11-23T05:00:00Z'), EpochTT('2017-11-23T06:00:00Z')];
    tint_user = [EpochTT('2018-02-03T12:00:00.0Z'), EpochTT('2018-02-04T00:00:00.0Z')]; % 3 boundaries
    % tint_user = [EpochTT('2018-01-08T00:00:00.0Z'), EpochTT('2018-01-08T12:00:00.0Z')]; % Interplanetary shock
    % tint_user = [EpochTT('2018-01-07T00:00:00.0Z'), EpochTT('2018-01-07T10:00:00.0Z')]; % Before interplanetary shock
    % tint_user = [EpochTT('2019-02-26T00:30:00.0Z'), EpochTT('2019-02-26T01:30:00.0Z')]; % quasi-para bowshock
    % tint_user = [EpochTT('2015-11-06T02:30:38.0Z'), EpochTT('2015-11-06T15:45:13.0Z')]; % Magnetopause
    % tint_user = [EpochTT('2018-10-01T00:00:00.0Z'), EpochTT('2018-10-04T00:00:00.0Z')]; % test long tint
    % tint_user = [EpochTT('2018-10-01T00:00:00.0Z'), EpochTT('2018-11-01T00:00:00.0Z')]; % test one month
    % tint_user = [EpochTT('2017-01-18T09:00:00.0Z'), EpochTT('2017-01-18T11:00:00.0Z')]; % noisy
    % tint_user = [EpochTT('2018-02-24T00:00:00.0Z'), EpochTT('2018-02-24T12:00:00.0Z')]; % 20180224
    % random_t = sample_t_from_tints(tints_active, 1);
    % dt = 60*60;
    % tint_user = [random_t + -dt/2, random_t + dt/2]; % random time

    tints_load = tint_user;

    % Load fgm
    fgm_filePrefix = 'mms1_fgm_srvy_l2';
    b = load_tints_MMS(fgm_filePrefix, 'mms1_fgm_b_gse_srvy_l2', tints_load);
    b_abs = irf_abs(b);

    % Load ion-fpi
    ion_fpi_filePrefix = 'mms1_fpi_fast_l2_dis-moms';
    v_ion = load_tints_MMS(ion_fpi_filePrefix, 'mms1_dis_bulkv_gse_fast', tints_load);
    v_ion_abs = irf_abs(v_ion);
    n_ion = load_tints_MMS(ion_fpi_filePrefix, 'mms1_dis_numberdensity_fast', tints_load);
    T_ion_para = load_tints_MMS(ion_fpi_filePrefix, 'mms1_dis_temppara_fast', tints_load);
    T_ion_perp = load_tints_MMS(ion_fpi_filePrefix, 'mms1_dis_tempperp_fast', tints_load);
    E_ion_spectr = load_tints_MMS(ion_fpi_filePrefix, 'mms1_dis_energyspectr_omni_fast', tints_load);
    E_ion = load_tints_MMS(ion_fpi_filePrefix, 'mms1_dis_energy_fast', tints_load);

    % Load ele-fpi
    ele_fpi_filePrefix = 'mms1_fpi_fast_l2_des-moms';
    n_ele = load_tints_MMS(ele_fpi_filePrefix, 'mms1_des_numberdensity_fast', tints_load);
    % T_ele_para = load_tints_MMS(ele_fpi_filePrefix, 'mms1_des_temppara_fast', tints_load);
    % T_ele_perp = load_tints_MMS(ele_fpi_filePrefix, 'mms1_des_tempperp_fast', tints_load);

    % Load mec
    mec_filePrefix = 'mms1_mec_srvy_l2_ephts04d';
    pos = load_tints_MMS(mec_filePrefix, 'mms1_mec_r_gse', tints_load)/6378;

    % T_ion_para_var = TSeries(T_ion_para.time, movvar(T_ion_para.data, 1*60*1e9, 1, 'SamplePoints', T_ion_para.time.epoch));
    % T_ion_perp_var = TSeries(T_ion_perp.time, movvar(T_ion_perp.data, 1*60*1e9, 1, 'SamplePoints', T_ion_perp.time.epoch));
    % T_ele_para_var = TSeries(T_ele_para.time, movvar(T_ele_para.data, 1*60*1e9, 1, 'SamplePoints', T_ele_para.time.epoch));
    % T_ele_perp_var = TSeries(T_ele_perp.time, movvar(T_ele_perp.data, 1*60*1e9, 1, 'SamplePoints', T_ele_perp.time.epoch));

    plt = modular_plot('title', 'MMS1');
    plt.lineplot('B', {b, b_abs}, 'ylabel', 'B_{GSE} (nT)');
    plt.lineplot('n', {n_ion, n_ele}, 'ylabel', 'n_i (cm^{-3})')
    plt.lineplot('T_ion', {T_ion_para, T_ion_perp}, 'ylabel', 'T_{i} (eV)', 'legend', {'T_{ipara}','T_{iperp}'})
    % plt.lineplot('T_ele', {T_ele_para, T_ele_perp}, 'ylabel', 'T_{e} (eV)', 'legend', {'T_{epara}','T_{eperp}'})
    % plt.lineplot('T_ion_var', {T_ion_para_var, T_ion_perp_var}, 'ylabel', 'T_{i} (eV)', 'legend', {'T_{ipara}','T_{iperp}'})
    % plt.lineplot('T_ele_var', {T_ele_para_var, T_ele_perp_var}, 'ylabel', 'T_{e} (eV)', 'legend', {'T_{epara}','T_{eperp}'})
    plt.spectrogram('E_ion_spectr', {E_ion_spectr, E_ion}, 'ylabel', 'W_i (eV)')

    % bins = E_ion.data;
    bins = (1:32)/32;
    tot = sum(E_ion_spectr.data, 2);
    X_e_m = sum(E_ion_spectr.data.*bins, 2)./tot;
    X_e_v = sum(E_ion_spectr.data.*(abs(bins - X_e_m)), 2)./tot;
    test1 = TSeries(E_ion_spectr.time, X_e_m);
    test2 = TSeries(E_ion_spectr.time, X_e_v);
    plt.lineplot('test1', test1, 'ylabel', 'test1')
    plt.lineplot('test2', test2, 'ylabel', 'test2')

    plt.lineplot('v_ion', {v_ion, v_ion_abs}, 'ylabel', 'v_{i} (kms^{-1})', 'legend', {'v_{ix}','v_{iy}','v_{iz}'})
    plt.lineplot('pos', pos, 'ylabel', 'R_{GSE} (R_E)', 'legend', {'x','y','z'})
    plt.show(tints_load);
end

