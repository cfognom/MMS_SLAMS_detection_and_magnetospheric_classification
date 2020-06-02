force_load_ri = [249];
% force_load_ri = [];
dt_view = 30; % Extra time in seconds before and after to view.
dt_load = 60*20; % Extra time in seconds before and after to load.

[tints, descs] = read_events('events/slams_burst.csv');
n_tints = length(tints)/2;

if isempty(force_load_ri)
    rng(100)
    random_order = randi(n_tints, n_tints, 1);
    ris = random_order;
else
    ris = force_load_ri;
end

n_ris = length(ris);
for i = 1:n_ris
    ri = ris(i);
    tint = select_tint(tints, ri);
    desc = descs{ri};
    file_name_data = ['labeled/data/', sprintf('%u_%u_%u', ri, tint(1).epoch, tint(2).epoch)];
    file_name_plot = ['labeled/plot/', sprintf('%u_%u_%u', ri, tint(1).epoch, tint(2).epoch)];
    try
        pre_marked = readmatrix(file_name_data);
        existing_file = true;
    catch
        pre_marked = [];
        existing_file = false;
    end
    if existing_file && isempty(force_load_ri) % File exists but it should not be forced loaded
        fprintf('Time interval with i = %u and ri = %u already done.\n', i, ri);
        continue;
    end
    try
        if ~mark_ri(tint, desc, ri, file_name_data, file_name_plot, dt_view, dt_load, existing_file, pre_marked)
            continue
        end
    catch
        break;
    end
end

function flag = mark_ri(tint, desc, ri, file_name_data, file_name_plot, dt_view, dt_load, existing_file, pre_marked)

    t_span = duration(0, 0, tint(2) - tint(1));
    tint_view = [tint(1) + -dt_view, tint(2) + dt_view];
    tint_load = [tint(1) + -dt_load, tint(2) + dt_load];
    
    % Load fgm
    fgm_filePrefix = 'mms1_fgm_srvy_l2';
    b = load_tints_MMS(fgm_filePrefix, 'mms1_fgm_b_gse_srvy_l2', tint_load);
    b_abs = irf_abs(b);
    
    % Load ion-fpi
    ion_fpi_filePrefix = 'mms1_fpi_fast_l2_dis-moms';
    v_ion = load_tints_MMS(ion_fpi_filePrefix, 'mms1_dis_bulkv_gse_fast', tint_load);
    v_ion_abs = irf_abs(v_ion);
    n_ion = load_tints_MMS(ion_fpi_filePrefix, 'mms1_dis_numberdensity_fast', tint_load);
    T_ion_para = load_tints_MMS(ion_fpi_filePrefix, 'mms1_dis_temppara_fast', tint_load);
    T_ion_perp = load_tints_MMS(ion_fpi_filePrefix, 'mms1_dis_tempperp_fast', tint_load);
    E_ion_spectr = load_tints_MMS(ion_fpi_filePrefix, 'mms1_dis_energyspectr_omni_fast', tint_load);
    E_ion = load_tints_MMS(ion_fpi_filePrefix, 'mms1_dis_energy_fast', tint_load);
    
    % Load ele-fpi
    ele_fpi_filePrefix = 'mms1_fpi_fast_l2_des-moms';
    n_ele = load_tints_MMS(ele_fpi_filePrefix, 'mms1_des_numberdensity_fast', tint_load);
    
    % Load mec
    mec_filePrefix = 'mms1_mec_srvy_l2_ephts04d';
    pos = load_tints_MMS(mec_filePrefix, 'mms1_mec_r_gse', tint_load)/6378;
    
    if isempty(b) || isempty(v_ion)
        fprintf('Warning: No files found for time interval with ri = %u. Skipping.\n', ri);
        flag = false;
        return
    end

    plot_title = sprintf('MMS1, timespan = %s', t_span);
    plt = modular_plot('title', plot_title);
    plt.lineplot('B', {b, b_abs}, 'ylabel', 'B_{GSE} (nT)');
    plt.lineplot('n', {n_ion, n_ele}, 'ylabel', 'n_i (cm^{-3})')
    plt.lineplot('T_ion', {T_ion_para, T_ion_perp}, 'ylabel', 'T_{i} (eV)', 'legend', {'T_{ipara}','T_{iperp}'})
    plt.spectrogram('E_ion_spectr', {E_ion_spectr, E_ion}, 'ylabel', 'W_i (eV)')
    plt.lineplot('v_ion', {v_ion, v_ion_abs}, 'ylabel', 'v_{i} (kms^{-1})', 'legend', {'v_{ix}','v_{iy}','v_{iz}'})
    plt.lineplot('pos', pos, 'ylabel', 'R_{GSE} (R_E)', 'legend', {'x','y','z'})
    % plt.lineplot('class_score', score_TS, 'ylabel', 'Score', 'legend', {'Msphere', 'SW','Msheath'}, 'colorOrder', [1 0 0; 0 0.5 0; 0 0 1])
    plt.show(tint_view);
    plt.mark(tint, 'target', 'B', 'color', 'b');
    fprintf('Info:\n\tri = %u\n', ri);
    fprintf('Description:\n\t%s\n', desc);
    fprintf('Controls:\n\tLMouse => Select point\n\tBackspace/RMouse => Delete nearest\n\tArrowKeys => Zoom and pan\n\tR => Restore view\n\tReturn => When done\n');
    tints_user = plt.ginput_mark('target', 'B', 'intervals', true, 'color', [1 0.5 0.5], 'preMark', pre_marked);
    if existing_file
        file_name_data = [file_name_data, '_mod'];
        file_name_plot = [file_name_plot, '_mod'];
    end
    plt.save([file_name_plot, '.png'], tint_view);
    fprintf('Plot saved to ''%s''\n', file_name_plot);
    close(plt.fig);
    if ~isempty(tints_user)
        writematrix(tints_user.epoch, file_name_data);
    else
        writematrix([], file_name_data)
    end
    fprintf('Data saved to ''%s''\n', file_name_data);
    flag = true;
end