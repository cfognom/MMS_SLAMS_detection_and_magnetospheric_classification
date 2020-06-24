function [tints_active, tints_fgm, tints_fpi, tints_mec] = get_tints_active(tints_user, name)
    
    tints_fgm = EpochTT(int64(get_cached_matrix(['tints_fgm_', name], @make_tints_fgm_matrix)));
    tints_fpi = EpochTT(int64(get_cached_matrix(['tints_fpi_', name], @make_tints_fpi_matrix)));
    tints_mec = EpochTT(int64(get_cached_matrix(['tints_mec_', name], @make_tints_mec_matrix)));

    tints_active = intersect_tints(tints_fgm, tints_fpi, tints_mec);
    
    function tints_fgm_matrix = make_tints_fgm_matrix()
        fgm_time_threshold = 0.125*2;
        tints_fgm_matrix = get_tints_instrument_active('mms1_fgm_srvy_l2', tints_user, fgm_time_threshold);
        tints_fgm_matrix = tints_fgm_matrix.epoch;
    end

    function tints_fpi_matrix = make_tints_fpi_matrix()
        fpi_time_threshold = 4.5*2;
        tints_fpi_matrix = get_tints_instrument_active('mms1_fpi_fast_l2_dis-moms', tints_user, fpi_time_threshold);
        tints_fpi_matrix = tints_fpi_matrix.epoch;
    end

    function tints_mec_matrix = make_tints_mec_matrix()
        mec_time_threshold = 2*30;
        tints_mec_matrix = get_tints_instrument_active('mms1_mec_srvy_l2_ephts04d', tints_user, mec_time_threshold);
        tints_mec_matrix = tints_mec_matrix.epoch;
    end
end

