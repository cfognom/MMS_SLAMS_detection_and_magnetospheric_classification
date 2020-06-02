function tints_active = get_tints_active(tints_user)
    
    tints_active = EpochTT(int64(get_cached_matrix('tints_active', @make_tints_active_matrix)));
    
    function tints_active_matrix = make_tints_active_matrix()
        fgm_time_threshold = 0.125*2;
        fpi_time_threshold = 4.5*2;
        mec_time_threshold = 2*30;
        tints_active_fgm = get_tints_instrument_active('mms1_fgm_srvy_l2',          tints_user, fgm_time_threshold);
        tints_active_fpi = get_tints_instrument_active('mms1_fpi_fast_l2_dis-moms', tints_user, fpi_time_threshold);
        tints_active_mec = get_tints_instrument_active('mms1_mec_srvy_l2_ephts04d', tints_user, mec_time_threshold);
        tints_active = intersect_tints(tints_active_fgm, tints_active_fpi, tints_active_mec);
        tints_active_matrix = tints_active.epoch;
    end
end

