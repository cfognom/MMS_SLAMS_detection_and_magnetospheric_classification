function MMS_load_edp_fast_l2_sc4(starttime,stoptime)

global E4gse;

%% FIX INPUT TIMES
startepochtt = irf_time(starttime,'vector>epochtt');
stopepochtt = irf_time(stoptime,'vector>epochtt');
tint  = irf.tint(startepochtt,stopepochtt);

%% INITIATE DATABASE
%mms.db_init('local_file_db','E:\Data\mms\data\');
mms.db_init('local_file_db','C:\Users\tomask\Documents\Space Physics\Data\mms\data\');

%% LOAD FGM DATA
E4gse = mms.db_get_ts('mms4_edp_fast_l2_dce','mms4_edp_dce_gse_fast_l2',tint);


return;