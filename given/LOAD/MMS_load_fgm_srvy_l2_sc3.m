function MMS_load_fgm_srvy_l2_sc3(starttime,stoptime)

global B3gse;

%% FIX INPUT TIMES
startepochtt = irf_time(starttime,'vector>epochtt');
stopepochtt = irf_time(stoptime,'vector>epochtt');
tint  = irf.tint(startepochtt,stopepochtt);

%% INITIATE DATABASE
%mms.db_init('local_file_db','E:\Data\mms\data\');
mms.db_init('local_file_db','C:\Users\tomask\Documents\Space Physics\Data\mms\data\');

%% LOAD FGM DATA
B3gse = mms.db_get_ts('mms3_fgm_srvy_l2','mms3_fgm_b_gse_srvy_l2',tint);


return;