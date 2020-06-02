function MMS_load_edp_fast_l2_sc1(starttime,stoptime)

global E1gse;

%% FIX INPUT TIMES
startepochtt = irf_time(starttime,'vector>epochtt');
stopepochtt = irf_time(stoptime,'vector>epochtt');
tint  = irf.tint(startepochtt,stopepochtt);

%% INITIATE DATABASE
%mms.db_init('local_file_db','D:\Data\mms\data\');
%mms.db_init('local_file_db','C:\Users\tomask\Documents\Space Physics\Data\mms\data\');
mms.db_init('local_file_db','E:\Data\mms\data\');

%% LOAD FGM DATA
E1gse = mms.db_get_ts('mms1_edp_fast_l2_dce','mms1_edp_dce_gse_fast_l2',tint);


return;