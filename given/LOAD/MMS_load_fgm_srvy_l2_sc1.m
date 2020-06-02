function MMS_load_fgm_srvy_l2_sc1(starttime,stoptime)

global B1gse;

%% FIX INPUT TIMES
startepochtt = irf_time(starttime,'vector>epochtt');
stopepochtt = irf_time(stoptime,'vector>epochtt');
tint  = irf.tint(startepochtt,stopepochtt);

%% INITIATE DATABASE
mms.db_init('local_file_db','C:\Users\carlh\Documents\MATLAB\Exjobb\data\');
%mms.db_init('local_file_db','D:\Data\mms\data\');


%% LOAD FGM DATA
B1gse = mms.db_get_ts('mms1_fgm_srvy_l2','mms1_fgm_b_gse_srvy_l2',tint);


return;