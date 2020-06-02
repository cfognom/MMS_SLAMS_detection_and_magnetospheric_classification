function MMS_load_fpi_l2_moms_sc1(starttime,stoptime)

global vi1gse ve1gse ni1 ne1 Tiperp1 Tipar1 Teperp1 Tepar1 pressuretensori1gse pressuretensore1gse Eispectomni1 Eespectomni1 Ei1 Ee1;


%% FIX INPUT TIMES
startepochtt = irf_time(starttime,'vector>epochtt');
stopepochtt = irf_time(stoptime,'vector>epochtt');
tint  = irf.tint(startepochtt,stopepochtt);

%% INITIATE DATABASE
mms.db_init('local_file_db','C:\Users\carlh\Documents\MATLAB\Exjobb\Carl\data\');
%mms.db_init('local_file_db','D:\Data\mms\data\');
%mms.db_init('local_file_db','F:\Data\mms\data\');
%mms.db_init('local_file_db','C:\Users\tomask\Documents\Space Physics\Data\mms\data\');

%% LOAD FPI1 DATA
pressuretensori1gse = mms.db_get_ts('mms1_fpi_fast_l2_dis-moms','mms1_dis_prestensor_gse_fast',tint);
pressuretensore1gse = mms.db_get_ts('mms1_fpi_fast_l2_des-moms','mms1_des_prestensor_gse_fast',tint);
vi1gse = mms.db_get_ts('mms1_fpi_fast_l2_dis-moms','mms1_dis_bulkv_gse_fast',tint);
ve1gse = mms.db_get_ts('mms1_fpi_fast_l2_des-moms','mms1_des_bulkv_gse_fast',tint);
ni1 = mms.db_get_ts('mms1_fpi_fast_l2_dis-moms','mms1_dis_numberdensity_fast',tint);
ne1 = mms.db_get_ts('mms1_fpi_fast_l2_des-moms','mms1_des_numberdensity_fast',tint);
Tipar1 = mms.db_get_ts('mms1_fpi_fast_l2_dis-moms','mms1_dis_temppara_fast',tint);
Tiperp1 = mms.db_get_ts('mms1_fpi_fast_l2_dis-moms','mms1_dis_tempperp_fast',tint);
Tepar1 = mms.db_get_ts('mms1_fpi_fast_l2_des-moms','mms1_des_temppara_fast',tint);
Teperp1 = mms.db_get_ts('mms1_fpi_fast_l2_des-moms','mms1_des_tempperp_fast',tint);
Eispectomni1 =  mms.db_get_ts('mms1_fpi_fast_l2_dis-moms','mms1_dis_energyspectr_omni_fast',tint);
Eespectomni1 =  mms.db_get_ts('mms1_fpi_fast_l2_des-moms','mms1_des_energyspectr_omni_fast',tint);
Ei1 =  mms.db_get_ts('mms1_fpi_fast_l2_dis-moms','mms1_dis_energy_fast',tint);
Ee1 =  mms.db_get_ts('mms1_fpi_fast_l2_des-moms','mms1_des_energy_fast',tint);

return;