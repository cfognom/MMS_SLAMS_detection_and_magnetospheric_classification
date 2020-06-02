function MMS_load_mec_srvy_l2_ephts04d_sc1(date)

global R1;


%% MOVE TO CORRECT DIRECTORY
old_dir = cd;
yr = int2str(date(1));
mo = date(2);
if mo < 10
    mo = ['0' int2str(mo)];
else
    mo = int2str(mo);
end;
da = date(3);
if da < 10
   da = ['0' int2str(da)];
else
    da = int2str(da);
end;


%% GOTO MEC1 DATA DIR
data_dir = ['C:\Users\carlh\Documents\MATLAB\Exjobb\Carl\data\mms1\mec\srvy\l2\ephts04d\' yr '\' mo];
mms.db_init('local_file_db','C:\Users\carlh\Documents\MATLAB\Exjobb\Carl\data\');


cd(data_dir);

%% LOAD MEC1 DATA
filename = ['mms1_mec_srvy_l2_ephts04d_' yr mo da '*.cdf'];
Mms1_mec_srvy_l2=dataobj(filename);
R1 = mms.variable2ts(get_variable(Mms1_mec_srvy_l2,'mms1_mec_r_gse'));


%% CLEANUP
cd(old_dir);


return;