function MMS_load_mec(date)

global R3;


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
data_dir = ['C:\Users\tomask\Documents\Space Physics\Data\mms\data\mms3\mec\srvy\l2\ephts04d\' yr '\' mo];
%data_dir = ['E:\Data\mms\data\mms1\mec\srvy\l2\ephts04d\' yr '\' mo];
cd(data_dir);

%% LOAD MEC1 DATA
filename = ['mms3_mec_srvy_l2_ephts04d_' yr mo da '*.cdf'];
Mms3_mec_srvy_l2=dataobj(filename);
R3 = mms.variable2ts(get_variable(Mms3_mec_srvy_l2,'mms3_mec_r_gse'));


%% CLEANUP
cd(old_dir);


return;