function MMS_SDC_getalldata_SC1(start_date_string,stop_date_string)

%Format for start and stop date strings: '2017-04-01'

MMS_SDC_get_files_TK(start_date_string,stop_date_string,'mms1','fgm','srvy','l2','');
MMS_SDC_get_files_TK(start_date_string,stop_date_string,'mms1','fpi','fast','l2','');
MMS_SDC_get_files_TK(start_date_string,stop_date_string,'mms1','mec','srvy','l2','ephts04d');

return;