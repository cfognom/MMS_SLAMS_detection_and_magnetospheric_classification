function MMS_SDC_get_B_from_lists(listpath)

filelist = dir('*.mat');
[nfiles x] = size(filelist);

%for n=1:nfiles
for n=1:1

    readstruct = load([filelist(n).folder '\' filelist(n).name]);
    jetlist = readstruct.jetlist;
    startutc = jetlist(1).startutc;
    yr = 
    
end;



return;