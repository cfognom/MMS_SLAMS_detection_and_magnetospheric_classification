    function output_file_name = MMS_SDC_get_files_TK(start_date,end_date,sc_id,instrument_id,data_rate_mode,data_level,descriptor)


%% Query SDC for mms data
%
% Written:
% 2015-11-16, Per-Arne Lindqvist
% Modified:
% 2015-11-17, Per-Arne Lindqvist
%
%% Run this section once to store username and password in "up_options"
up_options = weboptions('Username','tkarlsson','Password','Mynta_012');
% disp(up_options);
%
%% Run this section to setup what to ask for and build query
%
% enter comma-separated list of spacecraft (zero or more)
% (mms1, mms2, mms3, mms4)
% sc_id = 'mms1,mms2,mms3,mms4';
%
% enter comma-separated list of instruments
% (afg, asp1, asp2, aspoc, dfg, dsp, edi, edp, epd-eis, feeps, fields,
% fpi, hpca, mec, scm, sdp)
% instrument_id = 'fpi';
%
% enter comma-separated list of data rate mode (normally only one)
% (brst, comm, f128, fast, hk, slow, srvy)
% data_rate_mode = 'fast';
%
% enter comma-separated list of data level (normally only one)
% (l1a, l1b, l2, l2pre, ql, sitl)
% data_level = 'l2';
%
% enter comma-separated list of descriptors (zero or more)
% (101, 105, 10e, 14f, 173, 174, 175, 176, 177, 178, 179, 17a, 17b, 2015,
% ace, amb, beam, bpsd, cal, dce, dce128, dce2d, dce32, dce64, dce8,
% dcecomm, dcv, dcv128, dcv32, dcv64, dcvcomm, des-dist, des-moms,
% dis-dist, dis-moms, efield, electron, electron-bottom, electron-top,
% electronenergy, epht89d, epht89q, ephts04d, epsd, extof, hmfe, ion,
% ion-bottom, ion-top, logicals, moments, partenergy, phxtof, plot, sc128,
% sc256, sc32, scb, scf, schb, scm, scmcomm, scpot, scs, stat, swd,
% sweeps, tdn)
% descriptor = 'des-moms,dis-moms';
%
% enter start date in form yyyy-mm-dd
%start_date = '2017-05-28';
%
% enter end date in form yyyy-mm-dd
%end_date = '2017-05-29';
%
% build query
query = '';
separator = '?';
if length(sc_id) > 0
  query = [query separator 'sc_id=' sc_id];
  separator = '&';
end
if length(instrument_id) > 0
  query = [query separator 'instrument_id=' instrument_id];
  separator = '&';
end
if length(data_rate_mode) > 0
  query = [query separator 'data_rate_mode=' data_rate_mode];
  separator = '&';
end
if length(data_level) > 0
  query = [query separator 'data_level=' data_level];
  separator = '&';
end
if length(descriptor) > 0
  query = [query separator 'descriptor=' descriptor];
  separator = '&';
end
if length(start_date) > 0
  query = [query separator 'start_date=' start_date];
  separator = '&';
end
if length(end_date) > 0
  query = [query separator 'end_date=' end_date];
  separator = '&';
end
%
%% Run this section to get comma-separated list of available file names
%
url_init = 'https://lasp.colorado.edu/mms/sdc/sitl/files/api/v1/file_names/science';
%
% execute the query and put list of file names in variable "file_names"
file_names = webread([url_init query], up_options);
% "file_names" now contains a comma-separated string of all found files
disp(file_names)
%
%% Run this section to get data
%
% Select output file directory and name
%
% If an output directory is given, it must exist
output_directory = 'E:\Data\';
%output_directory = 'somedirectory/';
%
% If only one file is requested, the output file name should be <name>.cdf
% where <name> is the same as obtained by the file name query in the
% previous section.
% If multiple files are requested, the output file name should be
% <name>.zip, where <name> is arbitrary. Unzip this file to extract the
% data files to the mms standard hierarchical directory structure.
%output_file = ['download_ions_' start_date '--'  end_date(6:10) '.cdf'];
output_file = ['download_' instrument_id '_' start_date '--'  end_date(6:10) '.zip'];
%
url_init = 'https://lasp.colorado.edu/mms/sdc/sitl/files/api/v1/download/science';
%
% execute the query and receive the output in file "output directory""output_file"
output_file_name = websave([output_directory output_file], [url_init query], up_options);
% "output_file_name" now contains the full path to the output file
disp(output_file_name);
%
%
%% Unpack data file
%




