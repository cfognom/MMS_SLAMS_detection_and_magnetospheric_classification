function [tints, desc] = read_events(file_path)
%READ_EVENTS Summary of this function goes here
%   Detailed explanation goes here
    data = readcell(file_path, 'Delimiter', ',');
    dates = datetime(data(2:end, 2:3), 'InputFormat', 'yyyy-MM-dd HH:mm:ss Z', 'TimeZone', 'UTC');
    desc = data(2:end, 8);
    starts = datetimeToEpochTT(dates(:,1));
    stops = datetimeToEpochTT(dates(:,2));
    tints = create_tints(starts, stops);
end

