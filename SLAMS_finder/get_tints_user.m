function tints = get_tints_user()
%GET_TINT_BOWSHOCK Summary of this function goes here
%   Detailed explanation goes here
    %% bowshock dates
    % Dates when MMS is near bow shock
    % Taken form MMS science data center orbit plots
    % Phase 2 begins around 2017-03, so the probes might be at the 'edges' of the
    % bow shock onwards from this date
    dates = [
        "2015-10-01", "2016-01-30"; % near bowshock
        "2016-11-01", "2017-05-09"; % near bowshock
        "2017-09-08", "2018-05-23";
        "2018-09-23", "2019-06-16";
        "2019-09-26", "2020-05-14"];

    tints = sort(datetimeToEpochTT(datetime(dates)));
end

