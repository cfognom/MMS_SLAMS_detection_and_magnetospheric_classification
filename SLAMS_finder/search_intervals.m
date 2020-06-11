function tints = search_intervals()
    % Time intervals when MMS is near bow shock
    % Taken form MMS science data center orbit plots
    time_intervals = {
        '2015-10-01T00:00:00.000000000Z', '2016-01-30T00:00:00.000000000Z';
        '2016-11-01T00:00:00.000000000Z', '2017-05-09T00:00:00.000000000Z';
        '2017-09-08T00:00:00.000000000Z', '2018-05-23T00:00:00.000000000Z';
        '2018-09-23T00:00:00.000000000Z', '2019-06-16T00:00:00.000000000Z';
        '2019-09-26T00:00:00.000000000Z', '2020-06-08T00:00:00.000000000Z';
    };

    time_intervals = time_intervals';
    time_intervals = reshape(time_intervals, [], 1);
    time_intervals = char(time_intervals);
    tints = EpochTT(time_intervals);
end

