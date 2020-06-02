function diffed_ts = diff_ts(ts)
%DIFF_TS Summary of this function goes here
%   Detailed explanation goes here
    dd = diff_data(ts.data, 1);
    dt = diff_data(double(ts.time.epoch)*1e-9, 1);
    diffed_ts = TSeries(ts.time, dd./dt);
end

