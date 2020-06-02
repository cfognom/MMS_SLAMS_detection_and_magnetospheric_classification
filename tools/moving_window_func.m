function ts = moving_window_func(ts, time_window, func, sample_points)
%MOVING_WINDOW_FUNC Summary of this function goes here
%   Detailed explanation goes here
    if nargin > 3
        n_points = length(sample_points);
        points = sample_points;
    else
        n_points = length(ts.time);
        points = ts.time;
    end
    data = zeros(n_points, 1);
    dt = time_window/2;
    for i = 1:n_points
        tint = [points(i) + -dt, points(i) + dt];
        data(i) = func(ts.tlim(tint));
    end
    ts = TSeries(points, data);
end

