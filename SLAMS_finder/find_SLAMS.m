function [bg_ts, limit_ts, time] = find_SLAMS(tints)
%FIND_SLAMS Summary of this function goes here
%   Detailed explanation goes here

    [starts, stops, n_tints] = unpack_tints(tints);

    global b_smooth SLAMS_time limit_ts
    n_tints = length(tints_SW)/2
    SLAMS_time = [];
    for i = 1:n_tints
        tint = tints_SW(i*2 - 1:i*2);
        ts = irf_abs(MMS_load("mms1_fgm_srvy_l2", "b_gse", tint));
        [b_smooth_new, limit_ts_new, SLAMS_time_new] = find_SLAMS(ts);
        if i == 1
            b_smooth = b_smooth_new;
            limit_ts = limit_ts_new;
        else
            b_smooth = b_smooth.combine(b_smooth_new); 
            limit_ts = limit_ts.combine(limit_ts_new); 
        end
        size(SLAMS_time)
        size(SLAMS_time_new)
        SLAMS_time = [SLAMS_time; SLAMS_time_new];
    end
    % SLAMS_time = SLAMS_time'
    fprintf('Detected %u SLAMS', length(SLAMS_time)/2)

    data = ts.data;
    dif_ts = diff_ts(ts);
    variance = movvar(data, 1001);
    % data1 = [ts.data, abs(dif_ts.data)];
    data1 = [data, variance];

    % dif = diff(data)
    % m = data./movmean(data, 1001)
    % v = movvar(data, 1001)
    % scatter3(data(1:end-1), dif, v(1:end-1), 2)
    % xlabel('abs(b)');
    % ylabel('diff');
    % zlabel('variance')

    % data = ts;
    % bg_data = movmean(data, 1001);
    % bg_data = movmedian(data, 5001);
    % bg_data = movfun(@(x) bg_ransac(x), data, 1001);
    % bg_data = movfun(@(x) bg_LAR(x), data, 1001);
    % bg_data = movfun(@bg_weighted_mean, data1, 1, 1001);
    d = movmean(data, 101);
    mi = movmin(d, 501);
    c = 0.1;
    w = 1./(1 + c*(d - mi));
    bg_data = d.*w;
    % bg_data = movmean(bg_data, 501);

    % bg_data = variance
    % bg_data = movmean(abs(dif_ts.data), 101)
    % size(bg_data)
    % size(ts.time)
    bg_ts = TSeries(ts.time, bg_data);
    coeff = 2.5;
    % bg_ts = bg_ransac_TS(ts);
    limit_ts = TSeries(bg_ts.time, bg_data*coeff);
    % bg_ts = bg_LAR_TS(ts);

    rel = data./bg_ts.data;
    dif = diff(rel >= coeff);
    idx = find(dif);
    if ~isempty(idx)
        if dif(idx(1)) < 0
            idx = [1; idx];
        end
        if dif(idx(end)) > 0
            idx = [idx; length(data)];
        end
    end
    % time = arrayfun(@(x) irf_time(x, 'epochtt>epoch'), ts.time(idx));
    time = ts.time.epochUnix(idx);
    % time = reshape(time, [],2);

    % olier = isoutlier(data, 'mean');
    % nnz(olier)
    % time = arrayfun(@(x) irf_time(x, 'epochtt>epoch'), ts.time(olier == 1));

    function out = bg_weighted_mean(in)
        x = in(:, 1);
        v = in(:, 2);
        % d = abs(diff_data(x, 1));
        coeff = 4;
        w = 1./(1 + coeff*v);
        out = sum(w.*x)/sum(w);
        out = sum(x./(1 + v))/length(x);
        % out = sum(w.*x)/length(x);
    end

    function out = bg_ransac(y)
        n = length(y);
        x = 1:n;
        % ii = (1:numel(idx))';
        % a = accumarray([rem(ii-1,n)+1,ceil(ii/n)],idx,[],[],nan);
        % n = length(a)
        % a(1,1)
            points = double([x; y]);
            % sampleSize = 2; % number of points to sample per trial
            % maxDistance = 2; % max allowable distance for inliers
            
            % fitLineFcn = @(points) polyfit(points(:,1),points(:,2),1); % fit function using polyfit
            % evalLineFcn = ...   % distance evaluation function
            %   @(model, points) sum((points(:, 2) - polyval(model, points(:,1))).^2,2);
            
            % [modelRANSAC, inlierIdx] = ransac(points,fitLineFcn,evalLineFcn, ...
            %   sampleSize,maxDistance);

            [bestParameter1, bestParameter2] = ransac(points, 2, 1000, 0.3, 0);
            out = double(bestParameter1*ceil(n/2) + bestParameter2);
    end

    function out = bg_LAR(yx)
        n = length(yx);
        x = 1:n;

        y = @(b,x) b(1)*x + b(2);             % Objective function

        OLS = @(b) sum(abs(y(b,x) - yx));          % Ordinary Least Squares cost function
        opts = optimset('MaxFunEvals',50000, 'MaxIter',10000);
        B = fminsearch(OLS, rand(2,1), opts);       % Use ‘fminsearch’ to minimise the ‘OLS’ function

        out = y(B,ceil(n/2));
    end

    function out = bg_LAR_TS(ts)
        % n = length(y);
        yx = double(ts.data)';
        x = double(ts.time.epoch)';

        y = @(b,x) b(1)*x + b(2);             % Objective function

        OLS = @(b) sum(abs(y(b,x) - yx));          % Ordinary Least Squares cost function
        opts = optimset('MaxFunEvals',50000, 'MaxIter',10000);
        B = fminsearch(OLS, rand(2,1), opts)       % Use ‘fminsearch’ to minimise the ‘OLS’ function

        out = TSeries(ts.time, y(B,x)');
    end

    function out = bg_ransac_TS(ts)
        % n = length(y);
        y = ts.data';
        x = ts.time.epoch';
        points = double([x; y]);

        [bestParameter1, bestParameter2] = ransac(points, 2, 1000, 0.3, 0);
        y = double(bestParameter1*x + bestParameter2);
        out = TSeries(ts.time, y');
        % X = [ts.time.epoch, ts.data];
        % fitFcn = @(x) polyfit(x(:, 1), x(:, 2), 1);
        % distFcn = @(model, x) sum((x(:, 2) - polyval(model, x(:, 1))).^2, 2);
        % [model, inlierIdx] = ransac(X, fitFcn, distFcn, 2, 0.3);
        % disp("hej")
        % data = polyval(model, X(:, 1));
        % out = TSeries(ts.time, data);
    end
end

