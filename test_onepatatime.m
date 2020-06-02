function [label, bg] = test_onepatatime(data, w)
%TEST_ONEPATATIME Summary of this function goes here
%   Detailed explanation goes here
    % data = ts.data;
    n_points = length(data);
    label = zeros(1, n_points);
    bg = zeros(1, n_points);
    bg(1) = data(1);
    % root = 1;
    % state = 0;
    for i = 2:n_points
        bg_m = mean(data(max(i-w-1, 1):i-1));
        % batch = data(i:min(i+50-1, n_points));
        % before_batch = data(max(i-50-1, 1):i-1);
        % a = mean(before_batch)/mean(batch);
        % if max(batch) > ((a-1)*4+1)*bg_m
        if data(i) > 2*bg_m
            % state = state + 1;
            % root = i;
            label(i) = 1;
            bg(i) = 0.005*data(i) + 0.995*bg(i-1);
            % bg(i) = bg(i-1);
        else
            bg(i) = data(i);
            % bg(i) = mean(batch);
        end
        % if data(i) < 0.5*bg_m
        %     root = i;
        %     state = state - 1;
        % end
        % label(i) = state;
    end
end

