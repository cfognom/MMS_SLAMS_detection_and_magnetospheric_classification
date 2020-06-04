function tints_SLAMS = testfunc(tints_SLAMS, tint_allowed)
%TESTFUNC Summary of this function goes here
%   Detailed explanation goes here
    % if isempty(tints_SLAMS_SW)
    %     tints_SLAMS = [];
    %     return
    % end
    % mem = ismember(tints_SLAMS_SW.epoch, tints_SLAMS.epoch);
    % dif = diff(mem);
    % idx = dif(1:2:end) == 0;
    % idx = reshape([idx, idx]', [], 1);
    % tints_SLAMS = tints_SLAMS_SW(idx);

    start_idx = find(tints_SLAMS >= tint_allowed(1), 1);
    if mod(start_idx, 2) == 0
        start_idx = start_idx + 1;
    end
    stop_idx = find(tints_SLAMS <= tint_allowed(2), 1, 'last');
    if mod(stop_idx, 2) == 1
        stop_idx = stop_idx - 1;
    end
    tints_SLAMS = tints_SLAMS(start_idx:stop_idx);
end

