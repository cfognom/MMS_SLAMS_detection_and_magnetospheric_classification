function [center, duration] = get_center_and_duration(tints)
    if isempty(tints)
        center = [];
        duration = [];
        return
    end
    dif = diff(double(tints.epoch)*1e-9);
    duration = dif(1:2:end);
    center = tints(1:2:end) + duration/2;
end

