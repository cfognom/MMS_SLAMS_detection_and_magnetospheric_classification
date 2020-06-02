function tints_active = get_tints_instrument_active(filePrefix, tints, time_threshold)
%GET_TINTS_INSTRUMENT_ACTIVE Summary of this function goes here
%   Detailed explanation goes here

    n_tints = length(tints)/2;

    starts = [];
    stops = [];
    for i = 1:n_tints
        fprintf('Getting active time intervals for ''%s'' %u/%u\n', filePrefix, i, n_tints)
        tint = select_tint(tints, i);
        files = mms.db_list_files(filePrefix, tint);
        if ~isempty(files)
            starts = [starts; files.start]; %#ok<AGROW>
            stops = [stops; files.stop]; %#ok<AGROW>
        end
    end
    if isempty(starts)
        tints_active = [];
        return
    end
    
    n_tints_active = length(starts);

    tints_active = create_tints(starts, stops);

    dif = starts(2:end) - stops(1:end - 1);
    % a = sort(dif);
    % format long
    % a(end - 100:end)
    idx = 2*find(dif > time_threshold)';
    idx = [1, reshape([idx; idx + 1], 1, []), 2*n_tints_active];
    tints_active = tints_active(idx);
end

