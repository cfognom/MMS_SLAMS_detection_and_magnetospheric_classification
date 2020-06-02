%% Start
db_path = 'C:\Users\carlh\Documents\MATLAB\Exjobb\Carl\data\';
cache_path = 'cache\';
sc = 'mms1';

tints_search_user = get_tints_user();

tints_active = get_tints_active(tints_search_user);

finder = SLAMS_finder();

settings = {
    'Use_SW_classifier', true, ...
    'SW_smoothing', 3*60, ...
    'ExtraLoadTime', 3*60, ...
    'SLAMS_B_bg_method', 'median' ...
    'SLAMS_B_bg_window', 60 ...
    'SLAMS_Threshold', 2 ...
    'SLAMS_DeltaMaxMerge', 0.5 ...
    'SLAMS_MinDur', 1 ...
    'Compute_B_max', false
};

plot_prediction = true;

n_tints = length(tints_active)/2;

fileID = fopen('identified_SLAMS.txt', 'w');
fprintf(fileID, 'Searching in time intervals: ');
fprintf(fileID, '%s, %s, %s, %s, %s, %s, %s, %s, %s\n', 'id', 'start_UTC', 'end_UTC', 'delta_seconds', 'start_EpochTT', 'end_EpochTT', 'delta_EpochTT', 'B_max', 'B_max_relative');
current_id = 1;
for i = 9:n_tints
% for i = 700:n_tints
    fprintf('Looking for SLAMS in interval %u/%u\n', i, n_tints);
    tint = select_tint(tints_active, i);
    tint = [tint(1) + settings{6}, tint(2) + -settings{6}];
    [tints_SLAMS, B_max, B_max_rel] = finder.evaluate(tint, settings{:});
    if plot_prediction
        finder.plot_prediction(tint);
    end
    current_id = print_SLAMS(fileID, tints_SLAMS, B_max, B_max_rel, current_id);
end
disp('Done!')
fclose(fileID);

function current_id = print_SLAMS(fileID, tints_SLAMS, B_max, B_max_rel, current_id)
    if ~isempty(tints_SLAMS)
        n_SLAMS = length(tints_SLAMS)/2;
        fprintf('Found %u SLAMS!\n', n_SLAMS);
        SLAMS_starts = tints_SLAMS(1:2:end);
        SLAMS_ends = tints_SLAMS(2:2:end);
        delta_epoch = SLAMS_ends.epoch - SLAMS_starts.epoch;
        delta_seconds = double(delta_epoch)*1e-9;
        id = int64((current_id:(current_id + n_SLAMS - 1))');
        current_id = current_id + n_SLAMS;

        cel = cell(n_SLAMS, 9);
        arr = repelem(1, n_SLAMS);
        cel(:, 1) = mat2cell(id, arr, 1);
        cel(:, 2) = mat2cell(SLAMS_starts.utc, arr, 30);
        cel(:, 3) = mat2cell(SLAMS_ends.utc, arr, 30);
        cel(:, 4) = mat2cell(delta_seconds, arr, 1);
        cel(:, 5) = mat2cell(SLAMS_starts.epoch, arr, 1);
        cel(:, 6) = mat2cell(SLAMS_ends.epoch, arr, 1);
        cel(:, 7) = mat2cell(delta_epoch, arr, 1);
        cel(:, 8) = mat2cell(B_max, arr, 1);
        cel(:, 9) = mat2cell(B_max_rel, arr, 1);
        cel = cel';
        fprintf(fileID, '%u, %s, %s, %f, %u, %u, %u, %f, %f\n', cel{:});
    end
end