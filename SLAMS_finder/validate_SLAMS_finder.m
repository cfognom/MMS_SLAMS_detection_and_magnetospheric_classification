%% Start
db_path = 'C:\Users\carlh\Documents\MATLAB\Exjobb\Carl\data\';
cache_path = 'cache\';
sc = 'mms1';

tints_search_user = get_tints_user();

tints_active = get_tints_active(tints_search_user);

finder = SLAMS_finder('Show_train_progress', true);

labeled_files = dir('labeled/data/*.txt');
rng(100);
labeled_files = labeled_files(randperm(length(labeled_files)));
% labeled_files = labeled_files(1:2);

settings = {
    'Include_B_stats', true, ...
    'Include_region_stats', true, ...
    'Region_time_windows', [15, 30, 60, 2*60, 4*60, 8*60], ...
    'Include_GSE_coords', true, ...
    'Extra_load_time', 4*60, ...
    'SLAMS_B_bg_method', 'median', ...
    'SLAMS_B_bg_window', 60, ...
    'SLAMS_threshold', 2
};

all_results = [];
run_count = 0;

fileID = fopen('validation_results.txt', 'w');
for bgmtd = ["median", "mean", "harmmean"]
    settings = setting_setter(settings, 'SLAMS_B_bg_method', char(bgmtd));
    for bgtime = [60, 3*60, 10*60]
        settings = setting_setter(settings, 'SLAMS_B_bg_window', bgtime);
        % settings{6} = max(bgtime, settings{6});
        for thresh = [1.8, 2, 2.5]
            run_count = run_count + 1;
            settings = setting_setter(settings, 'SLAMS_threshold', thresh);
            results = validate(finder, labeled_files, settings, 'Plot', true);
            results.bgmtd = bgmtd;
            results.bgtime = bgtime;
            results.thresh = thresh;
            all_results = [all_results, results]; %#ok<AGROW>
            print_run(fileID, run_count, settings, results)
        end
    end
end
fclose(fileID);

plot_precision_recall(all_results)

function results = validate(finder, labeled_files, settings, varargin)

    p = inputParser;
    addRequired(p, 'finder');
    addRequired(p, 'labeled_files');
    addRequired(p, 'settings');
    addParameter(p, 'Plot', false);
    parse(p, finder, labeled_files, settings, varargin{:});
    r = p.Results;

    n_train_files = length(r.labeled_files);
    records = cell(n_train_files, 1);
    for i = 1:n_train_files
        fprintf('Validating interval %u/%u\n', i, n_train_files);
        labeled = load_labeled(r.labeled_files(i));
        records{i} = labeled;
        tint = labeled.tint;

        SLAMS = r.finder.evaluate(tint, settings{:});
        % if ~isempty(SLAMS)
        %     SLAMS.region_posterior
        %     SLAMS.region_posterior_windows
        %     df
        % end
        if r.Plot
            r.finder.plot_comparison(tint, labeled.tints_true_SLAMS);
        end
        records{i}.SLAMS = SLAMS;
    end
    records = vertcat(records{:});
    
    if ~isempty([records.SLAMS])
        [precision, recall, F1] = calc_performance(records);
        results.precision = precision;
        results.recall = recall;
        results.F1 = F1;
    else
        results.precision = nan;
        results.recall = nan;
        results.F1 = nan;
    end
end

function print_run(fileID, run_count, settings, results)
    fprintf(fileID, '--- run %u -----------------------------------------------', run_count);

    settings_str = construct_settings_string(settings);
    settings_str = vertcat({'Settings:'}, settings_str);
    settings_str = [strjoin(settings_str, '\n\t'), '\n'];
    fprintf(fileID, settings_str);

    results_str = {
        'Precision', mat2str(results.precision);
        'Recall', mat2str(results.recall);
        'F1', mat2str(results.F1)
    };
    results_str = join(results_str, ' = ');
    results_str = [strjoin(results_str, '\n\t'), '\n'];
    fprintf(fileID, results_str);
end

function plot_precision_recall(results)
    tick = linspace(0, 1);
    [X, Y] = meshgrid(tick, tick);
    Z = 2*(X.*Y)./(X + Y);
    figure;
    contour(X,Y,Z, 'ShowText', 'on')
    hold on
    grid on
    for result = results

        switch result.bgmtd
        case 'mean'
            symbol = 's';
        case 'median'
            symbol = 'd';
        case 'harmmean'
            symbol = 'o';
        end

        switch result.bgtime
        case 60*10
            marker_size = 10;
        case 60*3
            marker_size = 8;
        case 60
            marker_size = 6;
        end

        switch result.thresh
        case 2.5
            color_face = 'r';
        case 2
            color_face = 'g';
        case 1.8
            color_face = 'y';
        end

        plot(result.recall, result.precision, 'Marker', symbol, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', color_face, 'MarkerSize', marker_size);
        xlabel('Recall')
        ylabel('Precision')
    end
end

function labeled = load_labeled(file)
    
    file_name = file.name;
    tints_SLAMS = EpochTT(int64(readmatrix([file.folder, '\', file_name])));
    file_name = split(file_name, '.');
    file_name = split(file_name{1}, '_');
    epoch_start = int64(str2double(file_name{2}));
    epoch_end = int64(str2double(file_name{3}));
    tint = EpochTT([epoch_start, epoch_end]);
    
    tints_SLAMS = intersect_tints(tint, tints_SLAMS);
    if isempty(tints_SLAMS)
        tints_SLAMS = [];
    end

    labeled.tint = tint;
    labeled.tints_true_SLAMS = tints_SLAMS;
end

function [precision, recall, F1] = calc_performance(records)
    SLAMS = [records.SLAMS];
    tints_true_SLAMS = [records.tints_true_SLAMS];
    tints_pred_SLAMS = sort([[SLAMS.start]; [SLAMS.stop]]);
    [center_true_SLAMS, duration_true_SLAMS] = get_center_and_duration(tints_true_SLAMS);
    [center_pred_SLAMS, duration_pred_SLAMS] = get_center_and_duration(tints_pred_SLAMS);
    center_error_matrix = 1e-9*double(abs(center_true_SLAMS.epoch - center_pred_SLAMS.epoch'));
    acceptable_duration_offset = (0.5*duration_true_SLAMS + 0.5*duration_pred_SLAMS');
    touching = center_error_matrix < acceptable_duration_offset;
    recall_arr = max(touching, [], 2);
    recall = nnz(recall_arr)/length(recall_arr);
    % true_centers_inside_pred = center_error_matrix < 0.5*duration_pred_SLAMS';
    precision_arr = max(touching, [], 1);
    precision = nnz(precision_arr)/length(precision_arr);
    F1 = 2*(recall*precision)/(recall + precision);
end