% Calculates precision, recall and F1 of automatically detected SLAMS
% compared to manually detected SLAMS from SLAMS_marker.m

sc = 'MMS1';

plot_comparison = true;

show_train_progress = false;

save_file_name = 'validation_results.txt';

settings = {
    'Include_B_stats', false, ...
    'Include_region_stats', false, ...
    'Include_GSE_coords', false, ...
    'Extra_load_time', 30, ...
    'SLAMS_B_bg_method', 'median', ...
    'SLAMS_B_bg_window', 60, ...
    'SLAMS_threshold', 2, ...
    'SLAMS_min_duration', 0
};

main(sc, settings, save_file_name, plot_comparison, show_train_progress);

function main(sc, settings, save_file_name, plot_comparison, show_train_progress)

    finder = SLAMS_finder('Spacecraft', sc, 'Show_train_progress', show_train_progress);
    
    labeled_files = dir('labeled/data/*.txt');
    rng(100);
    labeled_files = labeled_files(randperm(length(labeled_files)));
    
    all_results = [];
    run_count = 0;
    
    fileID = fopen(save_file_name, 'w');
    for bgmtd = ["median", "mean", "harmmean"]
        settings = setting_set(settings, 'SLAMS_B_bg_method', char(bgmtd));
        for bgtime = [30, 60, 2*60]
            settings = setting_set(settings, 'SLAMS_B_bg_window', bgtime);
            settings = setting_set(settings, 'Extra_load_time', bgtime/2);
            for thresh = [2, 2.5]
                settings = setting_set(settings, 'SLAMS_threshold', thresh);
                for min_dur = [0, 1]
                    settings = setting_set(settings, 'SLAMS_min_duration', min_dur);
                    run_count = run_count + 1;
                    results = validate(finder, labeled_files, settings, 'Plot', plot_comparison);
                    results.bgmtd = bgmtd;
                    results.bgtime = bgtime;
                    results.thresh = thresh;
                    results.min_dur = min_dur;
                    all_results = [all_results, results]; %#ok<AGROW>
                    print_run(fileID, run_count, settings, results)
                end
            end
        end
    end
    fclose('all');
    
    plot_precision_recall(all_results)
end

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
    
    if ~isempty(vertcat(records.SLAMS))
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
        case 60*2
            marker_size = 10;
        case 60
            marker_size = 8;
        case 30
            marker_size = 6;
        end

        switch result.thresh
        case 2.5
            color_face = [1 0 0];
        case 2
            color_face = [0 1 0];
        case 1.8
            color_face = [1 1 0];
        end

        switch result.min_dur
        case 0
            color_face = color_face*0.6;
        case 1
            color_face = color_face*1;
        end

        plot(result.recall, result.precision, 'Marker', symbol, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', color_face, 'MarkerSize', marker_size);
        xlabel('Recall')
        ylabel('Precision')
        axis square
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
    SLAMS = vertcat(records.SLAMS);
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