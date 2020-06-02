function [tints_classes, score_TS] = eval_SW_classifier(tints, parameters, classifier)
%CLASSIFY_SW Summary of this function goes here
%   Detailed explanation goes here

    % Init
    n_tints = floor(length(tints)/2);
    n_instruments = length(parameters);
    n_variables = cellfun(@(x) length(x{2}), parameters);
    n_variables_tot = sum(n_variables);
    classifier_time_scale = 10; % seconds
    [n_classes, ~] = size(classifier);
    n_classes = 3;
    % T = delaunayn(C);

    tints_classes = cell(n_classes, 1);
    for c = 1:n_classes
        tints_classes{c} = [];
    end
    empty_tints_idx = zeros(1, n_tints*2);
    for it = 1:n_tints
        fprintf('Classifying interval %u/%u\n', it, n_tints);
        tint = tints(it*2 - 1:it*2);
        n_points = ceil((tint(2) - tint(1))/classifier_time_scale);
        c = 1;
        t = linspace(double(tint(1).epoch), double(tint(2).epoch), n_points);
        X = zeros(n_points, n_variables_tot);
        for i = 1:n_instruments
            p = parameters{i};
            try
                TS = MMS_load(p{1}, p{2}, tint);
                for v = 1:n_variables(i)
                    if iscell(TS)
                        ts = TS{v};
                    else
                        ts = TS;
                    end
                    if length(p) > 2 && p{3}(v) == "abs"
                        ts = irf_abs(ts);
                    end
                    [~, d] = size(ts.data);
                    X(:,c:c + d - 1) = interp1(double(ts.time.epoch), ts.data, t, 'linear', 'extrap');
                    c = c + d;
                end 
            catch
                empty_tints_idx(it*2 - 1:it*2) = 1;
            end
        end
        % Remove bad points, i.e. points that are missing some instrument data.
        % X = X(bad_points == 0, :);

        % y = dsearchn(classifier, X);
        if empty_tints_idx(it*2) == 0

            % [y, score] = predict(classifier, X);
            % y = classifier(X);
            % [~, score] = predict(classifier, X);
            % X = movmean(X, 1001, 1);
            [y, Score] = classifier(X);
            % score = movmean(score, 21);
            % y = (score < 0) + 1;
            if it == 1
                score_TS = TSeries(EpochTT(int64(t)), Score);
            else
                score_TS = score_TS.combine(TSeries(EpochTT(int64(t)), Score));
            end
            % y = classifier(X);
            % y = movfun(@majorityVote, y, 11);
            
            % X
            % diff(y)
            id = find(diff(y))';
            % id
            % size(X)
            % idx
            % dates = t(idx);
            d = (t(id) + t(id + 1))/2;
            dates = [t(1), reshape([d; d], 1, []), t(n_points)];
            idx = [1, reshape([id; id + 1], 1, []), n_points];
            clas = y(idx);
            for c = 1:n_classes
                EpochTT(int64(dates(clas == c)));
                tints_classes{c} = horzcat(tints_classes{c}, EpochTT(int64(dates(clas == c))));
            end

            % b = irf_abs(MMS_Load("mms1_fgm_srvy_l2", "b_gse", tint));
            % tmp = MMS_Load("mms1_fpi_fast_l2_dis-moms", ["bulkv_gse", "numberdensity"], tint);
            % v = irf_abs(tmp{1});
            % n = tmp{2};
            % X = [interp1(b.time.epoch, b.data, t), ...
            %     interp1(v.time.epoch, v.data, t), ...
            %     interp1(n.time.epoch, n.data, t)];
            % point2cluster = dsearchn(classifier, X);
        end
    end
    
    function v = majorityVote(x)
        [~, v] = max(accumarray(x',1));
    end
end

