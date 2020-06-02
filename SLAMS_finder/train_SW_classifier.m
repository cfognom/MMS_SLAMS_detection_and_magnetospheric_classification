function classifier = train_SW_classifier(X)
%TRAIN_SW_MODEL Summary of this function goes here
%   Detailed explanation goes here

    % n_variables = cellfun(@(x) length(x{2}), parameters);
    % id = '';
    % for i = 1:length(parameters)
    %     for v = 1:n_variables(i)
    %         id = [id, '-', char(parameters{i}{2}(v))]; %#ok<AGROW>
    %     end
    % end
    % folder = 'savedState/';
    % fileName_trainData = [folder 'trainData' id '.txt'];
    % fileName_cluster = [folder 'cluster' id '.txt'];
    % fileName_model = [folder 'classifier' id '.mat'];

    % if isfile(fileName_trainData)% && false
    %     X = readmatrix(fileName_trainData);
    % else
    %     rng(100)
    %     n_t = 10000;
    %     t_train = sample_tints(tint, n_t);
    %     X = get_train_data(t_train, parameters);
    %     writematrix(X, fileName_trainData);
    % end


    % figure;
    % plot3(X(:,1), X(:,2), X(:,3), '.');
    % xlabel('abs(b)');
    % ylabel('abs(bulkv_i)');
    % zlabel('n_i')

    % if isfile(fileName_model) && isfile(fileName_cluster) && false
    %     % C = readmatrix(fileName_model);
    %     classifier = loadCompactModel(fileName_model);
    %     y = readmatrix(fileName_cluster);
    % else
        % X = X(:,2:end);

    if false
        X = [X(:,5)./(1 + irf_abs(X(:,6:7), 1)), X(:,8)];
        [n, ~] = size(X);

        % Normalize data
        m = mean(X);
        s = std(X);
        X_norm = (X - m)./s;

        X_trimmed = rmoutliers(X_norm, 'mean');

        figure;
        % c = y == 1:max(y);
        % scatter3(X(:,1), X(:,2), X(:,3), 2, cVec(y,:))
        scatter(X(:,1), X(:,2), 2)
        xlabel('abs(bulkv_i)');
        ylabel('n_i')

        % y = dbscan(X_norm, 0.1, int32(n/60));
        figure;
        ksdensity(X_trimmed, 'BandWidth', 0.05)
        df
        % gm = fitgmdist(X_norm, 3, 'Replicates', 100);
        % y = cluster(gm, X_norm);

        % figure;
        % % cVec = [1 0 0; 0 1 0; 0 0 1];
        % % c = y == 1:max(y);
        % % scatter3(X(:,1), X(:,2), X(:,3), 2, cVec(y,:))
        % scatter(X(:,1), X(:,2), 2, y)
        % ylabel('abs(bulkv_i)');
        % zlabel('n_i')
            
    else
        % X = [irf_abs(X(:,2:4), 1), irf_abs(X(:,5:7), 1), X(:,8)];
        X = [irf_abs(X(:,2:4), 1), X(:,5) - (irf_abs(X(:,6:7), 1)), X(:,8)];
        % X = [irf_abs(X(:,2:4), 1), irf_dot(X(:,5:7), [1, 0, 0]), X(:,8)];
        % X = [irf_abs(X(:,5:7), 1), irf_abs(X(:,11:13), 1), X(:,15)];
        % X = [irf_abs(X(:,2:4), 1), irf_abs(X(:,8:10), 1), X(:,14)];
        [n, ~] = size(X);

        % Normalize data
        m = mean(X);
        s = std(X);
        X_norm = (X - m)./s;

        figure;
        % c = y == 1:max(y);
        % scatter3(X(:,1), X(:,2), X(:,3), 2, cVec(y,:))
        scatter3(X(:,1), X(:,2), X(:,3), 2)
        xlabel('abs(b)');
        ylabel('abs(bulkv_i)');
        zlabel('n_i')

        % Cluster
        % n_clusters = 3;
        % [y, C_norm] = kmeans(X_norm, n_clusters, 'Distance', 'cityblock');

        y = dbscan(X_norm, 0.3, int32(n/60));
        % unique(y)

        % Z = pdist(X_norm);
        % Z = linkage(X_norm, 'ward'); % Great: chebychev, mahalanobis, correlation. Good: spearman.
        % % [coeff,score,latent,tsquared,explained,mu] = pca(X);
        % % Z = linkage(score(:,1:3), 'ward', 'euclidean'); % Great: chebychev, mahalanobis, correlation. Good: spearman.
        % y = cluster(Z, 'Maxclust', 2);

        figure;
        % cVec = [1 0 0; 0 1 0; 0 0 1];
        % c = y == 1:max(y);
        % scatter3(X(:,1), X(:,2), X(:,3), 2, cVec(y,:))
        scatter3(X(:,1), X(:,2), X(:,3), 2, y)
        xlabel('abs(b)');
        ylabel('abs(bulkv_i)');
        zlabel('n_i')

        % dj
        scal = 20;
        model1 = fitcsvm(X(y == 1, :), nonzeros(y == 1), 'Kernelscale', scal, 'Standardize', true, 'KernelFunction', 'rbf');
        model2 = fitcsvm(X(y == 2, :), nonzeros(y == 2), 'Kernelscale', scal, 'Standardize', true, 'KernelFunction', 'rbf');
        model3 = fitcsvm(X(y == -1, :), nonzeros(y == -1), 'Kernelscale', scal, 'Standardize', true, 'KernelFunction', 'rbf', 'OutlierFraction', 0.1);

        classifier = @classifierFunction;

        y = classifier(X);

        figure;
        cVec = [1 0 0; 0 1 0; 0 0 1; 1 0 1];
        % c = y == 1:max(y);
        scatter3(X(:,1), X(:,2), X(:,3), 2, cVec(y,:))
        xlabel('abs(b)');
        ylabel('abs(bulkv_i)');
        zlabel('n_i')
    end

    % dj

    % hold on;
    % scatter3(C(:,1), C(:,2), C(:,3), 20, 'k', 'filled')
    % view(3)

    function [y, score] = classifierFunction(X, smoothWindowSize)
        [~, score1] = predict(model1, X);
        [~, score2] = predict(model2, X);
        [~, score3] = predict(model3, X);
        if nargin > 1 && smoothWindowSize > 1
            score1 = movmean(score1, smoothWindowSize);
            score2 = movmean(score2, smoothWindowSize);
            score3 = movmean(score3, smoothWindowSize);
        end
        score = [score1, score2, score3];
        [~, y] = max(score, [], 2);
    end
end

