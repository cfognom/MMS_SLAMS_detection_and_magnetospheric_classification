function SLAMS = filter_SLAMS(SLAMS, varargin)
    %Filter_SLAMS filters SLAMS using a the following name value pairs:
    %   'ID', [Numeric array of ids]
    %   'After_date', [UTC char array]
    %   'Before_date', [UTC char array]
    %   'Min_duration', [Numeric value in seconds]
    %   'Max_duration', [Numeric value in seconds]
    %   'GSE_filter', [Handle to a function with signature:
    %       bool = myFunc(GSE),
    %       where bool is a column vector with the same length as the number of SLAMS and
    %       GSE is a [n_SLAMS x 3] matrix (columns are x,y,z in GSE).
    %       If bool is true, SLAMS is included]
    %   'B_filter', [Handle to a function with signature:
    %       bool = myFunc(B_background_mean, B_mean, B_max, B_relative_mean, B_relative_max),
    %       where bool is a column vector with the same length as the number of SLAMS and
    %       inputs are also column vectors with the same size.
    %       If bool is true, SLAMS is included]
    %   'Region_filter', [Handle to a function with signature:
    %       bool = myFunc(region_posterior, region_mahaldist, region_logpdf),
    %       where bool is a column vector with the same length as the number of SLAMS and
    %       region_posterior and region_mahaldist is [n_SLAMS x n_classes x (1 + n_windows)] matrices.
    %       The first dimension is for each SLAMS entry, the second is for each region class and
    %       the third is for the instant values (first index) and the mean values of the region time windows (second to end indices).
    %       region_logpdf is a [n_SLAMS x 1 x (1 + n_windows)] matrix of the log of the probability density function.
    %       The first dimension separates the data for each SLAMS and
    %       the third separates instant values (first index) and mean values over the region time windows (second to end indices).
    %       If bool is true, SLAMS is included]

    p = inputParser;
    addRequired(p, 'SLAMS');
    addParameter(p, 'ID', [])
    addParameter(p, 'After_date', [])
    addParameter(p, 'Before_date', [])
    addParameter(p, 'Min_duration', [])
    addParameter(p, 'Max_duration', [])
    addParameter(p, 'Region_class', [])
    addParameter(p, 'GSE_filter', [])
    addParameter(p, 'B_filter', [])
    addParameter(p, 'Region_filter', [])
    parse(p, SLAMS, varargin{:})
    r = p.Results;

    if ~isempty(r.ID)
        idx = ismember([SLAMS.id], r.ID);
        SLAMS = SLAMS(idx);
    end
    
    if ~isempty(r.After_date)
        idx = [SLAMS.start] > EpochTT(r.After_date);
        SLAMS = SLAMS(idx);
    end

    if ~isempty(r.Before_date)
        idx = [SLAMS.stop] < EpochTT(r.Before_date);
        SLAMS = SLAMS(idx);
    end

    if ~isempty(r.Min_duration) || ~isempty(r.Max_duration)
        dur = [SLAMS.stop] - [SLAMS.start];
        if ~isempty(r.Min_duration)
            idx = dur > r.Min_duration;
            SLAMS = SLAMS(idx);
        end
        if ~isempty(r.Max_duration)
            idx = dur < r.Max_duration;
            SLAMS = SLAMS(idx);
        end
    end

    if ~isempty(r.GSE_filter)
        idx = r.GSE_filter(vertcat(SLAMS.pos_GSE));
        SLAMS = SLAMS(logical(idx));
    end

    if ~isempty(r.B_filter)
        idx = r.B_filter( ...
            vertcat(SLAMS.B_bg_mean), ...
            vertcat(SLAMS.B_mean), ...
            vertcat(SLAMS.B_max), ...
            vertcat(SLAMS.B_rel_mean), ...
            vertcat(SLAMS.B_rel_max) ...
        );
        SLAMS(logical(idx));
    end

    if ~isempty(r.Region_filter)
        posteriors = permute(cat(3, SLAMS.region_posterior), [3, 2, 1]);
        mahaldists = permute(cat(3, SLAMS.region_mahaldist), [3, 2, 1]);
        logpdfs = permute(cat(3, SLAMS.region_logpdf), [3, 2, 1]);

        idx = r.Region_filter(posteriors, mahaldists, logpdfs);
        % size(find(idx))
        SLAMS = SLAMS(logical(idx));
    end
end

