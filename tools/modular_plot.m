classdef modular_plot < handle
    %MODULAR_PLOT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        fig
        original_zoom
        plot_title
        plot_axes
        hcas
        current_hca_key
        execute_list
    end
    
    methods
        function obj = modular_plot(varargin)
            %MODULAR_PLOT Construct an instance of this class
            %   Detailed explanation goes here
            p = inputParser;
            addParameter(p, 'title', 'default title', @ischar);
            parse(p, varargin{:})
            r = p.Results;

            obj.plot_title = r.title;
            obj.execute_list = {};
            obj.hcas = containers.Map();
        end

        function hca_func = get_hca_func(obj, name)
            if isempty(name)
                key = obj.current_hca_key;
                hca_func = @() pick_hca(key);
            else
                if ~isKey(obj.hcas, name)
                    % hca = irf_panel(name);
                    obj.current_hca_key = name;
                    obj.hcas(name) = [];
                    hca_func = @() add_hca(name);
                else
                    obj.current_hca_key = name;
                    hca_func = @() pick_hca(name);
                end
            end

            function hca = add_hca(name)
                hca = irf_panel(name);
                obj.hcas(name) = hca;
            end

            function hca = pick_hca(name)
                hca = obj.hcas(name);
            end
        end
        
        function lineplot(obj, name, ts, varargin)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            p = inputParser;
            addRequired(p, 'name', @ischar);
            addRequired(p, 'ts');
            addParameter(p, 'colorOrder', [])
            addParameter(p, 'color', {})
            addParameter(p, 'fontSize', 16)
            addParameter(p, 'ylabel', 'default ylabel')
            addParameter(p, 'legend', [])
            parse(p, name, ts, varargin{:})
            r = p.Results;
            
            hca_func = obj.get_hca_func(r.name);
            obj.execute_list{end + 1} = @execute_lineplot;
            
            function execute_lineplot()
                hca = hca_func();
                if ~isempty(r.colorOrder)
                    hca.set('ColorOrder', r.colorOrder);
                end
                if ~isempty(r.color)
                    r.color = {'color', r.color};
                end

                hold(hca,'on');
                if isa(r.ts, 'cell')
                    n_data = length(r.ts);
                    for i = 1:n_data
                        irf_plot(hca, r.ts{i}, r.color{:});
                    end
                else
                    irf_plot(hca, r.ts, r.color{:});
                end
                hold(hca,'off');

                if ~isempty(r.legend)
                    irf_legend(hca, r.legend, [0.98 0.15], 'FontSize', r.fontSize*0.8);
                end
                ylabel(hca, r.ylabel, 'interpreter', 'tex', 'FontSize', r.fontSize);
                hca.FontSize = r.fontSize;
            end
        end

        function spectrogram(obj, name, ts, varargin)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            p = inputParser;
            addRequired(p, 'name', @ischar);
            addRequired(p, 'ts');
            addParameter(p, 'fontSize', 16)
            addParameter(p, 'ylabel', 'default ylabel')
            parse(p, name, ts, varargin{:})
            r = p.Results;
            
            hca_func = obj.get_hca_func(r.name);
            obj.execute_list{end + 1} = @execute_spectrogram;

            function execute_spectrogram()
                hca = hca_func();

                ispec = struct();
                ispec.t = r.ts{1}.time.epochUnix;
                ispec.p = r.ts{1}.data;
                ispec.p_label={'flux', r.ts{1}.units};
                ispec.f_label={'Energy', '[eV]'};
                ispec.f = single(r.ts{2}.data);

                [~, hcb] = irf_spectrogram(hca, ispec, 'log');
                % hold(h(4), 'on');
                % irf_plot(h(4), paraTi, 'LineWidth', 1.5);
                % irf_plot(h(4), perpTi, 'LineWidth', 1.5);
                % hold(h(4), 'off');
                hca.YScale = 'log';
                hca.YTick = 10.^[1 2 3 4];
                % irf_legend(hca,'T_{i, ||}',[0.75 0.7],'color','k')
                % irf_legend(hca,'T_{i, \perp}',[0.8 0.7],'color','w')
                ylabel(hca, r.ylabel, 'interpreter', 'tex', 'FontSize', r.fontSize);
                irf_zoom(hca,'y',[10 30000]);
                colormap(hca,'jet');
                caxis(hca,[4.5 7.5]);
                hcb.Label.String(1) = {'log Diff. energy flux'};
                hca.FontSize = r.fontSize;
            end
        end

        function out = mark(obj, tints, varargin)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            p = inputParser;
            addRequired(p, 'tints');
            addParameter(p, 'target', '', @ischar);
            addParameter(p, 'color', 'rgbcmy')
            addParameter(p, 'intervals', false)
            addParameter(p, 'edges', false)
            addParameter(p, 'edgeColor', [])
            parse(p, tints, varargin{:})
            r = p.Results;

            hca_func = obj.get_hca_func(r.target);
            % obj.execute_list{end + 1} = @execute_mark;

            % function execute_mark()
            hca = hca_func();

            if isa(r.color, 'char')
                r.color = r.color';
            end

            if isa(r.edgeColor, 'char')
                r.edgeColor = r.edgeColor';
            end

            if isa(tints, 'cell')
                n_tints_type = length(tints);
                out = cell(1, n_tints_type);
                for t = 1:n_tints_type
                    out{t} = mark_type(tints{t}, r.color(t,:), get_edge_color(t));
                end
            else
                out = mark_type(tints, r.color(1,:), get_edge_color(1));
            end

            function hdl = mark_type(tints, c, edgeColor)
                if r.intervals
                    n_tints = length(tints)/2;
                    hdl = gobjects(1, n_tints);
                    for i = 1:n_tints
                        tint = select_tint(tints, i);
                        hdl(i) = irf_pl_mark(hca, tint, c, 'EdgeColor', edgeColor);
                    end
                else
                    hdl = irf_pl_mark(hca, tints.epochUnix, c);
                end
            end

            function edgeColor = get_edge_color(i)
                if r.edges
                    if isempty(r.edgeColor)
                        edgeColor = r.color(i,:);
                    else
                        edgeColor = r.edgeColor(i,:);
                    end
                else
                    edgeColor = 'none';
                end
            end
        end

        function show(obj, tint, varargin)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            p = inputParser;
            addRequired(p, 'tint');
            addParameter(p, 'anchorPos', [10, 50])
            addParameter(p, 'panelDims', [600, 800])
            parse(p, tint, varargin{:})
            r = p.Results;

            n_panels = length(keys(obj.hcas));
            obj.plot_axes = irf_plot(n_panels, 'newfigure');
            % obj.plot_axes(1)
            title(obj.plot_axes(1), obj.plot_title);
            
            % Set figure pos
            obj.fig = gcf;
            set(obj.fig, 'Position', [r.anchorPos(1), r.anchorPos(2), r.panelDims(1), r.panelDims(2)]);
        
            obj.update();
            obj.original_zoom = tint;
            obj.zoom(r.tint);
        end

        function update(obj)
            n_execute = length(obj.execute_list);
            for i = 1:n_execute
                execute = obj.execute_list{i};
                execute();
            end
            obj.execute_list = {};
        end

        function zoom(obj, tint)
            % Changes to all figure
            irf_plot_axis_align
            irf_zoom(obj.plot_axes, 'x', tint);
            irf_zoom(obj.plot_axes, 'y');
            %irf_pl_number_subplots(h);
            irf_timeaxis(obj.plot_axes);
        end

        function pts = ginput_mark(obj, varargin)

            p = inputParser;
            addParameter(p, 'target', '', @ischar);
            addParameter(p, 'color', 'r')
            addParameter(p, 'intervals', false)
            addParameter(p, 'preMark', [])
            parse(p, varargin{:})
            r = p.Results;

            hca_func = obj.get_hca_func(r.target);
            hca = hca_func();
            % axes(obj.plot_axes);
            % obj.plot_axes(1).userdata
            % ref_pt = get(hca,'userdata').zoom_x(1);
            zoom_coef = 0.1;
            % userdata = get(obj.fig, 'userdata');
            t_ref = get(obj.fig, 'userdata').t_start_epoch;
            % ref_pt = irf_time(x_zoom(1), 'epoch>ttns');
            % hdl_mark = [];
            % hdl_start = [];
            pts = [];
            if r.intervals
                n_pt_batch = 2;
            else
                n_pt_batch = 1;
            end
            if ~isempty(r.preMark)
                count = length(r.preMark)/2 + 1;
                pts = reshape(irf_time(r.preMark, 'ttns>epoch'), n_pt_batch, [])';
                if n_pt_batch == 2
                    hdl_mark = irf_pl_mark(hca, pts, r.color, 'EdgeColor', r.color);
                else
                    hdl_mark = irf_pl_mark(hca, pts, r.color);
                end
            else
                count = 1;
            end
            while true
                pt = zeros(1, n_pt_batch);
                i_pt = 1;
                while i_pt <= n_pt_batch
                    [x, ~, button] = ginput(1);
                    user_zoom(button);
                    if isempty(button) % Return
                        if ~isempty(pts)
                            pts = sort(EpochTT(irf_time(pts, 'epoch>ttns')));
                        end
                        return
                    end
                    tmp_pt = t_ref + x;
                    if button == 1
                        pt(i_pt) = tmp_pt;
                        if n_pt_batch == 2 
                            if i_pt == 1
                                hdl_start = irf_pl_mark(hca, tmp_pt, r.color);
                            else
                                delete(hdl_start);
                            end
                        end
                        i_pt = i_pt + 1;
                    elseif button == 3 || button == 8 % Delete
                        if i_pt == 1 && count > 1 % Delete finshed
                            [~, idx] = min(min(abs(pts - tmp_pt), [], 2), [], 1);
                            pts(idx, :) = [];
                            delete(hdl_mark(idx))
                            hdl_mark(idx) = [];
                            count = count - 1;
                        else % Delete beginning
                            pt = zeros(1, n_pt_batch);
                            i_pt = 1;
                            delete(hdl_start);
                        end
                    end
                end
                pts(count, 1:n_pt_batch) = pt;
                if n_pt_batch == 2
                    fprintf('Selected %.1fs interval.\n', pt(2) - pt(1))
                    hdl_mark(count) = irf_pl_mark(hca, pt, r.color, 'EdgeColor', r.color);
                else
                    hdl_mark(count) = irf_pl_mark(hca, pt, r.color);
                end
                count = count + 1;
            end

            function user_zoom(button)
                if button == 28 % Left => Pan left
                    [zoom_x, dt] = get_zoom_x_and_dt(zoom_coef);
                    obj.zoom([zoom_x(1) + -dt, zoom_x(2) + -dt])
                end
                if button == 29 % Right => Pan right
                    [zoom_x, dt] = get_zoom_x_and_dt(zoom_coef);
                    obj.zoom([zoom_x(1) + dt, zoom_x(2) + dt])
                end
                if button == 30 % Up => Zoom in
                    [zoom_x, dt] = get_zoom_x_and_dt(zoom_coef/(1 + 2*zoom_coef));
                    obj.zoom([zoom_x(1) + dt, zoom_x(2) + -dt])
                end
                if button == 31 % Down => Zoom out
                    [zoom_x, dt] = get_zoom_x_and_dt(zoom_coef);
                    obj.zoom([zoom_x(1) + -dt, zoom_x(2) + dt])
                end
                if button == 114 % R => restore zoom
                    obj.zoom(obj.original_zoom)
                end

                function [zoom_x, dt] = get_zoom_x_and_dt(z_coef)
                    zoom_x = get(hca,'userdata').zoom_x;
                    dt = (zoom_x(2) - zoom_x(1))*z_coef;
                end
            end
        end

        function save(obj, filename, tint)

            if nargin > 2
                obj.zoom(tint)
            end
            saveas(obj.fig, filename)
        end
    end
end

