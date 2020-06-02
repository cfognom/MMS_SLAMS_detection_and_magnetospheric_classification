function MMS_PlotSLAMS_SC1(tint)

global sw B1gse vi1gse ve1gse ni1 ne1 Tiperp1 Tipar1 Teperp1 Tepar1 pressuretensori1gse pressuretensore1gse Eispectomni1 Eespectomni1 Ei1 Ee1;
global R1 R2;
global b_smooth SLAMS_time score_TS limit_ts;

%% Initialize
delta_epochtt = 10*60; % 10 min
win = 13*20; % 20 min

% %% Pick out tint
% tint = [irf_time(starttime) irf_time(stoptime)];

%% Calculate
mp = 1.67e-27;
mu0 = 4*pi*1e-7;
%kB
vi1abs = irf_abs(vi1gse);
Wkin1_data = ni1.data.*vi1abs.data.^2*mp/2*1e21;
Wkin1_smooth_data = movmean(Wkin1_data,win);
Wkin1_ratio_data = Wkin1_data./Wkin1_smooth_data;
Wkin1 = TSeries(ni1.time,Wkin1_data);
Wkin1_smooth = TSeries(ni1.time,Wkin1_smooth_data);
Wkin1_ratio = TSeries(ni1.time,Wkin1_ratio_data);
B1abs = irf_abs(B1gse);
RE1 = R1;
RE1.data = R1.data/6378;




%% Initialize plotting
FS1 = 16;
FS2 = 0.8*FS1;
%  Prepare omnidirectional ion energy flux
ispec=struct('t', Eispectomni1.time.epochUnix);
ispec.p = Eispectomni1.data;
ispec.p_label={'flux', Eispectomni1.units};
ispec.f_label={'Energy', '[eV]'};
ispec.f = single(Ei1.data);


%% Plot
nplot = 7;
%h = irf_plot(nplot,'newfigure');
h = irf_plot(nplot, 'newfigure');
fig = gcf;
set(fig, 'Position', [10, 10, 600, 800]);

% Plot Bgse
hca = irf_panel('B1');
% irf_plot(hca,B1gse);
% hold(hca,'on');
irf_plot(hca,B1abs);
hold(hca,'on');
irf_plot(hca,b_smooth)
irf_plot(hca,limit_ts)
if ~isempty(SLAMS_time)
    irf_pl_mark(hca, SLAMS_time, 'r')
    for i = 1:2:length(SLAMS_time)
        irf_pl_mark(hca, [SLAMS_time(i), SLAMS_time(i + 1)], 'r')
    end
end
hold(hca,'off');
hca.FontSize = FS1;
l1 = ylabel(hca,'B_{GSE} (nT)','FontSize',FS1);
set(l1,'interpreter','tex');
hca.FontSize = FS1;

% Plot ni
hca = irf_panel('ni1');
irf_plot(hca,ni1);
hold on;
irf_plot(hca,ne1);
hold off;
ylabel(hca,'n_i (cm^{-3})','FontSize',FS1);
hca.FontSize = FS1;

% Plot ion temperature
hca = irf_panel('Ti1');
irf_plot(hca,Tipar1);
hold(hca,'on');
irf_plot(hca,Tiperp1);
hold(hca,'off');
irf_legend(hca,{'T_{ipar}','T_{iperp}'},[0.98 0.15],'FontSize',FS2);
l1 = ylabel(hca,'T_{i} (eV)','FontSize',FS1);
set(l1,'interpreter','tex');
hca.FontSize = FS1;

% Plot omnidirectional ion energy flux
hca = irf_panel('Ei_omni');
[hspect hcb] = irf_spectrogram(hca, ispec, 'log');
% hold(h(4), 'on');
% irf_plot(h(4), paraTi, 'LineWidth', 1.5);
% irf_plot(h(4), perpTi, 'LineWidth', 1.5);
% hold(h(4), 'off');
hca.YScale = 'log';
hca.YTick = 10.^[1 2 3 4];
% irf_legend(hca,'T_{i, ||}',[0.75 0.7],'color','k')
% irf_legend(hca,'T_{i, \perp}',[0.8 0.7],'color','w')
ylabel(hca,{'W_i (eV)'},'Interpreter','tex','FontSize',FS1);
irf_zoom(hca,'y',[10 30000]);
colormap(hca,'jet');
caxis(hca,[4.5 7.5]);
hcb.Label.String(1) = {'log Diff. energy flux'};
hca.FontSize = FS1;

% Plot vi
hca = irf_panel('vi1');
irf_plot(hca,vi1gse);
hold(hca,'on');
irf_plot(hca,vi1abs);
hold(hca,'off');
%irf_zoom(hca,'y',yrange);
irf_legend(hca,{'v_{ix}','v_{iy}','v_{iz}'},[0.98 0.1],'FontSize',FS2);
l1 = ylabel(hca,'v_{i} (kms^{-1})','FontSize',FS1);
set(l1,'interpreter','tex');
hca.FontSize = FS1;

% Plot position
hca = irf_panel('R1');
irf_plot(hca,RE1);
%irf_zoom(hca,'y',yrange);
irf_legend(hca,{'x','y','z'},[0.98 0.1],'FontSize',FS2);
l1 = ylabel(hca,'R_{GSE} (R_E)','FontSize',FS1);
set(l1,'interpreter','tex');
hca.FontSize = FS1;

% % Plot classification
hca = irf_panel('Classification');
irf_plot(hca,score_TS);
hca.set('ColorOrder', [1 0 0; 0 0.5 0; 0 0 1])
%irf_zoom(hca,'y',yrange);
irf_legend(hca,{'Msphere', 'SW','Msheath'},[0.98 0.1],'FontSize',FS2);
l1 = ylabel(hca,'Score','FontSize',FS1);
set(l1,'interpreter','tex');
hca.FontSize = FS1;


% Changes to all figure
irf_plot_axis_align
irf_zoom(h,'x',tint);
irf_zoom(h,'y');
%irf_pl_number_subplots(h);
irf_timeaxis(h);


title(h(1), 'MMS1');
%irf_plot_axis_align(h(1:nplot))



return;