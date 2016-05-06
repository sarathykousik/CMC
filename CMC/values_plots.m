

%% CMC


maxCMC_PD_res    =  [0.0817    0.0794    0.1418];
maxCMCmax_PD_scr =  [0.0538    0.0741    0.0750];
    
maxCMC_ctrl_res  =  [0.1357    0.1670    0.1651];
maxCMC_ctrl_scr  =  [0.0643    0.0666    0.0788];

% PD #18, 19, 21
updrs = [13 7 ;
         12 14;
         12 13];
     

%% Only DBS & MED

maxCMC_PD_res    =  [0.0817   0.1418];
maxCMCmax_PD_thresh =  [0.0538   0.0750];
    
maxCMC_ctrl_res  =  [0.1357   0.1651];
maxCMC_ctrl_thresh  =  [0.0643   0.0788];

%%
figure,
bar([maxCMC_PD_res; maxCMC_ctrl_res])
ylim([0 0.2])

%% Get max freq/CMC


set_1 = {'MEG0412+0413', 'MEG0422+0423', 'MEG0432+0433', 'MEG0442+0443',...
                'MEG0632+0633', 'MEG0712+0713'}
set_2 = {'MEG0412+0413', 'MEG0422+0423', 'MEG0432+0433', 'MEG0442+0443', ...
                'MEG0632+0633', 'MEG0712+0713', 'MEG0742+0743'}
set_3 = {'MEG0432+0433', 'MEG0442+0443', 'MEG0712+0713', 'MEG0742+0743', ...
                'MEG1812+1813', 'MEG1822+1823'}
set_4 = {'MEG0432+0433', 'MEG0712+0713', 'MEG0742+0743', 'MEG1822+1823'};

L_par_ch = union(union(union(set_1,set_2),set_3),set_4);

%% calculate total no of trials

for subjLoop=1:10
    for condLoop=1:7
    
        trlLenPD(subjLoop, condLoop) = length(PDdata_trl{subjLoop,condLoop}.trl);
        trlLenCtrl(subjLoop, condLoop) = length(Ctrldata_trl{subjLoop,condLoop}.trl);
        
        
    end
end

PD_coh_thresh_1 = 1-(0.01^(1./(mean(trlLenPD(6:10,1),1)-1)))
PD_coh_thresh_6 = 1-(0.01^(1./(mean(trlLenPD(5:10,6),1)-1)))
PD_coh_thresh_7 = 1-(0.01^(1./(mean(trlLenPD(5:10,7),1)-1)))


%% Ch from GGA: focal
cond_1 = {'MEG0432+0433', 'MEG0442+0443','MEG0232+0233', 'MEG1622+1623'}
cond_2 = {'MEG0432+0433', 'MEG0442+0443', 'MEG1812+1813','MEG0232+0233'}
cond_3 = {'MEG0432+0433', 'MEG0442+0443', 'MEG1812+1813'}
cond_4 = {'MEG0432+0433', 'MEG0442+0443', 'MEG1812+1813'}
cond_5 = {'MEG0432+0433', 'MEG0712+0713'}
cond_6 = {'MEG0432+0433', 'MEG0442+0443', 'MEG1812+1813', 'MEG1822+1823'}
cond_7 = {'MEG0432+0433', 'MEG0442+0443'}

ch_GGA = unique(horzcat(cond_1,cond_2,cond_3,cond_4,cond_5,cond_6,cond_7));

%%

clrOp1 = rgb('Crimson').*255
clrOp2 = rgb('DimGray').*255 %[0.4 0.4 0.4];
clrOp3 = rgb('DarkMagenta').*255

clrOp4 = rgb('DarkCyan').*255
clrOp5 = rgb('Navy').*255%[0.7 0.2 0.2];
clrOp6 = rgb('OrangeRed').*255

%%
figure(1)
cfg                 = [];
cfg.parameter       = 'powspctrm_z';
cfg.layout          = 'neuromag306cmb.lay';
cfg.channel         = ch_GGA;
cfg.graphcolor      =  [clrOp4; clrOp5; clrOp6];%[clrOp1; clrOp2; clrOp3];%
cfg.linewidth       = 4;
cfg.xlim            = [1 45];
cfg.ylim            = [0.03 0.1]
ft_singleplotER(cfg, Ctrl_CMC_grndAvg{[1,6,7]})


line([0 50], [PD_coh_thresh_1 PD_coh_thresh_1(1)], 'color', clrOp1,...
                'linewidth', 2), hold on
line([0 50], [PD_coh_thresh_6  PD_coh_thresh_6], 'color', clrOp2,...
                'linewidth', 2)
line([0 50], [PD_coh_thresh_7 PD_coh_thresh_7], 'color', clrOp3,...
    'linewidth', 2)
            
set(gca,'XTick',[0:10:50],'Fontsize',18, 'FontWeight', 'bold')
set(gca,'YTick',[0:0.05:0.15],'Fontsize',18, 'FontWeight', 'bold')
            
htitle = title('Ctrl');
set(htitle,'Fontsize',20);
set(htitle,'FontWeight','bold');
set(htitle,'Fontname','Helvetica');

hlegend = legend('PD-DBS ON', 'PD-DBS OFF (120min)','PD-MED ON')
set(hlegend,'Fontsize',18);
set(hlegend,'Fontname','Helvetica');
set(hlegend,'FontWeight','bold');
set(hlegend,'Location','Northeast');
legend boxoff

set(gca, 'Color', 'None'); % white bckgr
set(gcf, 'Color', 'None'); % white bckgr
% export_fig( gcf, 'Ctrl_CMC_singleplotGA','-painters','-pdf', '-r250' ); 
    
    

%%            
figure(2)
cfg                 = [];
cfg.parameter       = 'powspctrm';
cfg.layout          = 'neuromag306cmb.lay';
cfg.channel         = ch_GGA;
cfg.graphcolor      = [clrOp3; clrOp4];
cfg.linewidth           = 4;
cfg.xlim            = [1 45];
cfg.ylim            = [0.03 0.15]
ft_singleplotER(cfg, Ctrl_CMC_grndAvg{[1,3]})
legend('Control - condition #1', 'Control - condition #1')
line([0 50], [maxCMC_ctrl_thresh(1) maxCMC_ctrl_thresh(1)], 'color', clrOp3,...
                'linewidth', 2), hold on
line([0 50], [maxCMC_ctrl_thresh(2) maxCMC_ctrl_thresh(2)], 'color', clrOp4,...
                'linewidth', 2)
            
set(gca,'XTick',[0:10:50],'Fontsize',18, 'FontWeight', 'bold')
set(gca,'YTick',[0:0.05:0.15],'Fontsize',18, 'FontWeight', 'bold')

htitle = title('Control');
set(htitle,'Fontsize',20);
set(htitle,'FontWeight','bold');
set(htitle,'Fontname','Helvetica');

hlegend = legend('Control - DBS ON', 'Control - DBS OFF')
set(hlegend,'Fontsize',18);
set(hlegend,'Fontname','Helvetica');
set(hlegend,'FontWeight','bold');
set(hlegend,'Location','Northeast');
legend boxoff
set(gca, 'Color', 'None'); % white bckgr
set(gcf, 'Color', 'None'); % white bckgr
% export_fig( gcf, 'Ctrl_CMC_plot',...
%         '-painters','-pdf', '-r250' ); 

%% multiplotER

% figure
% 
% cfg                     = [];
% cfg.parameter           = 'powspctrm';
% cfg.layout              = 'neuromag306cmb.lay';
% cfg.graphcolor          = 'brkcmg'
% cfg.linewidth           = 1.5;
% cfg.graphcolor          = [clrOp1; clrOp2; clrOp3; clrOp4];
% % cfg.channel           = L_par_ch;

% cfg.ylim            = [0.03 0.15]
% % ft_multiplotER(cfg, ft_combineplanar([],freq_fourier))
% ft_multiplotER(cfg, PD_CMC_grndAvg{[1 3]})

figure
cfg                 = [];
cfg.parameter       = 'powspctrm_z';
cfg.layout          = 'neuromag306cmb.lay';
cfg.hlim                = [1 50];
cfg.linewidth           = 1.5;
cfg.vlim            = [0.03 0.15]
cfg.graphcolor          = [clrOp4; clrOp5; clrOp6];
ft_multiplotER(cfg, Ctrl_CMC_grndAvg{[1,6,7]});%, Ctrl_CMC_grndAvg{[1,3]})
title('Ctrl')
set(gca, 'Color', 'None'); % white bckgr
set(gcf, 'Color', 'None'); % white bckgr
export_fig( gcf, 'multiplotER-Ctrl_CMC_plot','-painters','-pdf', '-r250' ); 

%%
summ_plotVal=[mean(PD_CMC_summary([6:10],1)),mean(PD_CMC_summary([5:10],[6,7])); ...
        mean(Ctrl_CMC_summary(:,[1,6,7]))];
figure(1),
h1=bar(summ_plotVal);
% ylim([0 0.08])
set(gca,'XTick',[1 2])
set(gca,'XTickLabel',{'PD', 'Control'})
set(gca,'YTick',[0:0.05:0.2])
legend('DBS ON', 'DBS OFF(120min)','MED ON'), legend boxoff
hold on

hold on
plot([1:3], mean(updrs(:,[1 6 7])))


%%
x_vec =[1:3,6:8];
[ax,b,p] = plotyy(x_vec,reshape(summ_plotVal',1,[]),x_vec,...
    [mean(updrs(:,[1 6 7])) NaN NaN NaN],'bar','plot');

set(ax(2),'XTick',[2 7]); 
set(ax(2),'XtickLabel', [])
set(ax(2),'XTickLabel',{'PD', 'Control'})

% xlim([0 3])
ylim(ax(1), [0 0.1])
set(ax(1),'YTick',[0:0.02:0.1]); 
set(ax(1),'YTickLabel',[0:0.02:0.1])

ylim(ax(2), [-20 50])
set(ax(2),'YTick',[0:10:40]); 
set(ax(2),'YTickLabel',[0:10:40])

set(gca, 'Color', 'None'); % white bckgr
set(gcf, 'Color', 'None'); % white bckgr
export_fig( gcf, 'barplot_all',...
        '-painters','-pdf', '-r250' ); 

%% Topo Coh - Alpha/beta/gamma

cfg                 = [];
cfg.parameter       = 'powspctrm_z';
cfg.layout          = 'neuromag306cmb.lay';
% cfg.graphcolor      = 'brkcmgy';
cfg.channel         = {'all', '-MEG1722+1723'};%ch_GGA;
% cfg.zlim = [2e-24 10e-24] %beta
% cfg.zlim = [2.5e-24 5e-23] %gamma
% cfg.zlim = [1e-24 4e-24] %High gamma
% cfg.hlim = [13:30];
cfg.contournum          = 4;
 cfg.gridscale = 200;
cfg.zlim            = [0.02 0.06]
cfg.colorbar = 'no';
% loop=1;
% hFig=figure;
% set(hFig, 'Position', [100 50 1556 1000]);


for condLoop= 1
       hFig=figure(condLoop);
        set(hFig, 'Position', [100 50 1000 1000]);  
        ft_topoplotER(cfg, Ctrl_CMC_beta_grndAvg{condLoop})
        title(['Ctrl-',num2str(condLoop)])
        set(gcf, 'Color', 'None'); % white bckgr
        set(gca, 'Color', 'None'); % white bckgr

        export_fig( gcf, ['Ctrlcghk-',num2str(condLoop)],'-painters','-png', '-r100' ); 

  end
% suptitle('GGA-beta-CMC-Norm')


%% CMC Summary from individual

cfg_sel      = [];
cfg_sel.frequency = 1;

chan_pos_L = match_str(PDcoh_beta_cmb{1}.label, ch_GGA);
% chan_pos_R = match_str(PDcoh_cmb{1}.label, R_par_ch);

freq_range = 'beta'
foi=[13,30]; supText = [freq_range,'FT-CMC-mtm '];
% foi_pts = find(Ctrlcoh_beta_cmb{1,1}.freq>=foi(1) & Ctrlcoh_beta_cmb{1,1}.freq<=foi(2));
time_vector = [1 2 3 4 5 6 8];

for subjLoop = 1:10
    for condLoop = 1:7
    
%         [temp_maxCMCPD, PDmaxFreq]    = max(PDcoh_cmb{subjLoop, condLoop}.powspctrm(chan_pos_L,foi_pts), [], 2);
%         PD_CMC_summary(subjLoop,condLoop)      = mean(temp_maxCMCPD);
%         
%         [temp_maxCMCPD_Ipsi, PDmaxFreq_Ipsi]    = max(PDcoh_cmb{subjLoop, condLoop}.powspctrm(chan_pos_R,foi_pts), [], 2);
%         PD_CMC_summary_Ipsi(subjLoop,condLoop)      = mean(temp_maxCMCPD_Ipsi);
% 
%         [temp_maxCMCCtrl, CtrlmaxFreq]    = max(Ctrlcoh_cmb{subjLoop, condLoop}.powspctrm(chan_pos_L,foi_pts), [], 2);
%         Ctrl_CMC_summary(subjLoop,condLoop)      = mean(temp_maxCMCCtrl);
%         
%         [temp_maxCMCCtrl_Ipsi, CtrlmaxFreq_Ipsi]    = max(Ctrlcoh_cmb{subjLoop, condLoop}.powspctrm(chan_pos_R,foi_pts), [], 2);
%         Ctrl_CMC_summary_Ipsi(subjLoop,condLoop)      = mean(temp_maxCMCCtrl_Ipsi);
%         
        PD_CMC_summary(subjLoop,condLoop)    = ...
            eval(['mean(PDcoh_', freq_range,'_cmb{subjLoop, condLoop}.powspctrm_z(chan_pos_L));']);
%         PD_CMC_summary_Ipsi(subjLoop,condLoop) =...
%             eval(['mean(PDcoh_', freq_range,'_cmb{subjLoop, condLoop}.powspctrm(chan_pos_R));']);
        Ctrl_CMC_summary(subjLoop,condLoop)      = ...
            eval(['mean(Ctrlcoh_', freq_range,'_cmb{subjLoop, condLoop}.powspctrm_z(chan_pos_L));']);
%         Ctrl_CMC_summary_Ipsi(subjLoop,condLoop)      = ...
%             eval(['mean(Ctrlcoh_', freq_range,'_cmb{subjLoop, condLoop}.powspctrm(chan_pos_R));']);

    end
end

%%

clrOp1 = rgb('dodgerBlue');%[0.7 0.2 0.2];
clrOp2 = rgb('DimGray');%[0.4 0.4 0.4];
time_vector = [1:6,8]
axis_set = [0.5 8.5 0.02 0.07]
% axis_set = [0.5 7.5 0 1.5]
conds = [1:7];
PD_subjs = [5:10];
hold on
% subplot(2,1,1)
% p1=plot((PD_CMC_summary(PD_subjs,conds))', '-o'); hold on


% subplot(2,1,2)
% p2=plot((Ctrl_CMC_summary(:,conds))', '-o'); hold on
e2=errorbar(time_vector,mean(Ctrl_CMC_summary(:,conds)),sem(Ctrl_CMC_summary(:,conds)), '-s');
set(e2,'LineWidth', 4.5 )
set(e2, 'MarkerSize', 18,'Color', clrOp2 , 'MarkerFaceColor', clrOp2 , ...
    'MarkerEdgeColor', clrOp2);
axis(axis_set), %title('Control-Contra-mean');


e1 = errorbar(time_vector,mean(PD_CMC_summary(PD_subjs,conds)),sem(PD_CMC_summary(PD_subjs,conds)),'-o');
axis(axis_set), %title('PD-Contra-mean');
set(e1,'LineWidth', 4.5 )
set(e1, 'MarkerSize', 18,'Color', clrOp1 , 'MarkerFaceColor', clrOp1 , ...
    'MarkerEdgeColor', clrOp1);


set(gca,'XTick',[1:6 8])
% set(gca,'XTickLabel',[])
set(gca,'Fontsize', 18, 'fontname', 'Helvetica','FontWeight','bold');


xTick = [1:6 8];
set(gca,'xtick',xTick);
yTick = get(gca,'ytick');
set(gca,'xticklabel',[]);
xTickLabel = {{'DBS ON'}, {'DBS OFF';'(0min)'}, {'DBS OFF';'(30min)'},{'DBS OFF','(60min)'},...
    {'DBS OFF','(90min)'},{'DBS OFF','(120min)'}, 'Med ON'};
for k = 1:length(xTick)
    text(xTick(k),yTick(1)-0.04*(yTick(end)-yTick(1)-0.5),xTickLabel{k},...
        'HorizontalAlignment','center', 'Fontsize', 14, 'fontname', 'Helvetica',...
        'FontWeight','bold');
end
% set(gca,'XTickLabel',...
%     {'DBS ON', 'DBS OFF(0min)', 'DBS OFF(30min)', 'DBS OFF(60min)', 'DBS OFF(90min)'...
%     'DBS OFF(120min)', '','MED ON'},'fontsize', 10, 'fontname', 'Helvetica') ;

set(gca,'YTick',[0:0.02:0.08])
set(gca,'YTickLabel',[0:0.02:0.08])
set(gca,'Fontsize', 18)
set(gca,'Fontweight','bold');
set(gca,'Fontname','Helvetica');

hxlabel = xlabel('Conditions')
set(hxlabel,'Fontsize',24);
set(hxlabel,'FontWeight','bold');
set(hxlabel,'Fontname','Helvetica');


hylabel = ylabel('Beta CMC')
set(hylabel,'Fontsize',24);
set(hylabel,'FontWeight','bold');
set(hylabel,'Fontname','Helvetica');

hlegend = legend({'Control','PD'});
% hlegend = legend({'UPDRS III'},'Interpreter','latex');
set(hlegend,'Fontsize',18);
set(hlegend,'Fontname','Helvetica');
set(hlegend,'Location','Northeast');
legend boxoff


hxlabel = xlabel('Conditions')
set(hxlabel,'Fontsize',24);
set(hxlabel,'FontWeight','bold');
set(hxlabel,'Fontname','Helvetica');

htitle = title('Beta - CMC');
set(htitle,'Fontsize',24);
set(htitle,'FontWeight','bold');
set(htitle,'Fontname','Helvetica');
% 
export_fig -painters -r600 -q101 beta_CMC_trend.pdf
set(gcf, 'Color', 'None'); % white bckgr
set(gca, 'Color', 'None'); % white bckgr

%%  UPDRS


%% UPDRS- Mean=/- SEM

figure
% contS = [1:10];
subjP = [1:6];
time_vector = [1 2 4 6 8];
clrOpUpdrs = rgb('dodgerBlue');
e1=errorbar(time_vector, mean(updrs(subjP,[1 2 4 6 7]),1), ...
    sem(updrs(subjP,[1 2 4 6 7])), '-ro' );
set(e1, 'MarkerSize', 18,'Color', clrOpUpdrs , 'MarkerFaceColor', clrOpUpdrs , ...
    'MarkerEdgeColor', clrOpUpdrs);
set(e1,'LineWidth', 4.5 )
% hold on
% e2=errorbar(time_vector, mean(gammaControlmean(contS,:),1),  sem(gammaControlmean(contS,:)), '-b.')
% set(e2,'LineWidth', 2.5 ) 
xlim([0.5 8.5])
ylim([0 40])

set(gca,'XTick',[1 2 3 4 5 6 8])
% set(gca,'XTickLabel',[])
set(gca,'YTick',[0:10:40])
% set(gca,'YTickLabel',[])
set(gca,'XTickLabel',...
    {'Stim ON', 'StimOFF(0min)', 'StimOFF(30min)', 'StimOFF(60min)', 'StimOFF(90min)'...
    'StimOFF(120min)', '','MedON'},'fontsize', 10, 'fontname', 'Georgia') ;
xTick = time_vector;
set(gca,'xtick',xTick);
yTick = get(gca,'ytick');
set(gca,'xticklabel',[]);
xTickLabel = {{'DBS ON'}, {'DBS OFF';'(0min)'}, {'DBS OFF';'(60min)'},{'DBS OFF';'(120min)'},...
    'Med ON'};
for k = 1:length(xTick)
    text(xTick(k),yTick(1)-0.04*(yTick(end)-yTick(1)+0.5),xTickLabel{k},...
        'HorizontalAlignment','center', 'Fontsize', 14, 'fontname', 'Helvetica',...
        'FontWeight','bold');
end


hxlabel = xlabel('Conditions')
set(hxlabel,'Fontsize',24);
set(hxlabel,'FontWeight','bold');
set(hxlabel,'Fontname','Helvetica');

hylabel = ylabel('UPDRS III  (Mean \pm SEM)')
set(hylabel,'Fontsize',24);
set(hylabel,'FontWeight','bold');
set(hylabel,'Fontname','Helvetica');
% 
set(gca,'YTickLabel',[0:10:40])
set(gca,'Fontsize', 18)
set(gca,'Fontweight','bold');
set(gca,'Fontname','Helvetica');
% 
htitle = title('Clinical motor state');
set(htitle,'Fontsize',24);
set(htitle,'FontWeight','bold');
set(htitle,'Fontname','Helvetica');


% axis([0.5 8.5 0.5 1]);
% hlegend = legend('PD (n=10)');
% % hlegend = legend({'UPDRS III'},'Interpreter','latex');
% set(hlegend,'Fontsize',18);
% set(hlegend,'Fontname','Helvetica');
% set(hlegend,'Location','Northeast');
% legend boxoff

set(0,'defaultAxesFontName', 'Calibri')
set(gca, 'Color', 'None'); % white bckgr
set(gcf, 'Color', 'None'); % white bckgr
% export_fig -painters -r600 -q101 UPDRS-SEM.pdf
box off
export_fig( gcf, 'UPDRS_new','-transparent', ...
        '-painters','-pdf', '-r250' ); 
% saveas(gcf, 'UPDRS-trans.pdf', 'pdf') 

%% plotyy

fig = figure(1);
time_vector = [1 2 3 4 5 6 8];
PDsubjs = [5:10];
set(fig, 'Position', [10 80 1024 800]);
clf(fig)
hold on;
[AX,H1,H2] = plotyy([time_vector; time_vector]' , ...
    [mean(PD_CMC_summary(PDsubjs,:)), mean(PD_CMC_summary(PDsubjs,:))],...
    time_vector, updrs(:,[1 2 3 4 5 6 7]));

set(H1(1), 'linestyle', '-', 'Marker', 'x', 'color', 'r', 'linewidth', 2);
set(H1(2), 'linestyle', '-', 'Marker', 'o', 'color', 'b', 'linewidth', 2);
set(H2, 'linestyle', '--','Marker', 'v',  'color', 'k', 'linewidth', 2);

errorbar(time_vector, mean(PDITC_meanNeigh([1:10],:),1),  sem(PDITC_meanNeigh([1:10],:)), 'r.');
errorbar(time_vector+0.1, ...
    mean(conITC_meanNeigh([1:10],:),1),  sem(conITC_meanNeigh([1:10],:)), 'b.');

set(fig, 'CurrentAxes', AX(2));
hold on
errorbar([1 2 4 8], mean(updrs(:,[1 2 4 6 7]),1), ...
    sem(updrs(:,[1 2 4 6 7])), '.k');












