
clrOp1 = rgb('Crimson');
clrOp2 = rgb('MidnightBlue');
clrOp3 = rgb('DarkKhaki');

clrOp4 = rgb('dodgerBlue');
clrOp5 = rgb('DimGray');
clrOp6 = rgb('DarkGreen');
clrOp7 = rgb('Maroon');

%% grand avg CMC

cfg = []; 
cfg.parameter = 'powspctrm'

for loop = 1:7
   
    PD_CMC_grndAvg{loop} = ft_freqgrandaverage(cfg, PDcoh_cmb{5:10,loop});
    Ctrl_CMC_grndAvg{loop} = ft_freqgrandaverage(cfg, Ctrlcoh_cmb{:,loop});
    
end

%% CMC 
% Left
L_par_ch = {'MEG0412+0413', 'MEG0422+0423', 'MEG0432+0433', 'MEG0442+0443',...
    'MEG0632+0633', 'MEG0712+0713', 'MEG0742+0743', 'MEG1812+1813', 'MEG1822+1823'};

% Right
R_par_ch = {'MEG0722+0723', 'MEG0732+0733', 'MEG1042+1043', 'MEG1112+1113',...
    'MEG1122+1123', 'MEG1132+1133', 'MEG1142+1143', 'MEG2212+2213', 'MEG2222+2223'};

%%
figure,
bar([maxCMC_PD_res; maxCMC_ctrl_res])
ylim([0 0.2])


%% Singleplot CMC PD

% clrOp1 = rgb('Crimson');
% clrOp2 = rgb('MidnightBlue');
% clrOp3 = rgb('dodgerBlue');%[0.7 0.2 0.2];
% clrOp4 = rgb('DimGray');%[0.4 0.4 0.4];

figure(1)
cfg                 = [];
cfg.parameter       = 'powspctrm';
cfg.layout          = 'neuromag306cmb.lay';
cfg.channel         = L_par_ch;
cfg.graphcolor      = [clrOp1; clrOp2; clrOp4];
cfg.linewidth           = 4;
cfg.xlim            = [1 45];
cfg.ylim            = [0.03 0.15]
ft_singleplotER(cfg, PD_CMC_grndAvg{[1 6 7]})


line([0 50], [0.03 0.03], 'color', clrOp1,...
                'linewidth', 2), hold on
% line([0 50], [0.03 0.03], 'color', clrOp2,...
%                 'linewidth', 2)
            
set(gca,'XTick',[0:10:50],'Fontsize',18, 'FontWeight', 'bold')
set(gca,'YTick',[0:0.05:0.15],'Fontsize',18, 'FontWeight', 'bold')
            
htitle = title('PD');
set(htitle,'Fontsize',20);
set(htitle,'FontWeight','bold');
set(htitle,'Fontname','Helvetica');
    
hlegend = legend('PD - DBS ON', 'PD - DBS OFF', 'PD - MED ON')
set(hlegend,'Fontsize',18);
set(hlegend,'Fontname','Helvetica');
set(hlegend,'FontWeight','bold');
set(hlegend,'Location','Northeast');
legend boxoff

% set(gca, 'Color', 'None'); % white bckgr
% set(gcf, 'Color', 'None'); % white bckgr
% export_fig( gcf, 'PD_CMC_plot',...
%         '-painters','-pdf', '-r250' ); 
%     
    

%% Single plot CMC Contorl
figure(2)
cfg                 = [];
cfg.parameter       = 'powspctrm';
cfg.layout          = 'neuromag306cmb.lay';
cfg.channel         = L_par_ch;
cfg.graphcolor      = [clrOp4; clrOp5; clrOp6];
cfg.linewidth           = 4;
cfg.xlim            = [1 45];
cfg.ylim            = [0.01 0.15]
ft_singleplotER(cfg, Ctrl_CMC_grndAvg{[1 6 7]})
legend('Control - condition #1', 'Control - condition #1')
line([0 50], [0.03 0.03], 'color', clrOp3,...
                'linewidth', 2), hold on
line([0 50], [0.03 0.03], 'color', clrOp4,...
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
export_fig( gcf, 'Ctrl_CMC_plot',...
        '-painters','-pdf', '-r250' ); 

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
% cfg.xlim                = [1 50];
% cfg.ylim            = [0.03 0.15]
% % ft_multiplotER(cfg, ft_combineplanar([],freq_fourier))
% ft_multiplotER(cfg, PD_CMC_grndAvg{[1 3]})

figure
cfg                 = [];
cfg.parameter       = 'powspctrm';
cfg.layout          = 'neuromag306cmb.lay';
cfg.linewidth           = 1.5;
cfg.ylim            = [0.03 0.15]
cfg.graphcolor          = [clrOp1; clrOp2; clrOp3; clrOp4; clrOp5; clrOp6];
ft_multiplotER(cfg, PD_CMC_grndAvg{[1 6 7]}, Ctrl_CMC_grndAvg{[1 6 7]})

set(gca, 'Color', 'None'); % white bckgr
set(gcf, 'Color', 'None'); % white bckgr
% export_fig( gcf, 'multiplotERCtrl_CMC_plot',...
%         '-painters','-pdf', '-r250' ); 

%%

figure(1),
bar([maxCMC_PD_res; maxCMC_ctrl_res])
ylim([0 0.2])
set(gca,'XTick',[1 2])
set(gca,'XTickLabel',{'PD', 'Control'})
set(gca,'YTick',[0:0.05:0.2])
legend('DBS ON', 'MED ON'), legend boxoff
hold on

%%
[ax,b,p] = plotyy([1:2],[maxCMC_ctrl_res],[1:2], mean(updrs),'bar','plot');

set(gca,'XTick',[1.5]); 
set(gca,'XtickLabel', [])
set(gca,'XTickLabel',{'PD'})
xlim([0 3])
ylim(ax(1), [0 0.2])
set(ax(1),'YTick',[0:0.05:0.2]); 
set(ax(1),'YTickLabel',[0:0.05:0.2])

ylim(ax(2), [0 20])
set(ax(2),'YTick',[0:10:20]); 
set(ax(2),'YTickLabel',[0:10:20])

set(gca, 'Color', 'None'); % white bckgr
set(gcf, 'Color', 'None'); % white bckgr
export_fig( gcf, 'barplot_ctrl',...
        '-painters','-pdf', '-r250' ); 

%%

maxCMC_PD_hann_res =

    0.2608    0.2961    0.3007


maxCMC_ctrl_hann_res =

    0.2980    0.3148    0.3162