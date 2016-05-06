
clear all;  close all;  clc;

data_folder = 'J:\MEG_Research\SEF\SEF-TSSS\Patients\oldHPI_03';
addpath('C:\Program Files\MATLAB\R2012b\toolbox_add_on\tSSS_chk')
cd(data_folder)
filenames      = dir('*.mat');

%%
for loop = 1:3 %sort([[1:7:70],[2:7:70]])
    load(filenames(loop).name)
    [a b c]     =   fileparts(filenames(loop).name);

    disp('#######################################')
    disp(['****    ', b ,'        ****'])
    disp('#######################################')
    
    pow{loop}= powCalc2(saveData.dataArtRej);
    pow{loop}.filename = filenames(loop).name;
    
end

% save pow pow -v7.3

%% Power plots

clrOp1 = rgb('dodgerBlue');%[0.7 0.2 0.2];
clrOp2 = rgb('Black');%[0.4 0.4 0.4];
clrOp3 = rgb('Crimson');
clrOp4 = rgb('Olive');
lw=3;

hFig = figure(1)
set(hFig, 'Position', [10 80 1596 1024]);

for loop=1:3
    
   
    [a b c]     =   fileparts(pow{1,loop}.filename)
    subplot(1,3,loop)
    semilogy(pow_new{1,loop}.freq, mean(pow_new{1,loop}.powspctrm), 'color',clrOp1,'linewidth', lw), hold on
    semilogy(pow{1,loop}.freq, mean(pow{1,loop}.powspctrm),'-', 'color',clrOp3,'linewidth', lw);
    % semilogy(pow{3}.freq, mean(pow{3}.powspctrm),'-', 'color',clrOp3,'linewidth', lw);
    % semilogy(pow{4}.freq, mean(pow{4}.powspctrm),'-', 'color',clrOp4,'linewidth', lw);
%     hlegend = legend({'DBS OFF','DBS ON'});
axis([0 100 10e-25 40-20])

    % hlegend = legend({'UPDRS III'},'Interpreter','latex');
%     set(hlegend,'Fontsize',18);
%     set(hlegend,'Fontname','Helvetica');
%     set(hlegend,'Location','Northeast');
%     legend boxoff
%     box off
    % xlim([0 45])
%     set(gca,'XTick',[0:50:200]);
%     set(gca,'Fontsize', 18)
%     set(gca,'Fontweight','bold');
%     set(gca,'Fontname','Helvetica');

%     hxlabel = xlabel('Frequency Hz');
    % set(hxlabel,'Fontsize',24);
    % set(hxlabel,'FontWeight','bold');
    % set(hxlabel,'Fontname','Helvetica');
%     title(b(1:3))
%     hylabel = ylabel(['Spectral density']);%  fT/cm- $$\sqrt{Hz}$$'],'Interpreter', 'latex');
%     export_fig( gcf, ['PowPlot-',b(1:3)],'-transparent', ...
%             '-painters','-pdf', '-r250' ); 


end

%%
set(gcf, 'Color', 'None'); % white bckgr
set(gca, 'Color', 'None'); % white bckgr
export_fig( gcf, 'SEF-DBSON-topo130Hz','-transparent', ...
        '-painters','-pdf', '-r250' ); 
    
%%    

cfg_sel = [];
cfg_sel.trials = [240];
trlData{1,1} = ft_selectdata(cfg_sel, dat_fT{1,1}.dataArtRej);
trlData{1,2} = ft_selectdata(cfg_sel, dat_fT{1,2}.dataArtRej);


%%
cfg                 = [];
cfg.layout          = 'neuromag306cmb.lay';
% cfg.channel         = 'MEGGRAD';
% cfg.hlim            = [-0.1 0.2];
% cfg.vlim            = [-2000 2000];
figure,
ft_multiplotER(cfg,ft_combineplanar([],pow{2}));

%%
hlegend = legend({'tSSS', 'no tSSS'});
set(hlegend,'Fontsize',18);
set(hlegend,'Fontname','Helvetica');
set(hlegend,'Location','Northeast');
legend boxoff
box off

export_fig( gcf, 'multiER-full-LH','-transparent', ...
        '-painters','-pdf', '-r250' ); 

%% Timelock analysis

cfgtlck = [];
cfgtlck.channel = 'MEGGRAD';
cfgtlck.removemean = 'yes';

for loop = 1:2
    tlck{loop}     = ft_timelockanalysis(cfgtlck,dat_fT{loop}.dataArtRej);
    tlckCmb{loop}  = ft_combineplanar([], tlck{loop});
end

%% TFR 














