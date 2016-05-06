clc;    clear all;    close all;

%% Setup folders and files

proc.data_folder    = 'J:\MEG_Research\CMC\raw_ArtRej';
proc.results_folder = 'J:\MEG_Research\CMC\raw_ArtRej\results';

cd(proc.data_folder)
filenames      = dir('*raw.fif');

% for loop = 1: length(filenames)
%     file_sub(loop) = {filenames(loop).name};
% end
% 
% file_sub = reshape(file_sub, 3, [])';
% [row col] = size(file_sub);


%%
% cfg_math=[];
% cfg_math.operation='log';
% cfg_math.parameter='powspctrm';
% freq_log= ft_math(cfg_math,freq_desc);
par_ch = {'MEG0212+0213', 'MEG0232+0233', 'MEG0242+0243', 'MEG0412+0413',...
    'MEG0422+0423', 'MEG0432+0433', 'MEG0442+0443', 'MEG1522+1523', ...
    'MEG1612+1613', 'MEG1622+1623', 'MEG1632+1633', 'MEG1812+1813',...
    'MEG1822+1823', 'MEG1842+1843'};

cfg                  = [];
cfg.parameter        = 'powspctrm';
% cfg.ylim             = [0 20e-24];
cfg.maskstyle               = 'saturation';	
cfg.shading                 = 'interp';
cfg.layout      = 'neuromag306cmb.lay';
figure; ft_multiplotER(cfg,  ft_combineplanar([],freq_desc{1}))



%%

cfg                 = [];
cfg.parameter       = 'powspctrm';
cfg.layout          = 'neuromag306cmb.lay';
% cfg.channel         = 'MEGGRAD'
% cfg.xlim            = [13 30];
% cfg.zlim            = [0 0.15];
% cfg.refchannel = {'EMG'};
figure(10)
for loop = 1:9
    subplot(3,3,loop)
    ft_topoplotER(cfg, coh_sel{loop})
%     title(num2str(loop))
end
suptitle('PD')

figure(11)
for loop = 10:18
    subplot(3,3,loop-9)
    ft_topoplotER(cfg, coh_sel{loop})
%     title(num2str(loop))
end
suptitle('Control')

% export_fig( gcf, 'asdfControl','-transparent', ...
%         '-painters','-pdf', '-r250' ); 

%%

cfg                 = [];
cfg.parameter       = 'cohspctrm';
cfg.channel         = {'MEG0211', 'MEG0221', 'MEG0231', 'MEG0241', 'MEG0411', ...
                       'MEG0421', 'MEG0431', 'MEG0441', 'MEG1611', 'MEG1621', ...
                       'MEG1811', 'MEG1821'};
cfg.layout          = 'neuromag306all.lay';
cfg.refchannel      = {'EMG'};
sel_ch              = ft_selectdata(cfg, coh{1});
sel_cmb             = ft_combineplanar([],sel_ch);

%% EMG spectra

for loop = 1:18
    
    subplot(6,3, loop)
    semilogy(freq_desc{1,loop}.freq, mean(freq_desc{1,loop}.powspctrm(chan_pos_L,:)))
    title(freq_desc{1,loop}.filename)
    
end
 
% export_fig( gcf, 'MEG spectra','-transparent', ...
%         '-painters','-pdf', '-r250' ); 


%% grand avg CMC
% for loop = 1:18
%     coh_cmb_hann{loop}.idx = loop;
% end

% coh_cmb_res = reshape(coh_cmb_hann,3,6)'
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


%% MultiplotER - CMC
clrOp(1,:) = rgb('DimGray');%[0.4 0.4 0.4];
clrOp(2,:) = rgb('Red');%[0.4 0.4 0.4];
clrOp(3,:) = rgb('Cyan');%[0.4 0.4 0.4];
clrOp(4,:) = rgb('LimeGreen');%[0.4 0.4 0.4];
clrOp(5,:) = rgb('Blue');%[0.4 0.4 0.4];
% clrOp(6,:) = rgb('Teal');%[0.4 0.4 0.4];
clrOp(7,:) = rgb('Indigo');%[0.4 0.4 0.4];

cfg                 = [];
cfg.parameter       = 'powspctrm';
cfg.layout          = 'neuromag306cmb.lay';
cfg.graphcolor      = [clrOp1; clrOp2; clrOp3]
cfg.channel         = L_par_ch;
% cfg.highlight = 'on'
% cfg.highlightchannel = [L_par_ch;R_par_ch];
% cfg.highlightsymbol    = 'x';
% cfg.highlightsize     = 12
cfg.xlim            = [13 30];
cfg.ylim            = [0.08 0.13]
cfg.zlim            = [0.05 0.12]
cfg.colorbar = 'yes';
figure
% for loop=1:7
%     subplot(3,3,loop)
    ft_multiplotER(cfg, PD_CMC_grndAvg{[1 6 7]})
% end
title('2Hz-maxGrad')

%% Get max freq/CMC
% L_par_ch = {'MEG0232+0233', 'MEG0242+0243', 'MEG0422+0423', 'MEG0432+0433',...
%     'MEG0442+0443', 'MEG1612+1613', 'MEG1622+1623', 'MEG1632+1633', 'MEG1812+1813',...
%     'MEG1822+1823', 'MEG1842+1843'};
% 
% 
% R_par_ch = {'MEG1112+1113', 'MEG1132+1133', 'MEG1142+1143', 'MEG1332+1333', ...
%     'MEG1342+1343', 'MEG2212+2213', 'MEG2222+2223', 'MEG2232+2233', 'MEG2412+2413',...
%     'MEG2422+2423', 'MEG2442+2443'};
% 
% L_par_2= {'MEG0222+0223', 'MEG0232+0233', 'MEG0332+0333', 'MEG0412+0413', ...
%     'MEG0422+0423', 'MEG0432+0433', 'MEG0442+0443', 'MEG0622+0623', 'MEG0632+0633',...
%     'MEG0642+0643', 'MEG0712+0713', 'MEG0722+0723', 'MEG0732+0733', 'MEG0742+0743',...
%     'MEG1042+1043', 'MEG1622+1623', 'MEG1812+1813', 'MEG1822+1823'};
% new = {'MEG0412+0413', 'MEG0422+0423', 'MEG0432+0433', 'MEG0442+0443', 'MEG0712+0713'};

set_1 = {'MEG0412+0413', 'MEG0422+0423', 'MEG0432+0433', 'MEG0442+0443',...
                'MEG0632+0633', 'MEG0712+0713'};
set_2 = {'MEG0412+0413', 'MEG0422+0423', 'MEG0432+0433', 'MEG0442+0443', ...
                'MEG0632+0633', 'MEG0712+0713', 'MEG0742+0743'};
set_3 = {'MEG0432+0433', 'MEG0442+0443', 'MEG0712+0713', 'MEG0742+0743', ...
                'MEG1812+1813', 'MEG1822+1823'};
set_4 = {'MEG0432+0433', 'MEG0712+0713', 'MEG0742+0743', 'MEG1822+1823'};

L_par_ch_small = union(union(union(set_1,set_2),set_3),set_4);


%% CMC Summary from individual

cfg_sel      = [];
cfg_sel.frequency = 1;

chan_pos_L = match_str(PDcoh_cmb{1}.label, L_par_ch);
chan_pos_R = match_str(PDcoh_cmb{1}.label, R_par_ch);
foi=[13,30]
supText = 'Beta-CMC'
foi_pts = find(Ctrlcoh_cmb{1,1}.freq>=foi(1) & Ctrlcoh_cmb{1,1}.freq<=foi(2));

for subjLoop = 1:10
    for condLoop = 1:7
    
        [temp_maxCMCPD, PDmaxFreq]    = max(PDcoh_cmb{subjLoop, condLoop}.powspctrm(chan_pos_L,foi_pts), [], 2);
        PD_CMC_summary(subjLoop,condLoop)      = mean(temp_maxCMCPD);
        
        [temp_maxCMCPD_Ipsi, PDmaxFreq_Ipsi]    = max(PDcoh_cmb{subjLoop, condLoop}.powspctrm(chan_pos_R,foi_pts), [], 2);
        PD_CMC_summary_Ipsi(subjLoop,condLoop)      = mean(temp_maxCMCPD_Ipsi);

        [temp_maxCMCCtrl, CtrlmaxFreq]    = max(Ctrlcoh_cmb{subjLoop, condLoop}.powspctrm(chan_pos_L,foi_pts), [], 2);
        Ctrl_CMC_summary(subjLoop,condLoop)      = mean(temp_maxCMCCtrl);
        
        [temp_maxCMCCtrl_Ipsi, CtrlmaxFreq_Ipsi]    = max(Ctrlcoh_cmb{subjLoop, condLoop}.powspctrm(chan_pos_R,foi_pts), [], 2);
        Ctrl_CMC_summary_Ipsi(subjLoop,condLoop)      = mean(temp_maxCMCCtrl_Ipsi);

    end
end

%%
axis_set = [0.5 7.5 0.07 0.2]
conds = [1:7];
figure,
subplot(2,2,1)
p1=plot((PD_CMC_summary([5:10],conds))', '-o'); hold on
errorbar(mean(PD_CMC_summary([5:10],conds))',sem(PD_CMC_summary([5:10],conds))', '-or', 'linewidth', 4)
axis(axis_set), title('PD-Contra')
xlabel('Conditons'); ylabel('CMC')

subplot(2,2,2)
p2=plot((Ctrl_CMC_summary(:,conds))', '-o'); hold on
errorbar(mean(Ctrl_CMC_summary(:,conds))',sem(Ctrl_CMC_summary(:,conds)) , '-ob', 'linewidth', 4)
axis(axis_set), title('Control-Contra')
xlabel('Conditons'); ylabel('CMC')

subplot(2,2,3)
p2=plot(PD_CMC_summary_Ipsi([5:10],conds)', '-o'); hold on
errorbar(mean(PD_CMC_summary_Ipsi([5:10],conds))',sem(PD_CMC_summary_Ipsi([5:10],conds)),'-or', 'linewidth', 4)
axis(axis_set), title('PD-ipsi')
xlabel('Conditons'); ylabel('CMC')

subplot(2,2,4)
p2=plot((Ctrl_CMC_summary_Ipsi(:,conds))', '-o'); hold on
errorbar(mean(Ctrl_CMC_summary_Ipsi(:,conds))',sem(Ctrl_CMC_summary_Ipsi(:,conds)),  '-ob', 'linewidth', 4)
axis(axis_set), title('Control-ipsi')
xlabel('Conditons'); ylabel('CMC')
suptitle(supText)

%%
export_fig( gcf, 'Beta-CMC','-transparent', ...
        '-painters','-pdf', '-r250' ); 


%%
figure,
p1=plot(mean(PD_CMC_summary([5:10],:))', '-or'), hold on, title('grip')
% figure
p2=plot(mean(PD_CMC_summary([1:4],:))', '-om'), hold on, title('tap')
% figure
p3=plot(mean(Ctrl_CMC_summary)', '-db'), hold on
legend([p1(1) p2(1) p3(1)],{'PD-grip','PD-tap','Control',})
axis([0.5 7.5 0.05 0.2]), title('beta')
xlabel('Conditons'); ylabel('CMC')

% export_fig( gcf, 'CMC-alpha-power-move','-transparent', ...
%         '-painters','-pdf', '-r250' ); 


%% Power summary values from individual

cfg_sel      = [];
cfg_sel.frequency = 1;

chan_pos_L = match_str(PDfreq_cmb{1}.label, L_par_ch);
chan_pos_R = match_str(PDfreq_cmb{1}.label, R_par_ch);
foi=[31,45]
foi_pts = find(Ctrlfreq_cmb{1,1}.freq>=foi(1) & Ctrlfreq_cmb{1,1}.freq<=foi(2));

for subjLoop = 1:10
    for condLoop = 1:7
    
        [temp_maxfreqPD, PDmaxFreq]    = max(PDfreq_cmb{subjLoop, condLoop}.powspctrm(chan_pos_L,foi_pts), [], 2);
        PD_freq_summary(subjLoop,condLoop)      = mean(temp_maxfreqPD);
        
        [temp_maxfreqPD_Ipsi, PDmaxFreq_Ipsi]    = max(PDfreq_cmb{subjLoop, condLoop}.powspctrm(chan_pos_R,foi_pts), [], 2);
        PD_freq_summary_Ipsi(subjLoop,condLoop)      = mean(temp_maxfreqPD_Ipsi);

        [temp_maxfreqCtrl, CtrlmaxFreq]    = max(Ctrlfreq_cmb{subjLoop, condLoop}.powspctrm(chan_pos_L,foi_pts), [], 2);
        Ctrl_freq_summary(subjLoop,condLoop)      = mean(temp_maxfreqCtrl);
        
        [temp_maxfreqCtrl_Ipsi, CtrlmaxFreq_Ipsi]    = max(Ctrlfreq_cmb{subjLoop, condLoop}.powspctrm(chan_pos_R,foi_pts), [], 2);
        Ctrl_freq_summary_Ipsi(subjLoop,condLoop)      = mean(temp_maxfreqCtrl_Ipsi);

    end
end

%
figure,
axis_set = []; %[0.5 7.5 -54.2 -52]
time_vector = [1 2 3 4 5 6 7];

subplot(2,2,1)
p1=plot(log(PD_freq_summary([5:10],:))', '-o'); hold on
errorbar(time_vector,mean(log(PD_freq_summary([5:10],:)))',sem(log(PD_freq_summary([5:10],:)))',...
    '-dr', 'linewidth', 4)
axis(axis_set), title('PD-Contra')
xlabel('Conditons'); ylabel('Power')

subplot(2,2,2)
plot(log((Ctrl_freq_summary))', '-o'); hold on
errorbar(time_vector,mean(log(Ctrl_freq_summary))',sem(log(PD_freq_summary))',...
    '-sb', 'linewidth', 4)
axis(axis_set), title('Control-Contra')
xlabel('Conditons'); ylabel('Power')

subplot(2,2,3)
plot(log((PD_freq_summary_Ipsi([5:10],:)))', '-o'); hold on
errorbar(time_vector,mean(log(PD_freq_summary_Ipsi([5:10],:)))',sem(log(PD_freq_summary_Ipsi([5:10],:)))',...
    '-dr', 'linewidth', 4)
axis(axis_set), title('PD-ipsi')
xlabel('Conditons'); ylabel('Power')

subplot(2,2,4)
plot(log((Ctrl_freq_summary_Ipsi))', '-o'); hold on
errorbar(time_vector,mean(log(Ctrl_freq_summary_Ipsi))',sem(log(Ctrl_freq_summary_Ipsi))',...
    '-sb', 'linewidth', 4)
axis(axis_set), title('Control-ipsi')
xlabel('Conditons'); ylabel('Power')

suptitle('Gamma-power-mean')

%%
export_fig( gcf, 'Gamma-power-mean','-transparent', ...
        '-painters','-pdf', '-r250' ); 

%% grand avg Beta power

cfg_m = []
cfg_m.operation = 'log10';
cfg_m.parameter = 'powspctrm';

cfg = []; 
cfg.parameter = 'powspctrm'
cfg.channel = 'EMG*';

for loop = 1:7
   
    PD_pow_grndAvg{loop} = ft_freqgrandaverage(cfg, PDfreq_cmb{[5:10],loop});
%     PD_pow_grnd_log{loop} = ft_math(cfg_m, PD_pow_grndAvg{loop});
    
%     Ctrl_pow_grndAvg{loop} = ft_freqgrandaverage(cfg, Ctrlfreq_cmb{:,loop});
%     Ctrl_pow_grnd_log{loop} = ft_math(cfg_m, Ctrl_pow_grndAvg{loop});
end

%%

for ch=1:10
    for cond = 1:7
%         ch_mat(ch, cond) =  ft_channelselection('EMG*', Ctrlfreq_cmb{ch,cond}.label);
        EMGpos = find(strcmp(PDfreq_cmb{ch,cond}.label,...
            ft_channelselection('EMG*', PDfreq_cmb{ch,cond}.label)))
        PDfreq_cmb{ch,cond}.label{EMGpos} = 'EMG003';
    end
end



%% MultiplotER 

cfg                 = [];
cfg.parameter       = 'powspctrm';
cfg.layout          = 'neuromag306cmb.lay';
cfg.graphcolor      = 'brkcmgy';
% cfg.channel         = R_par_ch;
% cfg.highlight = 'on'
% cfg.highlightchannel = [L_par_ch;R_par_ch];
% cfg.highlightsymbol    = 'x';
% cfg.highlightsize     = 12;
cfg.xlim            = [13 30];
cfg.zlim            = [2e-24 5e-24]
cfg.colorbar = 'yes';
figure
for loop=1:7
    subplot(3,3,loop)
    ft_topoplotER(cfg, Ctrl_pow_grndAvg{loop})
    title(num2str(loop))
end



%%
figure
cfg                 = [];
cfg.parameter       = 'powspctrm';
cfg.layout          = 'neuromag306cmb.lay';
cfg.graphcolor      = 'brkcmgy'
cfg.channel = {'all','-MEG1432+1433', '-MEG2622+2623','-MEG1332+1333', ...
    '-MEG1412+1413', '-MEG1422+1423', '-MEG1442+1443', '-MEG2612+2613',...
    '-MEG2632+2633', '-MEG2642+2643'}
cfg.highlight = 'on'
cfg.highlightchannel = [L_par_ch;R_par_ch];
cfg.highlightsymbol    = 'x';
cfg.highlightsize     = 12
cfg.xlim            = [13 30];
cfg.zlim            = [1e-24 5e-24]

for loop=1:7
    
    subplot(3,3,loop)
    ft_topoplotER(cfg, PD_pow_grndAvg{loop})
    title(num2str(loop)) 
    
end

L_par_ch_avg = {'MEG0212+0213', 'MEG0222+0223', 'MEG0232+0233', 'MEG0242+0243', ...
    'MEG0412+0413', 'MEG0442+0443', 'MEG1522+1523', 'MEG1612+1613',...
    'MEG1622+1623', 'MEG1812+1813'}

all_ch = {'MEG0212+0213', 'MEG0222+0223', 'MEG0232+0233', 'MEG0242+0243', ...
    'MEG0412+0413', 'MEG0442+0443', 'MEG1122+1123', 'MEG1132+1133', 'MEG1312+1313',...
    'MEG1322+1323', 'MEG1342+1343', 'MEG1522+1523', 'MEG1612+1613', 'MEG1622+1623',...
    'MEG1812+1813', 'MEG2222+2223', 'MEG2412+2413', 'MEG2422+2423'}

R_par_ch_avg = setdiff(all_ch, L_par_ch_avg);

%%
for subjLoop = 1:6
    hFig = figure(subjLoop)
    set(hFig, 'Position', [800 80 1000 968]);
    
    for condLoop = 1:7
        subplot(3,3,condLoop)
        ft_topoplotER(cfg, coh_cmb{subjLoop, condLoop})
        title(num2str(condLoop))

    end
    suptitle(coh_cmb{subjLoop, 1}.filename(1:7))
end

