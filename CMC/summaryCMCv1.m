
%% Setup folders and files

% restoredefaultpath
% addpath('/usr/local/common/matlab_toolbox/fieldtrip/r6152')
%addpath('C:\Program Files\MATLAB\R2012b\toolbox_add_on\export_fig')
addpath('C:\Program Files\MATLAB\R2012b\toolbox_add_on\cmc\')
%addpath('C:\Program Files\MATLAB\R2012b\toolbox_add_on\')
ft_defaults
 
proc.data_folder_control    = 'J:\MEG_Research\CMC\control_move_artRej';
proc.data_folder_PD    = 'J:\MEG_Research\CMC\PD_move_artRej';
proc.results_folder = 'J:\MEG_Research\CMC\plots\';

mkdir(proc.results_folder)
% cd(proc.data_folder_PD)
filenames   = dir('*-raw.fif');

for loop = 1: length(filenames)
    file_sub(loop) = {filenames(loop).name};
end

file_sub    = reshape(file_sub, 7, [])';
[row col]   = size(file_sub);
cmc

save_folder='j:/meg_research/cmc/plots/';

%% Spatial norm of coh values=> (X-mean)./std

for subjLoop=1:6
    for condLoop=1:3
%         coh_cmb{subjLoop,condLoop}.powspctrm_z = atanh(coh_cmb{subjLoop,condLoop}.powspctrm);
        condLoop
    coh_beta_cmb{subjLoop,condLoop}.powspctrm_z = atanh(coh_beta_cmb{subjLoop,condLoop}.powspctrm);
%     PDcoh_beta_cmb{subjLoop,condLoop}.powspctrm_z = atanh(coh_beta_cmb{subjLoop,condLoop}.powspctrm);
%     Ctrlcoh_cmb{subjLoop,condLoop}.powspctrm_z = atanh(Ctrlcoh_cmb{subjLoop,condLoop}.powspctrm);
%     PDcoh_cmb{subjLoop,condLoop}.powspctrm_z = atanh(PDcoh_cmb{subjLoop,condLoop}.powspctrm);
% Ctrlcoh_beta_cmb{subjLoop,condLoop}.powspctrm_z = atanh(Ctrlcoh_beta_cmb{subjLoop,condLoop}.powspctrm);
%     PDcoh_beta_cmb{subjLoop,condLoop}.powspctrm_z = atanh(PDcoh_beta_cmb{subjLoop,condLoop}.powspctrm);
% Ctrlcoh_beta_cmb{subjLoop,condLoop}.powspctrmnorm = ...
%     (Ctrlcoh_beta_cmb{subjLoop,condLoop}.powspctrm- ...
%         mean(mean(Ctrlcoh_beta_cmb{subjLoop,condLoop}.powspctrm)))./ ...
%          std(Ctrlcoh_beta_cmb{subjLoop,condLoop}.powspctrm);
%     PDcoh_beta_cmb{subjLoop,condLoop}.powspctrmnorm = ...
%         (PDcoh_beta_cmb{subjLoop,condLoop}.powspctrm- ...
%             mean(mean(PDcoh_beta_cmb{subjLoop,condLoop}.powspctrm)))./ ...
%              std(PDcoh_beta_cmb{subjLoop,condLoop}.powspctrm);        
    end
end
    
%% grand avg CMC/power

cfg = []; 
cfg.parameter = 'powspctrm_z'

for loop = 1:3
   
%     PD_pow_grndAvg_HG{loop} = ft_freqgrandaverage(cfg, PDfreq_cmb{5:10,loop});
%     PD_pow_grndAvg_FT{loop} = ft_freqgrandaverage(cfg, PDfreq_cmb{1:4,loop});
%     PD_pow_grndAvg_all{loop} = ft_freqgrandaverage(cfg, PDfreq_cmb{:,loop});
%     Ctrl_pow_grndAvg{loop} = ft_freqgrandaverage(cfg, Ctrlfreq_cmb{:,loop});
     
%     PD_freq_beta_GA{loop} = ft_freqgrandaverage(cfg, PDfreq_beta_cmb{:,loop});
%     if loop==1
%         PD_CMC_grndAvg_HG{loop}  = ft_freqgrandaverage(cfg, PDcoh_cmb{[6 7 8],loop});
%         PD_CMC_beta_grndAvg_z{loop} = ft_freqgrandaverage(cfg,PDcoh_beta_cmb{6:10,loop});

%     else
%         PD_CMC_grndAvg_HG{loop}  = ft_freqgrandaverage(cfg, PDcoh_cmb{[6 7 8],loop});
%         PD_CMC_beta_grndAvg_z{loop} = ft_freqgrandaverage(cfg,PDcoh_beta_cmb{5:10,loop});
%     end
%     PD_CMC_grndAvg_FT{loop}  = ft_freqgrandaverage(cfg, PDcoh_cmb{1:4,loop});
%     PD_CMC_grndAvg_all{loop} = ft_freqgrandaverage(cfg, PDcoh_cmb{:,loop});
    
%     PD_CMC_alpha_grndAvg_all{loop} = ft_freqgrandaverage(cfg,PDcoh_alpha_cmb{:,loop});
    PD_CMC_beta_grndAvg_all{loop} = ft_freqgrandaverage(cfg,coh_beta_cmb{:,loop});
%     PD_CMC_gamma_grndAvg_all{loop} = ft_freqgrandaverage(cfg,PDcoh_gamma_cmb{:,loop});
%     PD_CMC_grndAvg{loop} = ft_freqgrandaverage(cfg, coh_cmb{[1 2 3],loop});
%     Ctrl_CMC_grndAvg{loop} = ft_freqgrandaverage(cfg, coh_cmb{[4 5 6],loop});
%     Ctrl_CMC_alpha_grndAvg_all{loop} = ft_freqgrandaverage(cfg,Ctrlcoh_alpha_cmb{:,loop});
%     Ctrl_CMC_beta_grndAvg{loop} = ft_freqgrandaverage(cfg,Ctrlcoh_beta_cmb{:,loop});
%     Ctrl_CMC_gamma_grndAvg_all{loop} = ft_freqgrandaverage(cfg,Ctrlcoh_gamma_cmb{:,loop});
    
    
end

%%
cfg_cmb = []; cfg_cmb.combinemethod = 'cmc';

for subjLoop=1:10
    for condLoop=1:7
   
        cfg                      = [];
        cfg.jackknife            = 'yes';
        freq_desc                = ft_freqdescriptives(cfg,  PDfreq_beta{subjLoop,condLoop});
        PDfreq_beta{subjLoop, condLoop}
        PDfreq_beta_cmb{subjLoop, condLoop} = ft_combineplanar_cmc(cfg_cmb,freq_desc);
        
    end
end



%% GGA
cfg = []; 
cfg.parameter = 'powspctrmnorm'

for loop = 1:7

    GGA_CMC_beta{loop} = ft_freqgrandaverage(cfg,Ctrl_CMC_beta_grndAvg_all{loop},...
                                                PD_CMC_beta_grndAvg_all{loop});
    
end

%%
save PD_pow_grndAvg_HG PD_pow_grndAvg_HG -v7.3
save PD_pow_grndAvg_FT PD_pow_grndAvg_FT -v7.3
save PD_pow_grndAvg_all PD_pow_grndAvg_all -v7.3
save Ctrl_pow_grndAvg Ctrl_pow_grndAvg -v7.3

save PD_CMC_grndAvg_HG PD_CMC_grndAvg_HG -v7.3
save PD_CMC_grndAvg_FT PD_CMC_grndAvg_FT -v7.3
save PD_CMC_grndAvg_all PD_CMC_grndAvg_all -v7.3
save Ctrl_CMC_grndAvg Ctrl_CMC_grndAvg -v7.3


%% SOI 
% Left
L_par_ch = {'MEG0412+0413', 'MEG0422+0423', 'MEG0432+0433', 'MEG0442+0443',...
    'MEG0632+0633', 'MEG0712+0713', 'MEG0742+0743', 'MEG1812+1813', 'MEG1822+1823'};

% Right
R_par_ch = {'MEG0722+0723', 'MEG0732+0733', 'MEG1042+1043', 'MEG1112+1113',...
    'MEG1122+1123', 'MEG1132+1133', 'MEG1142+1143', 'MEG2212+2213', 'MEG2222+2223'};

%% Plot chan pos
% sample = Ctrlcoh_cmb{1,1};
% sample.powspctrm=sample.powspctrm.*0;
figure
cfg                     = [];
cfg.parameter           = 'powspctrm_z';
cfg.layout              = 'neuromag306cmb.lay';
cfg.highlight           = 'on'
cfg.highlightchannel    = [L_par_ch];
cfg.highlightsymbol     = 'x';
cfg.highlightsize       = 12;
cfg.zlim                = [0.02 0.08]
% cfg.colorbar            = 'yes';
for loop=1:3
    subplot(1,3,loop)
    ft_topoplotER(cfg, PD_CMC_grndAvg{loop})
end
% ft_topoplotER(cfg,coh_cmb{1,:})

%%
loop=1
for subjLoop=4:6
    for condLoop=1:3
        subplot(3,3,loop)
        ft_topoplotER(cfg, coh_cmb{subjLoop, condLoop})
        title(coh_cmb{subjLoop, condLoop}.filename)
        loop=loop+1;
    end
end

%%

cfg = [];
 cfg.channel = {'MEG1132','MEG1133', 'MEG1312','MEG1313', 'MEG1322','MEG1323',...
 'MEG1332','MEG1333','MEG1342','MEG1343', 'MEG1432','MEG1433', 'MEG1442','MEG1443', ...
'MEG2412','MEG2413', 'MEG2422','MEG2423', 'MEG2432','MEG2433', 'MEG2442','MEG2443',...
'MEG2522','MEG2523', 'MEG2612','MEG2613', 'MEG2622','MEG2623', 'MEG2632','MEG2633',...
'MEG2642','MEG2643'}
% cfg.viewmode = 'vertical';
ft_databrowser(cfg,data_import)

%% Power summary values from individual 

cfg_sel      = [];
cfg_sel.frequency = 1;

chan_pos_L = match_str(PDfreq_cmb{1}.label, L_par_ch);
chan_pos_R = match_str(PDfreq_cmb{1}.label, R_par_ch);

foi=[31,45]; title_str='Gamma-FT-power';
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
figure,
axis_set = [0.5 7.5 -54.8 -53];%[0.5 7.5 -53.5 -51.2]
time_vector = [1 2 3 4 5 6 7];
PD_subjs = [1:4];
subplot(2,2,1)
p1=plot(log(PD_freq_summary(PD_subjs,:))', '-o'); hold on
errorbar(time_vector,mean(log(PD_freq_summary(PD_subjs,:)))',sem(log(PD_freq_summary(PD_subjs,:)))',...
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
plot(log((PD_freq_summary_Ipsi(PD_subjs,:)))', '-o'); hold on
errorbar(time_vector,mean(log(PD_freq_summary_Ipsi(PD_subjs,:)))',sem(log(PD_freq_summary_Ipsi(PD_subjs,:)))',...
    '-dr', 'linewidth', 4)
axis(axis_set), title('PD-ipsi')
xlabel('Conditons'); ylabel('Power')

subplot(2,2,4)
plot(log((Ctrl_freq_summary_Ipsi))', '-o'); hold on
errorbar(time_vector,mean(log(Ctrl_freq_summary_Ipsi))',sem(log(Ctrl_freq_summary_Ipsi))',...
    '-sb', 'linewidth', 4)
axis(axis_set), title('Control-ipsi')
xlabel('Conditons'); ylabel('Power')

suptitle(title_str)

%%
export_fig(gcf,['j:/meg_research/cmc/plots/','PD-HG-only-beta-mtm-GA'] ,...
         '-transparent', '-painters','-pdf', '-r250' )

%% Ch from GGA: focal
cond_1 = {'MEG0432+0433', 'MEG0442+0443','MEG0232+0233', 'MEG1622+1623'}
cond_2 = {'MEG0432+0433', 'MEG0442+0443', 'MEG1812+1813','MEG0232+0233'}
cond_3 = {'MEG0432+0433', 'MEG0442+0443', 'MEG1812+1813'}
cond_4 = {'MEG0432+0433', 'MEG0442+0443', 'MEG1812+1813'}
cond_5 = {'MEG0432+0433', 'MEG0712+0713'}
cond_6 = {'MEG0432+0433', 'MEG0442+0443', 'MEG1812+1813', 'MEG1822+1823'}
cond_7 = {'MEG0432+0433', 'MEG0442+0443'}

ch_GGA = unique(horzcat(cond_1,cond_2,cond_3,cond_4,cond_5,cond_6,cond_7))

%% CMC Summary from individual

cfg_sel      = [];
cfg_sel.frequency = 1;

% chan_pos_L = match_str(coh_cmb{1,1}.label, ch_GGA);
chan_pos_L = match_str(coh_cmb{1}.label, L_par_ch);

freq_range = 'beta'
foi=[13,30]; supText = [freq_range,'FT-CMC-mtm '];
foi_pts = find(coh_cmb{1,1}.freq>=foi(1) & coh_cmb{1,1}.freq<=foi(2));

for subjLoop = 1:6
    for condLoop = 1:3
%              CMC_summary(subjLoop,condLoop)    = ...
%              eval(['mean(coh_', freq_range,'_cmb{subjLoop, condLoop}.powspctrm(chan_pos_L));']);
    
        [temp_maxCMCPD, PDmaxFreq]    = max(coh_cmb{subjLoop, condLoop}.powspctrm(chan_pos_L,foi_pts), [], 2);
        CMC_summary(subjLoop,condLoop)      = mean(temp_maxCMCPD);
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
%         PD_CMC_summary(subjLoop,condLoop)    = ...
%             eval(['mean(PDcoh_', freq_range,'_cmb{subjLoop, condLoop}.powspctrm_z(chan_pos_L));']);
%         PD_CMC_summary_Ipsi(subjLoop,condLoop) =...
%             eval(['mean(PDcoh_', freq_range,'_cmb{subjLoop, condLoop}.powspctrm(chan_pos_R));']);
%         Ctrl_CMC_summary(subjLoop,condLoop)      = ...
%             eval(['mean(Ctrlcoh_', freq_range,'_cmb{subjLoop, condLoop}.powspctrm_z(chan_pos_L));']);
%         Ctrl_CMC_summary_Ipsi(subjLoop,condLoop)      = ...
%             eval(['mean(Ctrlcoh_', freq_range,'_cmb{subjLoop, condLoop}.powspctrm(chan_pos_R));']);

    end
end

%%
axis_set = [0.5 3.5 0.02 0.2]
% axis_set = [0.5 7.5 0 1.5]
conds = [1:3];
PD_subjs = [1:6];

% subplot(2,1,1)
p1=plot((PD_CMC_summary(1:3,conds))', '-o'); hold on
errorbar(mean(PD_CMC_summary(1:3,conds))',sem(PD_CMC_summary(PD_subjs,conds))', '-or', 'linewidth', 4);
p2=plot((PD_CMC_summary(4:6,conds))', '-o'); hold on
errorbar(mean(PD_CMC_summary(4:6,conds))',sem(PD_CMC_summary(PD_subjs,conds))', '-ob', 'linewidth', 4);
axis(axis_set), title('PD-Contra-mean');
xlabel('Conditons'); ylabel('CMC')

subplot(2,1,2)
p2=plot((Ctrl_CMC_summary(:,conds))', '-o'); hold on
errorbar(mean(Ctrl_CMC_summary(:,conds))',sem(Ctrl_CMC_summary(:,conds)) , '-ob', 'linewidth', 4);
axis(axis_set), title('Control-Contra-mean');
xlabel('Conditons'); ylabel('CMC')

% subplot(2,2,3)
% p2=plot(PD_CMC_summary_Ipsi(PD_subjs,conds)', '-o'); hold on
% errorbar(mean(PD_CMC_summary_Ipsi(PD_subjs,conds))',sem(PD_CMC_summary_Ipsi(PD_subjs,conds)),'-or', 'linewidth', 4);
% axis(axis_set), title('PD-ipsi');
% xlabel('Conditons'); ylabel('CMC')
% 
% subplot(2,2,4)
% p2=plot((Ctrl_CMC_summary_Ipsi(:,conds))', '-o'); hold on
% errorbar(mean(Ctrl_CMC_summary_Ipsi(:,conds))',sem(Ctrl_CMC_summary_Ipsi(:,conds)),  '-ob', 'linewidth', 4);
% axis(axis_set), title('Control-ipsi');
% xlabel('Conditons'); ylabel('CMC')
suptitle([supText,num2str(PD_subjs)])

%%
export_fig(gcf,[save_folder,supText,'-focalSOI'] ,...
         '-transparent', '-painters','-pdf', '-r250' )


%% CMC Summary from individual - MEAN diff

cfg_sel      = [];
cfg_sel.frequency = 1;

chan_pos_L = match_str(PDcoh_cmb{1}.label, L_par_ch);
chan_pos_R = match_str(PDcoh_cmb{1}.label, R_par_ch);

foi=[13,30]; supText = 'Beta-CMC';
foi_pts = find(Ctrlcoh_cmb{1,1}.freq>=foi(1) & Ctrlcoh_cmb{1,1}.freq<=foi(2));

for subjLoop = 1:10
    for condLoop = 1:7
    
        [temp_maxCMCPD]    = mean(PDcoh_cmb{subjLoop, condLoop}.powspctrm(chan_pos_L,foi_pts),1);%, [], 2);
        PD_CMC_summary(subjLoop,condLoop)      = max(temp_maxCMCPD);
        
        [temp_maxCMCPD_Ipsi]    = mean(PDcoh_cmb{subjLoop, condLoop}.powspctrm(chan_pos_R,foi_pts),1);%, [], 2);
        PD_CMC_summary_Ipsi(subjLoop,condLoop)      = max(temp_maxCMCPD_Ipsi);

        [temp_maxCMCCtrl]    = mean(Ctrlcoh_cmb{subjLoop, condLoop}.powspctrm(chan_pos_L,foi_pts),1);%, [], 2);
        Ctrl_CMC_summary(subjLoop,condLoop)      = max(temp_maxCMCCtrl);
        
        [temp_maxCMCCtrl_Ipsi]    = mean(Ctrlcoh_cmb{subjLoop, condLoop}.powspctrm(chan_pos_R,foi_pts),1);%, [], 2);
        Ctrl_CMC_summary_Ipsi(subjLoop,condLoop)      = max(temp_maxCMCCtrl_Ipsi);

    end
end


%%
axis_set = [0.5 4.5 0.04 0.15]
conds = [1 4 6 7];
PD_subjs = [2:4 6:10];
figure,
subplot(2,2,1)
p1=plot((PD_CMC_summary(PD_subjs,conds))', '-o'); hold on
errorbar(mean(PD_CMC_summary(PD_subjs,conds))',sem(PD_CMC_summary(PD_subjs,conds))', '-or', 'linewidth', 4)
axis(axis_set), title('PD-Contra')
xlabel('Conditons'); ylabel('CMC')

subplot(2,2,2)
p2=plot((Ctrl_CMC_summary(:,conds))', '-o'); hold on
errorbar(mean(Ctrl_CMC_summary(:,conds))',sem(Ctrl_CMC_summary(:,conds)) , '-ob', 'linewidth', 4)
axis(axis_set), title('Control-Contra')
xlabel('Conditons'); ylabel('CMC')

subplot(2,2,3)
p2=plot(PD_CMC_summary_Ipsi(PD_subjs,conds)', '-o'); hold on
errorbar(mean(PD_CMC_summary_Ipsi(PD_subjs,conds))',sem(PD_CMC_summary_Ipsi(PD_subjs,conds)),'-or', 'linewidth', 4)
axis(axis_set), title('PD-ipsi')
xlabel('Conditons'); ylabel('CMC')

subplot(2,2,4)
p2=plot((Ctrl_CMC_summary_Ipsi(:,conds))', '-o'); hold on
errorbar(mean(Ctrl_CMC_summary_Ipsi(:,conds))',sem(Ctrl_CMC_summary_Ipsi(:,conds)),  '-ob', 'linewidth', 4)
axis(axis_set), title('Control-ipsi')
xlabel('Conditons'); ylabel('CMC')
suptitle([supText,num2str(PD_subjs)])

%% Stats-pow
clc
conds = [1:7];
PD_subjs = [1:10];

cond_cmbs = nchoosek(conds,2);

for loop=1:size(cond_cmbs,1)
   
    [h(loop) p(loop) ci(loop,:) ~] = ttest2(PD_freq_summary(PD_subjs, cond_cmbs(loop,1))-PD_freq_summary(PD_subjs, cond_cmbs(loop,2)),...
        Ctrl_freq_summary(:, cond_cmbs(loop,1))-Ctrl_freq_summary(:, cond_cmbs(loop,2)));

end
cond_cmbs(find(p<0.05),:)
[p; [1:21]]

%% Stats CMC

conds = [1 6 7];
PD_subjs = [6:10];

cond_cmbs = nchoosek(conds,2);
p=[]; h=[], ci=[];
for loop=1:size(cond_cmbs,1)
   %interaction
%     [h(loop) p(loop) ci(loop,:) ~] = ttest2(PD_CMC_summary(PD_subjs, cond_cmbs(loop,1))-...
%         PD_CMC_summary(PD_subjs, cond_cmbs(loop,2)),...
%         Ctrl_CMC_summary(:, cond_cmbs(loop,1))-Ctrl_CMC_summary(:, cond_cmbs(loop,2)));
    %within PD
        [h(loop) p(loop) ci(loop,:) ~] = ttest(PD_CMC_summary(PD_subjs, cond_cmbs(loop,1))-...
        PD_CMC_summary(PD_subjs, cond_cmbs(loop,2)));
    
end
find(p<0.05)
cond_cmbs(find(p<0.1),:)
[p; [1:length(cond_cmbs)]]

%% CMC
%Results - interaction (22-2-2016: 11:19) - with 1s epochs
p_HG = [0.9437    0.3543    0.1184]
p_FT = [0.8408    0.6133    0.7274]
p_all= [0.8701    0.6235    0.3809]

% Results - within group (22-2-2016: 11:19) - with 1s epochs
PD_HG  = [0.5255    0.2537    0.1249]
PD_FT  = [0.7465    0.9447    0.7112]
PD_all = [0.4476    0.2635    0.3121]

ctrl_all = [0.3875    0.4339    0.9323]

%% Topo Coh - Alpha/beta/gamma

cfg                 = [];
cfg.parameter       = 'powspctrm';
cfg.layout          = 'neuromag306cmb.lay';
cfg.graphcolor      = 'brkcmgy';
cfg.channel         = {'MEG*', '-MEG1522+1523', '-MEG1612+1613'};
% cfg.channel         = R_par_ch;
% cfg.highlight = 'on'
% cfg.highlightchannel = [ch_GGA];
% cfg.highlightcolor = [0.7 0.5 0.7];
% cfg.highlightsymbol    = '*';
% cfg.highlightsize     = 8;
% cfg.highlightfontsize  =8;
% cfg.xlim            = [55 95];
% cfg.zlim = [1e-23 4e-23] %alpha
cfg.zlim = [2e-24 10e-24] %beta
% cfg.zlim = [2.5e-24 5e-23] %gamma
% cfg.zlim = [1e-24 4e-24] %High gamma
% cfg.zlim            = [0 2]
cfg.colorbar = 'yes';
loop=1;
hFig=figure;
set(hFig, 'Position', [100 50 1556 1000]);

% idx=[1 6 7 ];

%  for subjLoop=1:10
     
     for condLoop= 1:7
         subplot(3,3,condLoop)
%          subplot(3,10,subjLoop+(condLoop-1)*10)
%          sample=PDcoh_beta_cmb{subjLoop,idx(condLoop)};
%          sample.powspctrm = ztrans(sample.powspctrm);
%          ft_topoplotER(cfg, sample)
         ft_topoplotER(cfg, PD_freq_beta_GA{condLoop})
         

         title(num2str(condLoop))
%          title(file_sub{subjLoop,condLoop}(15:16))
%      end
%      loop=loop+1;

%     suptitle('PD-beta-CMC')
%     export_fig( gcf,[save_folder,file_sub{subjLoop,1}(1:21),'-alphaCMCTopo'] ,...
%         '-transparent', '-painters','-pdf', '-r250' );


  end
suptitle('GGA-beta-CMC-Norm')
%%
export_fig(gcf,[save_folder,'Ctrl-beta-CMC-norm'] ,...
         '-transparent', '-painters','-pdf', '-r250' )

%% Ch from GGA: focal
cond_1 = {'MEG0432+0433', 'MEG0442+0443','MEG0232+0233', 'MEG1622+1623'}
cond_2 = {'MEG0432+0433', 'MEG0442+0443', 'MEG1812+1813','MEG0232+0233'}
cond_3 = {'MEG0432+0433', 'MEG0442+0443', 'MEG1812+1813'}
cond_4 = {'MEG0432+0433', 'MEG0442+0443', 'MEG1812+1813'}
cond_5 = {'MEG0432+0433', 'MEG0712+0713'}
cond_6 = {'MEG0432+0433', 'MEG0442+0443', 'MEG1812+1813', 'MEG1822+1823'}
cond_7 = {'MEG0432+0433', 'MEG0442+0443'}

ch_GGA = unique(horzcat(cond_1,cond_2,cond_3,cond_4,cond_5,cond_6,cond_7))

%%
cfg                 = [];
cfg.parameter       = 'powspctrm';
cfg.layout          = 'neuromag306cmb.lay';
cfg.graphcolor      = 'brykcgm';
cfg.channel         = L_par_ch;
% cfg.highlight = 'on'
% cfg.highlightchannel = [L_par_ch;R_par_ch];
% cfg.highlightsymbol    = 'x';
% cfg.highlightsize     = 12;
cfg.xlim            = [1 100];
%cfg.zlim            = [2e-24 5e-24]
cfg.colorbar = 'yes';

for subjLoop=1:10
    
    hFig=figure(subjLoop);
    set(hFig, 'Position', [100 50 1556 1000]);
    
%      for condLoop=1:7
%          subplot(3,3,condLoop)
         ft_singleplotER(cfg, PDcoh_cmb{subjLoop, condLoop})
         title(file_sub{subjLoop,condLoop}(15:16))
%      end
    
    suptitle(num2str(subjLoop))
    export_fig( gcf,[save_folder,file_sub{subjLoop,1}(1:21),'-CMCmulti'] ,...
        '-transparent', '-painters','-pdf', '-r250' ); 
end

%%
cfg_neighb.method       = 'distance';
neighbours          = ft_prepare_neighbours(cfg_neighb, PDcoh_cmb{1,1});



%% Setup interaction data


conds = [1,7]
cfg_math=[]; cfg_math.operation='subtract';
cfg_math.parameter='powspctrm';

for subjLoop=1:10
    for condLoop=conds
        
        diff_PD{subjLoop} = ft_math(cfg_math, PDcoh_beta_cmb{subjLoop,conds(1)},...
                                        PDcoh_beta_cmb{subjLoop,conds(2)});
        diff_C{subjLoop} = ft_math(cfg_math, Ctrlcoh_beta_cmb{subjLoop,conds(1)},...
                                        Ctrlcoh_beta_cmb{subjLoop,conds(2)});
    
% diff_Effect{1}=cfg_math(diff_PD, diff_C);


    end
end

% diff_Effect=cfg_math(diff_PD, diff_C);

%% Design

% design=[];
% ntrials = size(PDfreq_fourier_beta{1,1}.powspctrm,1);
% design  = zeros(2,2*ntrials);
% design(1,1:ntrials) = 1;
% design(1,ntrials+1:2*ntrials) = 2;
% design(2,1:ntrials) = [1:ntrials];
% design(2,ntrials+1:2*ntrials) = [1:ntrials];

conds = [1,7]
design=[];
total=0;

for subjLoop=1%:10
%    for condLoop=1:7
       ntrials = size(freq_fourier_beta{1,conds(1)}.fourierspctrm,1) ;% size(PD_fourier_beta_cmb_1_7{subjLoop,conds(1)}.powspctrm,1);
       design(1,total+1:ntrials+total) = 1;
       design(2,1+total:ntrials+total) = [1:ntrials];
       total=total+ntrials;
%    end
end

for subjLoop=1%:10
%    for condLoop=1:7
       ntrials = size(freq_fourier_beta{1,conds(2)}.fourierspctrm,1) ;% size( PD_fourier_beta_cmb_1_7{subjLoop,conds(2)}.powspctrm,1);
       design(1,total+1:ntrials+total) = 2;
       design(2,1+total:ntrials+total) = [1:ntrials];
       total=total+ntrials;
%    end
end


plot(design(end,:));figure(gcf);
hold on,plot(design(1,:)*50,'r');figure(gcf);

%% Non-parametric stats
cfg                     = [];
% cfg.latency             = [0 0.1];
% cfg.frequency           = [13:30];
cfg.parameter           = 'fourierspctrm';
% cfg.channel             = {'MEG'};%, '-MEG1522+1523', '-MEG1612+1613'};
% cfg.channel = {'MEG0232+0233', 'MEG0412+0413', 'MEG0422+0423', 'MEG0432+0433',...
%     'MEG0632+0633', 'MEG0712+0713', 'MEG0742+0743', 'MEG1612+1613', 'MEG1622+1623',...
%     'MEG1632+1633', 'MEG1812+1813', 'MEG1822+1823', 'MEG1832+1833', 'MEG1842+1843'};
cfg.statistic           = 'ft_statfun_indepsamplesZcoh';%'ft_statfun_depsamplesT'; %
% cfg.statistic           = 'ft_statfun_depsamplesFunivariate';
% cfg.statistic           = 'ft_statfun_depsamplesregrT';
% cfg.method              = 'analytic';
% cfg.correctm            = 'no';
cfg.method              = 'montecarlo';
cfg.correctm            = 'fdr';
%        cfg.channelcmb  = {'MEG', 'EMG'};
cfg.clusterstatistic    = 'maxsum';
cfg.minnbchan           = 2;
cfg.tail                = 0;
cfg.clustertail         = 0;
cfg.alpha               = 0.05;
cfg.clusteralpha        = 0.05;
cfg.numrandomization    = 500;
cfg_neighb.method       = 'distance';
cfg.neighbours          = ft_prepare_neighbours(cfg_neighb, freq_fourier_beta{1,1});

% cfg.neighbours = neighbours;
% 
% design=[];
% ntrials = size(PDfreq_fourier_beta{1,1}.powspctrm,1);
% design  = zeros(2,2*ntrials);
% design(1,1:ntrials) = 1;
% design(1,ntrials+1:2*ntrials) = 2;
% design(2,1:ntrials) = [1:ntrials];
% design(2,ntrials+1:2*ntrials) = [1:ntrials];


subj                    = [1 2 3:10];
nsubj                   = length(subj);
% design                  = [];
% design(1,:)             = [ones(1,nsubj), 2*ones(1,nsubj)];% 3*ones(1,nsubj)...
% %                                 4*ones(1,nsubj) 5*ones(1,nsubj) 6*ones(1,nsubj)];
% design(2,:)             = [1:nsubj 1:nsubj];% 1:nsubj 1:nsubj 1:nsubj 1:nsubj ];
cfg.design              = design;
cfg.ivar                = 1;
%  cfg.uvar               = 2;
stat                    = ft_freqstatistics(cfg,freq_fourier_beta{1,1},...
                                freq_fourier_beta{1,7});

figure
plot(stat.prob)
figure
plot(stat.mask)
%%
stat.raweffect = PD_fourier_beta_cmb_1_7{1,1}.fourierspctrm-PD_fourier_beta_cmb_1_7{1,1}.fourierspctrm;

cfg = [];
cfg.alpha  = 0.05;
cfg.parameter = 'raweffect';
% cfg.zlim   = [-1e-27 1e-27];
cfg.layout              = 'neuromag306cmb.lay';
ft_clusterplot(cfg, stat);

%%

freq               = [];
freq.label         = {'1';'2'};
freq.freq          = 1:10;
freq.fourierspctrm = randn(20,2)+randn(20,2).*1i;
freq.cumtapcnt     = 2*ones(10,1);
freq.dimord        = 'rpttap_chan_freq';

% create a configuration for the statistics
cfg           = [];
cfg.design    = [ones(1,10) ones(1,10)*2];
cfg.parameter = 'fourierspctrm';
cfg.statistic = 'ft_statfun_indepsamplesZcoh';
% cfg.correctm = 'cluster';
cfg.method    = 'montecarlo';
cfg.numrandomization = 1;
cfg.label     = freq.label;
stat          = ft_freqstatistics(cfg, freq);





