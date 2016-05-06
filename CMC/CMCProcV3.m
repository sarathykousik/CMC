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

%% Loop over datasets

for loop=1:length(filenames)

    % Load data
    cfg                         = [];
    cfg.dataset                 = filenames(loop).name;
    cfg.channel                 = {'MEG', 'EMG'};
    data_import                 = ft_preprocessing(cfg);
    
    % Filter EMG
    cfg                         = [];
    cfg.lpfilter                = 'yes';
    cfg.lpfreq                  =  100;
    cfg.hpfilter                = 'yes';
    cfg.hpfreq                  =  1;
    cfg.dftfilter               = 'yes';
    cfg.channel                 =  {'EMG'}; 
    cfg.rectify                 = 'yes'
    cfg.detrend                 = 'yes';
    preproc_data_EMG            = ft_preprocessing(cfg, data_import);
%     plot(preproc_data_EMG.time{1}, preproc_data_EMG.trial{1})
    
    EMGpos = find(strcmp(data_import.label,ft_channelselection('EMG*', data_import.label)))
    data_import.trial{1}(EMGpos,:) = preproc_data_EMG.trial{1};

    % Epoch
    time_tot    = floor(data_import.time{1}(end));
    fs          = data_import.fsample;
    sample_tot  = time_tot*fs;
    trl         = [];
    trl         = [[0:fs:sample_tot-fs]' [fs:fs:sample_tot]'];
    trl(:,3)    = 0;%[fs.*size(trl,1)]';
    trl(:,1)    = trl(:,1)+1;

    cfg                         = [];
    cfg.channel                 = {'MEG'};
    cfg.trl                     = trl;
    data_epoched{loop}          = ft_redefinetrial(cfg, data_import);
    data_epoched{loop}.filename = filenames(loop).name;
    
end

data_epoched_trim = data_epoched;

% Combine datasets
cmb = [[2,3]; [8,9];[10,11]];

% cfg           = [];
% cfg.parameter = 'trial';

for cmbLoop=1:3
    data_epoched{cmb(cmbLoop,1)}= ...
        ft_appenddata(cfg, data_epoched{cmb(cmbLoop,1)},data_epoched{cmb(cmbLoop,2)})
    
    data_epoched{cmb(cmbLoop,2)} = {};
end

%%
for loop = 1:length(data_epoched)
   
    % Power calculation
    cfg                     = [];
    cfg.output              = 'pow';
    cfg.method              = 'mtmfft';
    cfg.taper               = 'hanning';
    cfg.foi                 = [1:2:100];
    cfg.keeptrials          = 'yes';
    freq{loop}              = ft_freqanalysis(cfg, data_epoched{loop});
    freq{loop}.filename     = data_epoched{loop}.filename;

    cfg                      = [];
    cfg.jackknife            = 'yes';
    freq_desc{loop}          = ft_freqdescriptives(cfg,  freq{loop});
    freq_desc{loop}.filename = data_epoched{loop}.filename;

    % Stability
    cfg                 = [];
    cfg.medianfilter    = 'yes';
    cfg.medianfiltord   = 50;
    cfg.channel         = 'EMG';
    EMG_ch              = ft_preprocessing(cfg, data_epoched{loop});
        for trialLoop=1:length(EMG_ch.trial)
            stab(trialLoop) = 1-(std(EMG_ch.trial{trialLoop})/mean(EMG_ch.trial{trialLoop}));
        end
    stability{loop}.stab           = stab;
    stability{loop}.filename       =  filenames(loop).name;
    stability_mean(loop)           = mean(stab);
    stability_std(loop)            = std(stab);
    
    % TFR
    cfg            = [];
    cfg.channel    = 'MEG';	                
    cfg.method     = 'wavelet';                
    cfg.width      = 7; 
    cfg.output     = 'pow';	
    cfg.foi        = 1:2:100;	                
    cfg.toi        = 0:.1:1;		              
    TFRwave{loop}  = ft_freqanalysis(cfg, data_epoched{loop});
   
end

%% Coh cmb

for loop = 1:9

    cfg             = [];
    cfg.method      = 'mtmfft';
    cfg.taper       = 'dpss';
    cfg.output      = 'powandcsd';
    cfg.keeptrials  = 'yes';
    cfg.jackknife   = 'yes';
    cfg.channelcmb  = {'MEG', 'EMG'};
    cfg.foi         = [1:45];
    cfg.tapsmofrq   = 5;
    freq_fourier    = ft_freqanalysis(cfg, data_scramble_control{loop});
    
    % CMC calc
    cfg                = [];
    cfg.method         = 'coh';
    cfg.channelcmb     = {'MEG', 'EMG'};
    cfg.jackknife      = 'yes';
    coh{loop}          = ft_connectivityanalysis(cfg, freq_fourier);
%     coh{loop}.filename = data_scramble_PD{loop}.filename;
    

    coh{loop}.label         = {coh{loop}.labelcmb{:,1}}';
    coh{loop}.powspctrm     = coh{loop}.cohspctrm;
    coh{loop}.dimord        = 'chan_freq';
    coh{loop}               = rmfield(coh{loop},'labelcmb');
    coh{loop}               = rmfield(coh{loop}, 'cohspctrm');
    coh_cmb{loop}           = ft_combineplanar([], coh{loop});


end


%% grand avg
% for loop = 1:9
%     coh_cmb_Ctrl_scr{loop}.idx = loop;
% end

% coh_cmb_res = reshape(coh_cmb_Ctrl_scr,3,3)'
cfg = []; 
cfg.parameter = 'powspctrm'

for loop = 1:7
   
    PD_CMC_grndAvg{loop}            = ft_freqgrandaverage(cfg, coh_cmb{[1:3],loop});
    PD_CMC_grndAvg{loop}.filename   = [coh_cmb{1,loop}.filename(1:7);  ...
                                       coh_cmb{2,loop}.filename(1:7);...
                                       coh_cmb{3,loop}.filename(1:7)];

    Ctrl_CMC_grndAvg{loop}          = ft_freqgrandaverage(cfg, coh_cmb{[4:6],loop});
    Ctrl_CMC_grndAvg{loop}.filename   = [coh_cmb{4,loop}.filename(1:7);  ...
                                       coh_cmb{5,loop}.filename(1:7);...
                                       coh_cmb{6,loop}.filename(1:7)];    
end

save Ctrl_CMC_grndAvg Ctrl_CMC_grndAvg -v7.3
save PD_CMC_grndAvg PD_CMC_grndAvg -v7.3

%% Get max freq/CMC


set_1 = {'MEG0412+0413', 'MEG0422+0423', 'MEG0432+0433', 'MEG0442+0443',...
                'MEG0632+0633', 'MEG0712+0713'}
set_2 = {'MEG0412+0413', 'MEG0422+0423', 'MEG0432+0433', 'MEG0442+0443', ...
                'MEG0632+0633', 'MEG0712+0713', 'MEG0742+0743'}
set_3 = {'MEG0432+0433', 'MEG0442+0443', 'MEG0712+0713', 'MEG0742+0743', ...
                'MEG1812+1813', 'MEG1822+1823'}
set_4 = {'MEG0432+0433', 'MEG0712+0713', 'MEG0742+0743', 'MEG1822+1823'};

L_par_ch = union(union(union(set_1, set_2),set_3),set_4);

%

cfg_sel      = [];
cfg_sel.frequency = 1;

chan_pos_L = match_str(Ctrl_CMC_grndAvg{1}.label, L_par_ch);
% chan_pos_R = match_str(coh_cmb{1}.label, R_par_ch);

for loop = 1:7
    
    
    [PDmaxCMC, PDmaxFreq]     = max(PD_CMC_grndAvg{loop}.powspctrm(chan_pos_L,13:30),...
                                                        [], 2);
     PDCMCsummary(loop)      = max(PDmaxCMC);
    
    [Ctrl_maxCMC, PDmaxFreq]     = max(Ctrl_CMC_grndAvg{loop}.powspctrm(chan_pos_L,...
                                                        13:30), [], 2);
     CtrlCMCsummary(loop)      = max(Ctrl_maxCMC);
    
    
end


% maxCMC_res = reshape(maxCMCmax',3,[])'
% maxCMC_summary = reshape(atanh(CtrlCMCsummary)',3,[])'
% maxCMC_ctrl_hann_res = reshape(atanh(maxCMCmax_Ctrl)',3,[])'

%
figure,
plot(PDCMCsummary, '-or'); hold on
plot(CtrlCMCsummary, '-db')
legend('PD', 'Ctrl')
axis([0.5 7.5 0.04 0.2])
title('CMC summary values from grand averages (max)')
ylabel('CMC'); xlabel('Conditions')

%%
export_fig( gcf, 'topos-Control','-transparent', ...
        '-painters','-pdf', '-r250' ); 


%%
figure
wat = [1:3];
cfg                 = [];
cfg.parameter       = 'powspctrm';
cfg.layout          = 'neuromag306cmb.lay';
cfg.graphcolor      = 'brkcmgy'
% cfg.channel         = 'MEGMAG'
% cfg.xlim            = [1 99];
cfg.xlim            = [13 30];
cfg.zlim = [0.04 0.1];
% ft_multiplotER(cfg, ft_combineplanar([],freq_fourier))

for loop=1:7
    subplot(3,3,loop)
    ft_topoplotER(cfg, Ctrl_CMC_grndAvg{loop})
    title(num2str(loop))
end
suptitle('Control')


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



%% Scramble CMC
data_scramble = data_control;
cfg_sel = [];

for subjLoop = 1:9
    for trialLoop = 1:length(data_control{1,subjLoop}.trial)-3
        
        data_scramble{1,subjLoop}.trial{trialLoop}(1,:) = data_scramble{1,subjLoop}.trial{trialLoop+3}(1,:);
        
    end
    cfg_sel.trials = 1:trialLoop;
    data_scramble{1,subjLoop} = ft_selectdata(cfg_sel, data_scramble{1,subjLoop});
    
end









