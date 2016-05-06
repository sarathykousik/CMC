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

for loop = 1:18

    cfg             = [];
    cfg.method      = 'mtmfft';
    cfg.taper       = 'dpss';
    cfg.output      = 'powandcsd';
    cfg.keeptrials  = 'yes';
    cfg.jackknife   = 'yes';
    cfg.channelcmb  = {'MEG', 'EMG'};
    cfg.foi         = [1:45];
    cfg.tapsmofrq   = 5;
    freq_fourier    = ft_freqanalysis(cfg, data_epoched{loop});
    
    % CMC calc
    cfg                = [];
    cfg.method         = 'coh';
    cfg.channelcmb     = {'MEG', 'EMG'};
    cfg.jackknife      = 'yes';
    coh{loop}          = ft_connectivityanalysis(cfg, freq_fourier);
%     coh{loop}.filename = data_scramble_PD{loop}.filename;
    
%     % Pow corr
%     cfg.method                     = 'powcorr';
%     powCorrSpc{loop}               = ft_connectivityanalysis(cfg, freq_fourier);
%     powCorrSpc{loop}.filename      = data_epoched{loop}.filename;
%     
%     % Pow corr ortho
%     cfg.method              = 'plv';
%     plv{loop}               = ft_connectivityanalysis(cfg, freq_fourier);
%     plv{loop}.filename      = data_epoched{loop}.filename;
    
    coh{loop}.label         = {coh{loop}.labelcmb{:,1}}';
    coh{loop}.powspctrm     = coh{loop}.cohspctrm;
    coh{loop}.dimord        = 'chan_freq';
    coh{loop}               = rmfield(coh{loop},'labelcmb');
    coh{loop}               = rmfield(coh{loop}, 'cohspctrm');
    coh_cmb{loop}           = ft_combineplanar([], coh{loop});

%     powCorrSpc{loop}.label         = {powCorrSpc{loop}.labelcmb{:,1}}';
%     powCorrSpc{loop}.powspctrm     = powCorrSpc{loop}.powcorrspctrm;
%     powCorrSpc{loop}.dimord        = 'chan_freq';
%     powCorrSpc{loop}               = rmfield(powCorrSpc{loop},'labelcmb');
%     powCorrSpc{loop}               = rmfield(powCorrSpc{loop}, 'powcorrspctrm');
%     powCorrSpc_cmb{loop}           = ft_combineplanar_itc([], powCorrSpc{loop});
% 
%     plv{loop}.label         = {plv{loop}.labelcmb{:,1}}';
%     plv{loop}.powspctrm     = plv{loop}.plvspctrm;
%     plv{loop}.dimord        = 'chan_freq';
%     plv{loop}               = rmfield(plv{loop},'labelcmb');
%     plv{loop}               = rmfield(plv{loop}, 'plvspctrm');
%     plv_cmb{loop}           = ft_combineplanar_itc([], plv{loop});

end

%%
figure
wat = [1:3];
cfg                 = [];
cfg.parameter       = 'powspctrm';
cfg.layout          = 'neuromag306cmb.lay';
cfg.graphcolor      = 'brkcmg'
% cfg.channel         = 'MEGMAG'
% cfg.xlim            = [1 99];
% cfg.ylim            = [0 0.2]
% ft_multiplotER(cfg, ft_combineplanar([],freq_fourier))
ft_multiplotER(cfg, PD_CMC_grndAvg{[1,3]}, Ctrl_CMC_grndAvg{[1,3]})

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

figure
cfg                         = [];
cfg.zlim                    = 'maxmin';
cfg.xlim                    = [0 1]; 
cfg.ylim                    = [1 100];
cfg.maskstyle               = 'saturation';	
cfg.masknans                = 'yes';
cfg.layout                  = 'neuromag306cmb.lay';
cfg.shading                 = 'interp';
cfg.colorbar                = 'yes';
cfg.axes                    = 'no';
% cfg.channel                 = {'MEG'};
ft_multiplotTFR(cfg, ft_combineplanar([],TFRwave));
    
%% EMG power

cfg                 = [];
cfg.output          = 'pow';
cfg.method          = 'mtmfft';
cfg.taper           = 'dpss';
cfg.foilim          = [1 200];
cfg.tapsmofrq       = 1.2*cfg.foilim;
% cfg.keeptrials      = 'yes';
cfg.channel         = {'EMG'};
freq          = ft_freqanalysis(cfg,preproc_data_EMG);

plot(freq.freq, freq.powspctrm)

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

%%

cfg = [];
cfg.parameter = 'cohspctrm';
figure;ft_connectivityplot(cfg, coh_temp{1});


%% Get max freq/CMC
L_par_ch = {'MEG0232+0233', 'MEG0242+0243', 'MEG0422+0423', 'MEG0432+0433',...
    'MEG0442+0443', 'MEG1612+1613', 'MEG1622+1623', 'MEG1632+1633', 'MEG1812+1813',...
    'MEG1822+1823', 'MEG1842+1843'};


R_par_ch = {'MEG1112+1113', 'MEG1132+1133', 'MEG1142+1143', 'MEG1332+1333', ...
    'MEG1342+1343', 'MEG2212+2213', 'MEG2222+2223', 'MEG2232+2233', 'MEG2412+2413',...
    'MEG2422+2423', 'MEG2442+2443'};

L_par_2= {'MEG0222+0223', 'MEG0232+0233', 'MEG0332+0333', 'MEG0412+0413', ...
    'MEG0422+0423', 'MEG0432+0433', 'MEG0442+0443', 'MEG0622+0623', 'MEG0632+0633',...
    'MEG0642+0643', 'MEG0712+0713', 'MEG0722+0723', 'MEG0732+0733', 'MEG0742+0743',...
    'MEG1042+1043', 'MEG1622+1623', 'MEG1812+1813', 'MEG1822+1823'};
new = {'MEG0412+0413', 'MEG0422+0423', 'MEG0432+0433', 'MEG0442+0443', 'MEG0712+0713'}

cfg_sel      = [];
cfg_sel.frequency = 1;

chan_pos_L = match_str(coh_cmb{1}.label, new);
chan_pos_R = match_str(coh_cmb{1}.label, R_par_ch);

for loop = 1:3
    
    coh_sel{loop} = ft_selectdata(cfg_sel, PD_CMC_grndAvg{loop});
    
    coh_sel{loop}.powspctrm = max(PD_CMC_grndAvg{loop}.powspctrm(:,13:30), [], 2);
    
    [maxCMC, maxFreq]     = max(PD_CMC_grndAvg{loop}.powspctrm(chan_pos_L,13:30), [], 2);
    
    maxFreqmax(loop)      = max(maxFreq+12);
    maxFreqstd(loop)      = std(maxFreq+12);
    
    
    maxCMCmax_PD(loop)      = mean(maxCMC);
    maxCMCstd(loop)      = std(maxCMC);
    
end


% maxCMC_res = reshape(maxCMCmax',3,[])'
maxCMC_PD_res = reshape(atanh(maxCMCmax_PD)',3,[])'

maxFreq_res = reshape(maxFreqmax',3,[])'

% reshape(maxFreqmean',3,[])'

% res = reshape(maxCMCmean',3,[])'

%%
figure,
bar([maxCMC_PD_res(:,[1 2 3]); maxCMC_C_res(:,[1 2 3])])



%% EMG spectra

for loop = 1:18
    
    subplot(6,3, loop)
    semilogy(freq_desc{1,loop}.freq, mean(freq_desc{1,loop}.powspctrm(chan_pos_L,:)))
    title(freq_desc{1,loop}.filename)
    
end
 
% export_fig( gcf, 'MEG spectra','-transparent', ...
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

%% grand avg
for loop = 1:18
    coh_cmb{loop}.idx = loop;
end

coh_cmb_res = reshape(coh_cmb,3,6)'
cfg = []; 
cfg.parameter = 'powspctrm'

for loop = 1:3
   
    PD_CMC_grndAvg{loop} = ft_freqgrandaverage(cfg, coh_cmb_res{1:3,loop});
    Ctrl_CMC_grndAvg{loop} = ft_freqgrandaverage(cfg, coh_cmb_res{4:6,loop});
    
end












