% test codes - CMC

%% Setup folders and files

% restoredefaultpath
addpath('C:\Program Files\MATLAB\R2012b\toolbox_add_on\fieldtrip-20150222')
addpath('C:\Program Files\MATLAB\R2012b\toolbox_add_on\export_fig')
addpath('C:\Program Files\MATLAB\R2012b\toolbox_add_on\CMC')
addpath('C:\Program Files\MATLAB\R2012b\toolbox_add_on\')
ft_defaults
 
proc.data_folder    = 'J:\MEG_Research\CMC\PD_move_artRej';
proc.results_folder = 'J:\MEG_Research\CMC\';

mkdir(proc.results_folder)
cd(proc.data_folder)
% filenames   = dir('*raw.fif');
filenames   = dir('*raw.fif');

for loop = 1: length(filenames)
    file_sub(loop) = {filenames(loop).name};
end

file_sub    = reshape(file_sub, 7, [])';
[row col]   = size(file_sub);


%%
% cd(proc.data_folder)
 subjLoop = 5;   condLoop = 2;
% for subjLoop=1:10 
%     for condLoop=1:7

        cfg                         = [];
        cfg.dataset                 = file_sub{subjLoop,condLoop};
        cfg.channel                 = {'EMG'};%, 'EOG'};
        data_import                 = ft_preprocessing(cfg);

        % Filter EMG
        cfg                         = [];
%         cfg.lpfilter                = 'yes';
%         cfg.lpfreq                  =  20;
%         cfg.hpfilter                = 'yes';
%         cfg.hpfreq                  =  2;
        cfg.bpfilter                = 'yes';
        cfg.bpfreq                  =  [10 100];
        cfg.dftfilter               = 'yes';
%         cfg.channel                 = {'EMG'}; 
%         cfg.rectify                 = 'yes'
%         cfg.hilbert                 = 'abs';
        cfg.detrend                 = 'yes';
        preproc_data                = ft_preprocessing(cfg, data_import);
        
        hFig = figure(1),
        len  = length(preproc_data.label);
        for nChan = 1:len
            subplot(len, 1, nChan)
            plot(preproc_data.time{1},preproc_data.trial{1}(nChan,:)), 
            title(preproc_data.label(nChan))
        end
        suptitle(file_sub{subjLoop,condLoop}(1:end-4))
        waitfor(hFig)

%     end
% end
% 
%% 
        hFig = figure(subjLoop+10), 
%         subplot(3,3,condLoop)
        figure
        plot(preproc_data_EMG.time{1},preproc_data_EMG.trial{1})

        thresh_EMG=preproc_data_EMG.trial{1};
        
        %%
        cfg          = [];
        cfg.method   = 'channel';
        cfg.ylim     = [0 0.01]; 
        cfg.blocksize = 20;
%         cfg.megscale = 1;
%         cfg.eogscale = 5e-8;
        dummy        = ft_databrowser(cfg,preproc_data);

        %%        

        cfg                                 = [];
        cfg.artfctdef.visual.artifact       = dummy.artfctdef.visual.artifact;%emg_trl;
        cfg.artfctdef.reject                = 'nan';
        EMGclean                          = ft_rejectartifact(cfg,preproc_data_EMG);
        EMGclean.trial{1}(1,isnan(EMGclean.trial{1}))=0;
        
        
        cfg          = [];
        cfg.method   = 'channel';
        cfg.ylim     = [0 0.01]; 
        cfg.blocksize = 20;
        ft_databrowser(cfg, EMGepoc)
        
        
        %%
        figure
        plot(preproc_data_EMG.time{1},preproc_data_EMG.trial{1})
        thresh_EMG = preproc_data_EMG.trial{1};
        thresh_EMG(1,find(thresh_EMG>15e-4))=0;
        hold on
        plot(preproc_data_EMG.time{1},thresh_EMG, 'r')
%     end
%     suptitle(num2str(subjLoop))
% end

%%

EMGpos = find(strcmp(data_import.label,ft_channelselection('EMG*', data_import.label)))
data_import.trial{1}(EMGpos,:) = preproc_data_EMG.trial{1};

cfg                                     = [];   
% cfg.datafile                            = file_sub{subjLoop,condLoop}
cfg.trialfun                            = 'trialfun_emgdetect_visual';
cfg.thresh                              = 1;
trl_list                                = ft_definetrial(cfg);


cfg         = [];
cfg.thresh  = 1;
cfg.shift   = 0.5; % 500ms
cfg.overlap = 0.1; % 10% overlap
cfg.savepdf = 'no';
emg_trl     = trialfun_emgdetect_visual(cfg,EMGclean);


cfg                                     = [];   
% cfg.datafile                            = file_sub{subjLoop,condLoop}
cfg.trialfun                            = 'trialfun_emgdetect_visual';
cfg.thresh                              = 1;
cfg.trials = 'all';

cfg.trl = emg_trl
data_epoched                            = ft_redefinetrial(cfg,data_import);

%%
cfg                     = [];
cfg.output              = 'pow';
cfg.method      = 'mtmfft';
cfg.channel = 'EMG*';
cfg.taper       = 'dpss';
cfg.foilim      = [1 45];
cfg.tapsmofrq   = 5;
cfg.keeptrials          = 'yes';
freqEMG = ft_freqanalysis(cfg, data_epoched);

figure,
semilogy(freqEMG.freq, squeeze(mean(freqEMG.powspctrm)))

%% EMG power spectra

cfg                 = [];
cfg.parameter       = 'powspctrm';
cfg.layout          = 'neuromag306cmb.lay';
% cfg.graphcolor      = [clrOp1];%clrOp2;clrOp3;clrOp4;clrOp5;clrOp6;clrOp7];
cfg.channel   = 'EMG*';
cfg.linewidth = 3;
cfg.colorbar = 'yes';
cfg.xlim = [10 45];
cfg.ylim = [0 3e-10 ];
figure
% subplot 211
% ft_singleplotER(cfg, PD_pow_grndAvg{:})%, coh_four_cmb)
% subplot 212

for condLoop=1:7
    figure(condLoop)
    subplot 211
    ft_singleplotER(cfg, PDfreq_cmb{[6:10],condLoop})%, coh_four_cmb)
    title('handgrips')

    subplot 212
    ft_singleplotER(cfg, PDfreq_cmb{[1:4],condLoop})%, coh_four_cmb)
    title('finger move')

end

%%
cfg             = [];
cfg.method      = 'mtmfft';
cfg.taper       = 'dpss';
cfg.output      = 'powandcsd';
cfg.keeptrials  = 'yes';
cfg.jackknife   = 'yes';
cfg.channelcmb  = {'MEG', 'EMG'};
cfg.foilim      = [1 45];
cfg.tapsmofrq   = 5;
freq_csd    = ft_freqanalysis(cfg, data_epoched);

cfg            = [];
cfg.output     = 'fourier';
cfg.keeptrials  = 'yes';
cfg.method     = 'mtmfft';
cfg.foilim     = [1 45];
cfg.tapsmofrq  = 5;
cfg.keeptrials = 'yes';
cfg.channel    = {'MEG' 'EMG'};
freqfourier    = ft_freqanalysis(cfg, data_epoched);

% CMC calc
cfg                             = [];
cfg.method                      = 'granger';

cfg.channelcmb                  = {'MEG*3','MEG*2', 'EMG'};
% cfg.complex                     = 'imag';
% cfg.jackknife                   = 'yes';
coh_csd                         = ft_connectivityanalysis(cfg, freq_csd);
% coh_four                         = ft_connectivityanalysis(cfg, freqfourier);

% coh_csd.label         = {coh_csd.labelcmb{:,1}}';
% coh_csd.powspctrm     = coh_csd.cohspctrm;
% coh_csd.dimord        = 'chan_freq';
% coh_csd               = rmfield(coh_csd,'labelcmb');
% coh_csd               = rmfield(coh_csd, 'cohspctrm');


coh_four.label         = {coh_four.labelcmb{:,1}}';
coh_four.powspctrm     = coh_four.cohspctrm;
coh_four.dimord        = 'chan_freq';
coh_four               = rmfield(coh_four,'labelcmb');
coh_four               = rmfield(coh_four, 'cohspctrm');


cfg_cmb = []; cfg_cmb.combinemethod = 'cmc';
% coh_csd_cmb          = ft_combineplanar_cmc(cfg_cmb, coh_four);
coh_four_cmb          = ft_combineplanar_cmc(cfg_cmb, coh_four);

% plots

cfg                 = [];
cfg.parameter       = 'powspctrm';
cfg.layout          = 'neuromag306cmb.lay';
cfg.graphcolor      = 'brkcmgy';
% cfg.channel         = R_par_ch;
% cfg.highlight = 'on'
% cfg.highlightchannel = [L_par_ch;R_par_ch];
% cfg.highlightsymbol    = 'x';
% cfg.highlightsize     = 12;
% cfg.xlim            = [13 30];
% cfg.zlim            = [2e-24 5e-24]
cfg.colorbar = 'yes';
figure
ft_multiplotER(cfg, coh_four_cmb)%, coh_four_cmb)

figure
cfg.xlim            = [13 30];
ft_topoplotER(cfg, coh_four_cmb)%, coh_four_cmb)




