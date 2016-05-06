clc;    clear all;    close all;

%% Setup folders and files

% restoredefaultpath
addpath('C:\Program Files\MATLAB\R2012b\toolbox_add_on\fieldtrip-20150222')
addpath('C:\Program Files\MATLAB\R2012b\toolbox_add_on\export_fig')
addpath('C:\Program Files\MATLAB\R2012b\toolbox_add_on\CMC')
addpath('C:\Program Files\MATLAB\R2012b\toolbox_add_on\')
ft_defaults
 
proc.data_folder    = 'J:\MEG_Research\CMC\raw_ArtRej\Control';
proc.results_folder = 'J:\MEG_Research\CMC\raw_ArtRej\results';

mkdir(proc.results_folder)
cd(proc.data_folder)
filenames   = dir('*raw.fif');

for loop = 1: length(filenames)
    file_sub(loop) = {filenames(loop).name};
end

file_sub    = reshape(file_sub, 7, [])';
[row col]   = size(file_sub);

%% Loop over datasets

for subjLoop=1:row
    for condLoop = 1:col
        
%         % Load data
        cfg                         = [];
        cfg.dataset                 = file_sub{subjLoop,condLoop};
        cfg.channel                 = {'EMG'};%{'MEG', 'EMG'};
        data_import                 = ft_preprocessing(cfg);

        % Filter EMG
        cfg                         = [];
%         cfg.lpfilter                = 'yes';
%         cfg.lpfreq                  =  100;
        cfg.hpfilter                = 'yes';
        cfg.hpfreq                  =  10;
        cfg.dftfilter               = 'yes';
        cfg.channel                 =  {'EMG'}; 
        cfg.rectify                 = 'yes'
        cfg.detrend                 = 'yes';
        preproc_data_EMG            = ft_preprocessing(cfg, data_import);
    %     plot(preproc_data_EMG.time{1}, preproc_data_EMG.trial{1})

        EMGpos = find(strcmp(data_import.label,ft_channelselection('EMG*', data_import.label)))
        data_import.trial{1}(EMGpos,:) = preproc_data_EMG.trial{1};

      % Epoch
        cfg                                     = [];   
        cfg.datafile                            = file_sub{subjLoop,condLoop}
        cfg.trialfun                            = 'trialfun_emgdetect';
        cfg.thresh                              = 1;
        cfg                                     = ft_definetrial(cfg);
        
        data_epoched                            = ft_redefinetrial(cfg,data_import);
        data_trl{subjLoop, condLoop}.trl        = data_epoched.cfg.trl;
        data_trl{subjLoop, condLoop}.filename   = file_sub{subjLoop,condLoop};

%         % Power calculation
%         cfg                     = [];
%         cfg.output              = 'pow';
%         cfg.method              = 'mtmfft';
%         cfg.taper               = 'hanning';
%         cfg.foi                 = [1:1:45];
%         cfg.keeptrials          = 'yes';
%         freq{subjLoop,condLoop}              = ft_freqanalysis(cfg, data_epoched);
%         freq{subjLoop,condLoop}.filename     = file_sub{subjLoop,condLoop};
% 
%         cfg                      = [];
%         cfg.jackknife            = 'yes';
%         freq_desc{subjLoop,condLoop}          = ft_freqdescriptives(cfg,  freq{subjLoop,condLoop});
%         freq_desc{subjLoop,condLoop}.filename = file_sub{subjLoop,condLoop};
%         
%         freq_cmb{subjLoop,condLoop}  = ft_combineplanar([],freq_desc{subjLoop,condLoop});
%         freq_cmb{subjLoop,condLoop}.filename = file_sub{subjLoop,condLoop};
% 
%         % Stability
%         cfg                 = [];
%         cfg.medianfilter    = 'yes';
%         cfg.medianfiltord   = 50;
%         cfg.channel         = 'EMG';
%         EMG_ch              = ft_preprocessing(cfg, data_epoched);
%         
%         for trialLoop=1:length(EMG_ch.trial)
%             stab(trialLoop) = 1-(std(EMG_ch.trial{trialLoop})/mean(EMG_ch.trial{trialLoop}));
%         end
%         
%         stability{subjLoop,condLoop}.stab           = stab;
%         stability{subjLoop,condLoop}.filename       =  file_sub{subjLoop,condLoop};
%         stability_mean(subjLoop,condLoop)           = mean(stab);
%         stability_std(subjLoop,condLoop)            = std(stab);
          
        cfg             = [];
        cfg.method      = 'mtmfft';
        cfg.taper       = 'dpss';
        cfg.output      = 'powandcsd';
        cfg.keeptrials  = 'yes';
        cfg.jackknife   = 'yes';
        cfg.channelcmb  = {'MEG', 'EMG'};
        cfg.foi         = [1:45];
        cfg.tapsmofrq   = 2;
        freq_fourier    = ft_freqanalysis(cfg, data_epoched);

        % CMC calc
        cfg                             = [];
        cfg.method                      = 'coh';
        cfg.channelcmb                  = {'MEG', 'EMG'};
        cfg.jackknife                   = 'yes';
        coh{subjLoop,condLoop}          = ft_connectivityanalysis(cfg, freq_fourier);
        coh{subjLoop,condLoop}.filename =  file_sub{subjLoop,condLoop};
% 
% 
%         coh{subjLoop,condLoop}.label         = {coh{subjLoop,condLoop}.labelcmb{:,1}}';
%         coh{subjLoop,condLoop}.powspctrm     = coh{subjLoop,condLoop}.cohspctrm;
%         coh{subjLoop,condLoop}.dimord        = 'chan_freq';
%         coh{subjLoop,condLoop}               = rmfield(coh{subjLoop,condLoop},'labelcmb');
%         coh{subjLoop,condLoop}               = rmfield(coh{subjLoop,condLoop}, 'cohspctrm');
        cfg_cmb = []; cfg_cmb.combinemethod = 'cmc';
        Ctrlcoh_cmb{subjLoop,condLoop}           = ft_combineplanar_cmc(cfg_cmb, Ctrlcoh{subjLoop,condLoop});
        Ctrlcoh_cmb{subjLoop,condLoop}.filename  =  Ctrlcoh{subjLoop,condLoop}.filename;

        cfg_cmb = []; cfg_cmb.combinemethod = 'cmc';
        PDcoh_cmb{subjLoop,condLoop}           = ft_combineplanar_cmc(cfg_cmb, PDcoh{subjLoop,condLoop});
        PDcoh_cmb{subjLoop,condLoop}.filename  = PDcoh{subjLoop,condLoop}.filename;

        
%         data_scramble = data_epoched;
%         for trialLoop = 1:length(data_scramble.trial)-3
% 
%             data_scramble.trial{trialLoop}(1,:) = ...
%                 data_scramble.trial{trialLoop+3}(1,:);
% 
%         end
%         
%         cfg_sel.trials = 1:trialLoop;
%         data_scramble = ft_selectdata(cfg_sel, data_scramble);
% 
%         cfg             = [];
%         cfg.method      = 'mtmfft';
%         cfg.taper       = 'dpss';
%         cfg.output      = 'powandcsd';
%         cfg.keeptrials  = 'yes';
%         cfg.jackknife   = 'yes';
%         cfg.channelcmb  = {'MEG', 'EMG'};
%         cfg.foi         = [1:45];
%         cfg.tapsmofrq   = 2;
%         freq_fourier_scr    = ft_freqanalysis(cfg, data_scramble);
% 
%         % CMC calc
%         cfg                = [];
%         cfg.method         = 'coh';
%         cfg.channelcmb     = {'MEG', 'EMG'};
%         cfg.jackknife      = 'yes';
%         coh_scr{subjLoop,condLoop}          = ft_connectivityanalysis(cfg, freq_fourier_scr);
%         coh_scr{subjLoop,condLoop}.filename =  file_sub{subjLoop,condLoop};
% 
% 
%         coh_scr{subjLoop,condLoop}.label         = {coh_scr{subjLoop,condLoop}.labelcmb{:,1}}';
%         coh_scr{subjLoop,condLoop}.powspctrm     = coh_scr{subjLoop,condLoop}.cohspctrm;
%         coh_scr{subjLoop,condLoop}.dimord        = 'chan_freq';
%         coh_scr{subjLoop,condLoop}               = rmfield(coh_scr{subjLoop,condLoop},'labelcmb');
%         coh_scr{subjLoop,condLoop}               = rmfield(coh_scr{subjLoop,condLoop}, 'cohspctrm');
%         coh_scr_cmb{subjLoop,condLoop}           = ft_combineplanar([], coh_scr{subjLoop,condLoop});
        
        
    end
end

save Ctrlcoh_cmb Ctrlcoh_cmb -v7.3
save PDcoh_cmb PDcoh_cmb -v7.3

%%

save Ctrlcoh Ctrlcoh -v7.3
save Ctrlcoh_cmb Ctrlcoh_cmb -v7.3
save Ctrlcoh_scr Ctrlcoh_scr -v7.3
save Ctrlcoh_scr_cmb Ctrlcoh_scr_cmb -v7.3
save Ctrlfreq Ctrlfreq -v7.3
save Ctrlfreq_desc Ctrlfreq_desc -v7.3
save Ctrlfreq_cmb Ctrlfreq_cmb -v7.3
save Ctrlstability Ctrlstability -v7.3
save Ctrlstability_mean Ctrlstability_mean
save Ctrldata_trl Ctrldata_trl -v7.3
save Ctrlfile_sub Ctrlfile_sub







