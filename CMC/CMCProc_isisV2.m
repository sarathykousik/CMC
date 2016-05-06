clc;    clear all;    close all;

%% Setup folders and files

% restoredefaultpath
% addpath('/usr/local/common/matlab_toolbox/fieldtrip/r6152')
%addpath('C:\Program Files\MATLAB\R2012b\toolbox_add_on\export_fig')
addpath('C:\Program Files\MATLAB\R2012b\toolbox_add_on\cmc\')
%addpath('C:\Program Files\MATLAB\R2012b\toolbox_add_on\')
ft_defaults
 
proc.data_folder_control    = 'J:\MEG_Research\CMC\redo';
proc.data_folder_PD    = 'J:\MEG_Research\CMC\redo';
proc.results_folder = 'J:\MEG_Research\CMC\redo\results\';

mkdir(proc.results_folder)
cd(proc.data_folder_control)
filenames   = dir('*-raw.fif');

for loop = 1: length(filenames)
    file_sub(loop) = {filenames(loop).name};
end

file_sub    = reshape(file_sub, 3, [])';
[row col]   = size(file_sub);

%% Loop over datasets

for subjLoop = 1:row
    for condLoop = 1:col
        
%         % Load data
        cfg                         = [];
        cfg.dataset                 = file_sub{subjLoop,condLoop};
        cfg.channel                 = {'MEG'};
        data_import_MEG                 = ft_preprocessing(cfg);
        
        cfg.channel                 = {'EMG'};
        data_import_EMG             = ft_preprocessing(cfg);

%        EMGpos = find(strcmp(data_import.label,ft_channelselection('EMG*', data_import.label)))

%         if (subjLoop==5 && condLoop~=1) || (subjLoop==1 && condLoop==7)
%             disp('######  Using custom trial defn')
%             load(['j:\MEG_Research\CMC\017_trial\',file_sub{subjLoop, condLoop}(1:end-4),...
%                 '-trlVisRej.mat'])
%             cfg                                 = [];
%             cfg.artfctdef.visual.artifact       = visRej.artfctdef.visual.artifact;%emg_trl;
%             cfg.artfctdef.reject                = 'nan';
%             EMGclean                          = ft_rejectartifact(cfg,data_import_EMG);
%             EMGclean.trial{1}(1,isnan(EMGclean.trial{1}))=0;
%             data_import_EMG.trial{1} = EMGclean.trial{1};
%         end
        
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
        preproc_data_EMG            = ft_preprocessing(cfg, data_import_EMG);
        
        cfg                         = [];
        cfg.lpfilter                = 'yes';
        cfg.lpfreq                  =  100;
        cfg.hpfilter                = 'yes';
        cfg.hpfreq                  =  1;
        cfg.dftfilter               = 'yes';
        cfg.channel                 =  {'MEG'}; 
        cfg.detrend                 = 'yes';
        preproc_data_MEG            = ft_preprocessing(cfg, data_import_MEG);
        
        preproc_data = ft_appenddata([],preproc_data_MEG,preproc_data_EMG);
        load('J:\MEG_Research\CMC\redo\data_trl.mat')
        % Epoch
        cfg                                     = [];   
        cfg.datafile                            = file_sub{subjLoop,condLoop};
%         cfg.trialfun                            = 'trialfun_emgdetect_visFB';
%         cfg.save_folder                         = proc.results_folder;
%         cfg.thresh                              = 1;
%         cfg.emgSignal                           = data_import_EMG.trial{1};
%         cfg.shift                               = 0.5; % in secs
%         cfg.overlap                             = 0.05;  % in % - here 5%
%         cfg                                     = ft_definetrial(cfg);
        cfg.trl                                 = data_trl{subjLoop,condLoop}.trl;
%         
        data_epoched                            = ft_redefinetrial(cfg,preproc_data);
%         data_trl{subjLoop, condLoop}.trl        = data_epoched.cfg.trl;
%         data_trl{subjLoop, condLoop}.filename   = file_sub{subjLoop,condLoop};

        % Power calculation
%         cfg                     = [];
%         cfg.output              = 'pow';
%         cfg.method              = 'mtmfft';
%         cfg.taper               = 'hanning';
        
        %alpha mtm
%         cfg.foi              = 9;
%         cfg.tapsmofrq        = 3;
%         cfg.keeptrials          = 'yes';
%         freq_alpha{subjLoop,condLoop}              = ft_freqanalysis(cfg, data_epoched);
%         freq_alpha{subjLoop,condLoop}.filename     = file_sub{subjLoop,condLoop};
%         
        %beta mtm
%         cfg.foi              = 22;
%         cfg.tapsmofrq        = 9;
%         cfg.keeptrials          = 'yes';
%         freq_beta{subjLoop,condLoop}              = ft_freqanalysis(cfg, data_epoched);
%         freq_beta{subjLoop,condLoop}.filename     = file_sub{subjLoop,condLoop};
        
        %gamma mtm
%         cfg.foi             = 38;
%         cfg.tapsmofrq       = 7;
%         cfg.keeptrials          = 'yes';
%         freq_gamma{subjLoop,condLoop}              = ft_freqanalysis(cfg, data_epoched);
%         freq_gamma{subjLoop,condLoop}.filename     = file_sub{subjLoop,condLoop};
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
        cfg                 = [];
        cfg.medianfilter    = 'yes';
        cfg.medianfiltord   = 50;
        cfg.channel         = 'EMG';
        EMG_ch              = ft_preprocessing(cfg, data_epoched);
        
        for trialLoop=1:length(EMG_ch.trial)
            stab(trialLoop) = 1-(std(EMG_ch.trial{trialLoop})/mean(EMG_ch.trial{trialLoop}));
        end
        
        stability{subjLoop,condLoop}.stab           = stab;
        stability{subjLoop,condLoop}.filename       =  file_sub{subjLoop,condLoop};
        stability_mean(subjLoop,condLoop)           = mean(stab);
        stability_std(subjLoop,condLoop)            = std(stab);

        cfg             = [];
        cfg.method      = 'mtmfft';
        cfg.taper       = 'dpss';
        cfg.output      = 'powandcsd';
        cfg.keeptrials  = 'yes';
        cfg.jackknife   = 'yes';
        cfg.channelcmb  = {'MEG', 'EMG'};
        cfg.foi         = [1:100];
        cfg.tapsmofrq   = 5;
        freq_fourier{subjLoop, condLoop}    = ft_freqanalysis(cfg, data_epoched);
        freq_fourier{subjLoop, condLoop}.filename   = file_sub{subjLoop,condLoop};
% % 
% %         cfg.foi              = 9;
% %         cfg.tapsmofrq        = 3;
% %         freq_fourier_alpha{subjLoop, condLoop}   = ft_freqanalysis(cfg, data_epoched);
% %         freq_fourier_alpha{subjLoop, condLoop}.filename   = file_sub{subjLoop,condLoop};
% 
        cfg.foi             = 22;
        cfg.tapsmofrq       = 9;
        freq_fourier_beta{subjLoop, condLoop}   = ft_freqanalysis(cfg, data_epoched);
        freq_fourier_beta{subjLoop, condLoop}.filename   = file_sub{subjLoop,condLoop};

%         cfg.foi             = 38;
%         cfg.tapsmofrq       = 7;
%         freq_fourier_gamma{subjLoop, condLoop}   = ft_freqanalysis(cfg, data_epoched);
%         freq_fourier_gamma{subjLoop, condLoop}.filename   = file_sub{subjLoop,condLoop};

         
%         % CMC calc
        cfg                             = [];
        cfg.method                      = 'coh';
        cfg.channelcmb                  = {'MEG', 'EMG'};
        cfg.jackknife                   = 'yes';
        coh{subjLoop,condLoop}          = ft_connectivityanalysis(cfg, freq_fourier{subjLoop, condLoop});
        coh{subjLoop,condLoop}.filename =  file_sub{subjLoop,condLoop};
%         
%         coh_alpha{subjLoop,condLoop}          = ft_connectivityanalysis(cfg, freq_fourier_alpha);
%         coh_alpha{subjLoop,condLoop}.filename =  file_sub{subjLoop,condLoop};
        coh_beta{subjLoop,condLoop}           = ft_connectivityanalysis(cfg,    freq_fourier_beta{subjLoop, condLoop});
        coh_beta{subjLoop,condLoop}.filename  =  file_sub{subjLoop,condLoop};
%         coh_gamma{subjLoop,condLoop}          = ft_connectivityanalysis(cfg, freq_fourier_gamma);
%         coh_gamma{subjLoop,condLoop}.filename =  file_sub{subjLoop,condLoop};
% 
        coh{subjLoop,condLoop}.label         = {coh{subjLoop,condLoop}.labelcmb{:,1}}';
        coh{subjLoop,condLoop}.powspctrm     = coh{subjLoop,condLoop}.cohspctrm;
        coh{subjLoop,condLoop}.dimord        = 'chan_freq';
        coh{subjLoop,condLoop}               = rmfield(coh{subjLoop,condLoop},'labelcmb');
        coh{subjLoop,condLoop}               = rmfield(coh{subjLoop,condLoop}, 'cohspctrm');
        coh{subjLoop,condLoop}.filename      = file_sub{subjLoop,condLoop};
%         
%         coh_alpha{subjLoop,condLoop}.label         = {coh_alpha{subjLoop,condLoop}.labelcmb{:,1}}';
%         coh_alpha{subjLoop,condLoop}.powspctrm     = coh_alpha{subjLoop,condLoop}.cohspctrm;
%         coh_alpha{subjLoop,condLoop}.dimord        = 'chan_freq';
%         coh_alpha{subjLoop,condLoop}               = rmfield(coh_alpha{subjLoop,condLoop},'labelcmb');
%         coh_alpha{subjLoop,condLoop}               = rmfield(coh_alpha{subjLoop,condLoop}, 'cohspctrm');
%         
        coh_beta{subjLoop,condLoop}.label         = {coh_beta{subjLoop,condLoop}.labelcmb{:,1}}';
        coh_beta{subjLoop,condLoop}.powspctrm     = coh_beta{subjLoop,condLoop}.cohspctrm;
        coh_beta{subjLoop,condLoop}.dimord        = 'chan_freq';
        coh_beta{subjLoop,condLoop}               = rmfield(coh_beta{subjLoop,condLoop},'labelcmb');
        coh_beta{subjLoop,condLoop}               = rmfield(coh_beta{subjLoop,condLoop}, 'cohspctrm');
%         
%         coh_gamma{subjLoop,condLoop}.label         = {coh_gamma{subjLoop,condLoop}.labelcmb{:,1}}';
%         coh_gamma{subjLoop,condLoop}.powspctrm     = coh_gamma{subjLoop,condLoop}.cohspctrm;
%         coh_gamma{subjLoop,condLoop}.dimord        = 'chan_freq';
%         coh_gamma{subjLoop,condLoop}               = rmfield(coh_gamma{subjLoop,condLoop},'labelcmb');
%         coh_gamma{subjLoop,condLoop}               = rmfield(coh_gamma{subjLoop,condLoop}, 'cohspctrm');
%         
        cfg_cmb = []; cfg_cmb.combinemethod = 'cmc';
        coh_cmb{subjLoop,condLoop}           = ft_combineplanar_cmc(cfg_cmb, coh{subjLoop,condLoop});
        coh_cmb{subjLoop,condLoop}.filename  =  coh{subjLoop,condLoop}.filename;
%        
%         coh_alpha_cmb{subjLoop,condLoop}           = ft_combineplanar_cmc(cfg_cmb, coh_alpha{subjLoop,condLoop});
%         coh_alpha_cmb{subjLoop,condLoop}.filename  =  coh_alpha{subjLoop,condLoop}.filename;
% 
        coh_beta_cmb{subjLoop,condLoop}           = ft_combineplanar_cmc(cfg_cmb, coh_beta{subjLoop,condLoop});
        coh_beta_cmb{subjLoop,condLoop}.filename  =  coh{subjLoop,condLoop}.filename;
% 
%         
%         coh_gamma_cmb{subjLoop,condLoop}           = ft_combineplanar_cmc(cfg_cmb, coh_gamma{subjLoop,condLoop});
%         coh_gamma_cmb{subjLoop,condLoop}.filename  =  coh{subjLoop,condLoop}.filename;
% 
%       
        
        data_scramble = data_epoched;
        for trialLoop = 1:length(data_scramble.trial)-3

            data_scramble.trial{trialLoop}(1,:) = ...
                data_scramble.trial{trialLoop+3}(1,:);

        end
        
        cfg_sel.trials = 1:trialLoop;
        data_scramble = ft_selectdata(cfg_sel, data_scramble);

        cfg                 = [];
        cfg.method          = 'mtmfft';
        cfg.taper           = 'dpss';
        cfg.output          = 'powandcsd';
        cfg.keeptrials      = 'yes';
        cfg.jackknife       = 'yes';
        cfg.channelcmb      = {'MEG', 'EMG'};
        cfg.foi             = [1:100];
        cfg.tapsmofrq       = 5;
        freq_fourier_scr    = ft_freqanalysis(cfg, data_scramble);

        % CMC calc
        cfg                                 = [];
        cfg.method                          = 'coh';
        cfg.channelcmb                      = {'MEG', 'EMG'};
        cfg.jackknife                       = 'yes';
        coh_scr{subjLoop,condLoop}          = ft_connectivityanalysis(cfg, freq_fourier_scr);
        coh_scr{subjLoop,condLoop}.filename =  file_sub{subjLoop,condLoop};


        coh_scr{subjLoop,condLoop}.label         = {coh_scr{subjLoop,condLoop}.labelcmb{:,1}}';
        coh_scr{subjLoop,condLoop}.powspctrm     = coh_scr{subjLoop,condLoop}.cohspctrm;
        coh_scr{subjLoop,condLoop}.dimord        = 'chan_freq';
        coh_scr{subjLoop,condLoop}               = rmfield(coh_scr{subjLoop,condLoop},'labelcmb');
        coh_scr{subjLoop,condLoop}               = rmfield(coh_scr{subjLoop,condLoop}, 'cohspctrm');
        coh_scr_cmb{subjLoop,condLoop}           = ft_combineplanar([], coh_scr{subjLoop,condLoop});
        
        
    end
end

%%
save coh coh -v7.3
save coh_cmb coh_cmb -v7.3
save coh_scr coh_scr -v7.3
save coh_scr_cmb coh_scr_cmb -v7.3
save data_trl data_trl -v7.3




%%

save Ctrlcoh Ctrlcoh -v7.3
save Ctrlcoh_cmb Ctrlcoh_cmb -v7.3
save Ctrlcoh_scr Ctrlcoh_scr -v7.3
save Ctrlcoh_scr_cmb Ctrlcoh_scr_cmb -v7.3
save Ctrlfreq Ctrlfreq -v7.3
save Ctrlfreq_desc Ctrlfreq_desc -v7.3
save Ctrlfreq_cmb Ctrlfreq_cmb -v7.3
save Ctrlstability Ctrlstability -v7.3
save Ctrlstability_mean Ctrlstability_mean -v7.3
save Ctrldata_trl Ctrldata_trl -v7.3
save Ctrlfile_sub Ctrlfile_sub -v7.3
save Ctrlcoh_alpha Ctrlcoh_alpha -v7.3
save Ctrlcoh_alpha_cmb Ctrlcoh_alpha_cmb -v7.3
save Ctrlcoh_beta Ctrlcoh_beta -v7.3
save Ctrlcoh_beta_cmb Ctrlcoh_beta_cmb -v7.3
save Ctrlcoh_gamma Ctrlcoh_gamma -v7.3
save Ctrlcoh_gamma_cmb Ctrlcoh_gamma_cmb -v7.3


%%

save PDcoh PDcoh -v7.3
save PDcoh_cmb PDcoh_cmb -v7.3
save PDcoh_scr PDcoh_scr -v7.3
save PDcoh_scr_cmb PDcoh_scr_cmb -v7.3
save PDfreq PDfreq -v7.3
save PDfreq_desc PDfreq_desc -v7.3
save PDfreq_cmb PDfreq_cmb -v7.3
save PDstability PDstability -v7.3
save PDstability_mean PDstability_mean -v7.3
save PDdata_trl PDdata_trl -v7.3
save PDfile_sub PDfile_sub -v7.3
save PDcoh_alpha PDcoh_alpha -v7.3
save PDcoh_alpha_cmb PDcoh_alpha_cmb -v7.3
save PDcoh_beta PDcoh_beta -v7.3
save PDcoh_beta_cmb PDcoh_beta_cmb -v7.3
save PDcoh_gamma PDcoh_gamma -v7.3
save PDcoh_gamma_cmb PDcoh_gamma_cmb -v7.3

%%
save PDfreq_alpha PDfreq_alpha -v7.3
save PDfreq_beta PDfreq_beta -v7.3
save PDfreq_gamma PDfreq_gamma -v7.3
save PDfreq_fourier PDfreq_fourier -v7.3
save PDfreq_fourier_alpha PDfreq_fourier_alpha -v7.3
save PDfreq_fourier_beta PDfreq_fourier_beta -v7.3
save PDfreq_fourier_gamma PDfreq_fourier_gamma -v7.3

%%

for subjLoop=1:10
    for condLoop=1:7
        
        
        cfg_cmb = []; cfg_cmb.combinemethod = 'cmc';
        PDcoh_cmb{subjLoop,condLoop}           = ft_combineplanar_cmc(cfg_cmb, PDcoh{subjLoop,condLoop});
        PDcoh_cmb{subjLoop,condLoop}.filename  =  PDcoh{subjLoop,condLoop}.filename;     
        
    end
end
save PDcoh_cmb PDcoh_cmb -v7.3


