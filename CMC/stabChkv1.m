clc; close all; clear all;

%% Setup folders and files

proc.data_folder    = 'J:\MEG_Research\CMC\Control_move_RAW';
% proc.results_folder = 'J:\MEG_Research\CMC\raw_ArtRej\results';

cd(proc.data_folder)
filenames      = dir('*.fif');

for loop = 1: length(filenames)
    file_sub(loop) = {filenames(loop).name};
end
% 
file_sub = reshape(file_sub, 7, [])';
[row col] = size(file_sub);

%%
for subjLoop = 1:row
    for condLoop = 1:col
        % Import only EMG channel, epoch it
        % Load data
        cfg                         = [];
        cfg.dataset                 = file_sub(subjLoop,condLoop);
        cfg.channel                 = {'EMG'};
        cfg.lpfilter                = 'yes';
        cfg.lpfreq                  =  100;
        cfg.hpfilter                = 'yes';
        cfg.hpfreq                  =  1;
        cfg.dftfilter               = 'yes';
        cfg.channel                 =  {'EMG'}; 
        cfg.rectify                 = 'yes'
        cfg.detrend                 = 'yes';

        cfg.medianfilter            = 'yes';
        cfg.medianfiltord           = 50;

        data_import                 = ft_preprocessing(cfg);

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
        data_epoched                = ft_redefinetrial(cfg, data_import);
%         data_epoched{subjLoop,condLoop}.filename = filenames(loop).name;


        % Preprocessing - Filter-notch-BP, rectify, median filter
        cfg                 = [];
        cfg.medianfilter    = 'yes';
        cfg.medianfiltord   = 50;
        cfg.channel         = 'EMG';
        EMG_ch              = ft_preprocessing(cfg, data_epoched);


        % Estimate stability for each subject

        % Stability

            for trialLoop=1:length(EMG_ch.trial)
                stab(trialLoop) = 1-(std(EMG_ch.trial{trialLoop})/mean(EMG_ch.trial{trialLoop}));
            end
        Control_stability{subjLoop,condLoop}.stab           = stab;
        Control_stability{subjLoop,condLoop}.filename       = file_sub{subjLoop,condLoop};
        Control_stability_mean(subjLoop,condLoop)           = mean(stab);
        Control_stability_std(subjLoop,condLoop)            = std(stab);
    
    end
end

run cmc.m
save Control_stability Control_stability -v7.3
save Control_stability_mean Control_stability_mean
save Control_stability_std Control_stability_std
    
figure, boxplot(PD_stability_mean([9,10,12],:)), title('PD Stab 9,10,12'), ylim([0 0.8])
figure, boxplot(Control_stability_mean([1,3,4],:)), title('C Stab 1,3,4'), ylim([0 0.8])

figure, boxplot(PD_stability_mean([1:5],:)), title('PD Stab 1-5'), ylim([0 0.8])

% export_fig( gcf, 'Control-stab-mean','-transparent', ...
%         '-painters','-pdf', '-r250' ); 

figure, boxplot(Control_stability_std), title('Control Stab std'), ylim([0 0.8])
% export_fig( gcf, 'Control-stab-std','-transparent', ...
%         '-painters','-pdf', '-r250' ); 

    
    
anova_rm({PD_stability_mean([9,10,12],:), Control_stability_mean([1 3 4],:)})    
anova_rm({PD_stability_mean([6:13],:), Control_stability_mean})    
anova_rm({PD_stability_mean, Control_stability_mean})    

anova2(PD_stability_mean([9,10,12],:))
anova1(Control_stability_mean)

    
    
    
    