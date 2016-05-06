%%
clc;    clear all;     close all;
ft_defaults;
addpath('C:\Program Files\MATLAB\R2012b\toolbox_add_on\SEF')
proc                             = [];
proc.data_folder                 = 'J:\MEG_Research\CMC\for test';
proc.save_folder                 = 'J:\MEG_Research\CMC\for test\cleanData';
cd(proc.data_folder)
filenames      = dir('*.fif');

%%

for condLoop  = 2:length(filenames)%[1,2,4:length(fileList)]
    
    [a b c]     =   fileparts(filenames(condLoop).name);

    disp('#######################################')
    disp(['****    ', b ,'        ****'])
    disp('#######################################')
    
    % Import file
    cfg                         = [];
    cfg.dataset                 = filenames(condLoop).name;%filenames(condLoop).name;
    proc.data_import            = ft_preprocessing(cfg);
    
    cfg                         = [];
    cfg.dataset                 = filenames(condLoop).name;
    proc.data_import            = ft_preprocessing(cfg);

    cfg                         = [];
    cfg.lpfilter                = 'yes';
    cfg.lpfreq                  =  200;
    cfg.hpfilter                = 'yes';
    cfg.hpfreq                  =  2;
    cfg.dftfilter               = 'yes';
    cfg.channel                 =  {'MEG'}; 
    cfg.detrend                 = 'yes';
    proc.preproc_data_MEG       = ft_preprocessing(cfg, proc.data_import);

    cfg.channel                 =  {'EOG*'}; 
    proc.preproc_data_EOG       = ft_preprocessing(cfg, proc.data_import);
    cfg.channel                 =  {'ECG*'}; 
    proc.preproc_data_ECG       = ft_preprocessing(cfg, proc.data_import);

    % Epoch data
    time_tot    = floor(proc.data_import.time{1}(end));
    fs          = proc.data_import.fsample;
    sample_tot  = time_tot*fs;
    trl         = [];
    trl         = [[0:fs:sample_tot-fs]' [fs:fs:sample_tot]'];
    trl(:,3)    = 0;%[fs.*size(trl,1)]';
    trl(:,1)    = trl(:,1)+1;

    cfg = [];
    cfg.trl = trl;
%     proc.data_epoched = ft_redefinetrial(cfg, proc.preproc_data_ECG);

    proc.data_epoched           = ft_redefinetrial(cfg,proc.preproc_data_MEG);
    proc.data_epoched_EOG       = ft_redefinetrial(cfg,proc.preproc_data_EOG);
    proc.data_epoched_ECG       = ft_redefinetrial(cfg,proc.preproc_data_ECG);
    proc.data_epoched           = ft_appenddata([],...
        proc.data_epoched, proc.data_epoched_EOG,proc.data_epoched_ECG);
    

    % Jump artefact
    cfg                                 = [];
    cfg.trl = trl;
    cfg.continuous                      = 'no';
    % cutoff and padding
    cfg.artfctdef.zvalue.channel        = 'MEG';
    cfg.artfctdef.zvalue.cutoff         = 30;
    cfg.artfctdef.zvalue.trlpadding     = 0;
    cfg.artfctdef.zvalue.artpadding     = 0;
    cfg.artfctdef.zvalue.fltpadding     = 0;
    % algorithmic parameters
    cfg.artfctdef.zvalue.cumulative     = 'yes';
    cfg.artfctdef.zvalue.medianfilter   = 'yes';
    cfg.artfctdef.zvalue.medianfiltord  = 9;
    cfg.artfctdef.zvalue.absdiff        = 'yes';
    [~, artifact_jump]                  = ft_artifact_zvalue(cfg, proc.data_epoched);

%     % EOG
%     cfg                                 = [];
%     cfg.trl = trl;
%     cfg.continuous                      = 'no'; 
%     % cutoff and padding
%     cfg.artfctdef.zvalue.channel        = 'EOG*';
%     cfg.artfctdef.zvalue.cutoff         = 4;
%     cfg.artfctdef.zvalue.trlpadding     = 0;
%     cfg.artfctdef.zvalue.artpadding     = 0.1;
%     cfg.artfctdef.zvalue.fltpadding     = 0;
%     % algorithmic parameters
%     cfg.artfctdef.zvalue.bpfilter       = 'yes';
%     cfg.artfctdef.zvalue.bpfilttype     = 'but';
%     cfg.artfctdef.zvalue.bpfreq         = [1 15];
%     cfg.artfctdef.zvalue.bpfiltord      = 4;
%     cfg.artfctdef.zvalue.hilbert        = 'yes';
%     [~, artifact_EOG]                   = ft_artifact_zvalue(cfg, proc.data_epoched);

    cfg                                 = [];
    cfg.artfctdef.jump.artifact         = artifact_jump;
%     cfg.artfctdef.eog.artifact          = artifact_EOG;
    cfg.artfctdef.reject                = 'complete';
    cfg.artfctdef.crittoilim            = [-0.100 0.300];
    dataArtRej                          = ft_rejectartifact(cfg,proc.data_epoched);

    %
    tok_name = tokenize(b, '_');

    saveData.file                       = filenames(condLoop).name;
    saveData.subject                    = [tok_name{1},'_',tok_name{2}]; 
    saveData.condName                   = condLoop; 
    saveData.dataArtRej                 = dataArtRej;
    saveData.origTrialLen               = length(proc.data_epoched.trial);

    % Saving 
    %     save ([proc.save_folder,'\', subj,'\',b, '-SaveData'], 'saveData', '-v7.3');
    %     disp('*********************************')
    %     disp(['Saved data to:   ', [proc.save_folder,'\', subj,'\',b, '-SaveData']]);
    %     disp('*********************************')

    save ([proc.save_folder,'\',b, '-SaveData'], 'saveData', '-v7.3');
    disp('*********************************')
    disp(['Saved data to:   ', [proc.save_folder,'\',b, '-SaveData']]);
    disp('*********************************')
    proc.data_import = [];
    proc.data_trial= [];
    proc.preproc_data_MEG= [];
    proc.preproc_data_EOG= [];
    proc.data_epoched    = [];
    proc.data_epoched_EOG= [];
    artifact_jump= [];
    artifact_EOG= [];
    dataArtRej= [];
    saveData= [];

end


