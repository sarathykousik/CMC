
cd('J:\MEG_Research\CMC\for test\visClean')

filenames           = dir('*.mat');
 
for loop = 1:length(filenames)
    [a b c]     =   fileparts(filenames(loop).name);
    disp(['######## ', b])
    load(filenames(loop).name)

cfg              = [];
cfg.output       = 'pow';
cfg.method       = 'mtmfft';
cfg.foilim       = [1 200];
cfg.tapsmofrq    = 2;
cfg.keeptrials   = 'yes';
cfg.channel      = {'MEGGRAD'};
freq_on          = ft_freqanalysis(cfg, visClean);

cfg = [];
cfg.jackknife  = 'yes';
freq_desc    = ft_freqdescriptives(cfg,  freq_on);

freq_cmb{loop} = ft_combineplanar([], freq_desc);

end

%%
cfg                  = [];
cfg.parameter        = 'powspctrm';
cfg.ylim             = [0 20e-24];
cfg.maskstyle               = 'saturation';	
cfg.shading                 = 'interp';
cfg.layout      = 'neuromag306cmb.lay';
figure; ft_multiplotER(cfg,  freq_cmb{[1 2 3]})

%%

for loop = 1:6
   
    freq_arr(loop,:) = mean(freq_cmb{loop}.powspctrm,1);
    
    
end

%%

semilogy(freq_cmb{1}.freq, freq_arr([1 2 3],:))
legend('DBS ON', 'DBS OFF','MED ON')
xlabel('Frequency'), ylabel('Power')

title('Subj 1')

%%
Fs = proc.preproc_MEG.fsample;
chan = 100;
[Pxx,freq] = pwelch(proc.preproc_MEG.trial{1}(chan,:),...
                        hanning(ceil(Fs)),ceil(Fs*0.2),ceil(Fs), Fs);

                    
[y,r,vr]=ssan(proc.preproc_MEG.trial{1}(chan,:), 20);
[Pxx_ssa,freq_ssa] = pwelch(y,...
                        hanning(ceil(Fs)),ceil(Fs*0.2),ceil(Fs), Fs);


figure(10)                    
semilogy(freq(1:200),Pxx(1:200));
hold on
semilogy(freq_ssa(1:200),Pxx_ssa(1:200), 'r');







