function pow = powCalc(data)

% cfg_l             = [];
% cfg_l.method      = 'mtmfft';
% cfg_l.foi         = [1:300];
% cfg_l.taper       = 'hanning';
% cfg_l.channel     = 'MEGGRAD'
% powL  = ft_freqanalysis(cfg_l, data);

cfg             = [];
cfg.method      = 'mtmfft';
cfg.foi         = [40:200];
cfg.taper       = 'dpss';
cfg.tapsmofrq   = 0.1*cfg.foi;
cfg.channel     = 'MEGGRAD'
pow_H  = ft_freqanalysis(cfg, data);

cfg_l             = [];
cfg_l.method      = 'mtmfft';
cfg_l.foi         = [1:39];
cfg_l.taper       = 'hanning';
cfg_l.channel     = 'MEGGRAD'
pow_L  = ft_freqanalysis(cfg_l, data);

cfg = [];
cfg.parameter = 'powspctrm';
pow=ft_appendfreq(cfg,pow_L, pow_H)

return