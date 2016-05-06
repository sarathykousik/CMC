function pow = powCalc2(data)

% cfg_l             = [];
% cfg_l.method      = 'mtmfft';
% cfg_l.foi         = [1:300];
% cfg_l.taper       = 'hanning';
% cfg_l.channel     = 'MEGGRAD'
% powL  = ft_freqanalysis(cfg_l, data);



cfg_l             = [];
cfg_l.method      = 'mtmfft';
cfg_l.foi         = [1:200];
cfg_l.taper       = 'hanning';
cfg_l.channel     = 'MEGGRAD'
pow  = ft_freqanalysis(cfg_l, data);

% cfg = [];
% cfg.parameter = 'powspctrm';
% pow=ft_appendfreq(cfg,pow_L, pow_H)

return