
data_folder = 'J:\MEG_Research\ArtRejtSSS\results';


cd(data_folder)
filenames      = dir('*.mat');




%%
for loop = 1:length(filenames)
    load(filenames(loop).name)
        [a b c]     =   fileparts(filenames(condLoop).name);

    disp('#######################################')
    disp(['****    ', b ,'        ****'])
    disp('#######################################')
    
    data = saveData.dataArtRej;
    powCalc()
    
    for trlLoop = 1:length(data{1,loop}.dataArtRej.trial)
        dat_fT{1,loop}.dataArtRej.trial{trlLoop} = data{1,loop}.dataArtRej.trial{trlLoop}./1e-12.*10;
    end
    
end

%% Power analysis
function output = powCalc()

cfg             = [];
cfg.method      = 'mtmfft';
cfg.foi         = [40:200];
cfg.taper       = 'dpss';
cfg.tapsmofrq   = 0.1*cfg.foi;
cfg.channel     = 'MEGGRAD'

cfg_l             = [];
cfg_l.method      = 'mtmfft';
cfg_l.foi         = [1:39];
cfg_l.taper       = 'hanning';
cfg_l.channel     = 'MEGGRAD'



for loop = 1:2
    powG{loop} = ft_freqanalysis(cfg,   dat_fT{1,loop}.dataArtRej);
    powL{loop} = ft_freqanalysis(cfg_l, dat_fT{1,loop}.dataArtRej);
    dat_fT{1,loop}.file   
end

%%
cfg = [];
cfg.parameter = 'powspctrm';
for loop = 1:2
    pow{loop} = ft_appendfreq(cfg, powL{loop}, powG{loop});
end

%%
clrOp1 = rgb('dodgerBlue');%[0.7 0.2 0.2];
clrOp2 = rgb('DimGray');%[0.4 0.4 0.4];
% figure,
semilogy(pow{1}.freq, pow{1}.powspctrm(20,:), 'color',clrOp1,'linewidth', 4), hold on
semilogy(pow{2}.freq, pow{2}.powspctrm(20,:), 'color',clrOp2,'linewidth', 4);

hlegend = legend({'tSSS', 'no tSSS'});
% hlegend = legend({'UPDRS III'},'Interpreter','latex');
set(hlegend,'Fontsize',18);
set(hlegend,'Fontname','Helvetica');
set(hlegend,'Location','Northeast');
legend boxoff
box off

set(gca,'XTick',[0:50:200]);
set(gca,'Fontsize', 18)
set(gca,'Fontweight','bold');
set(gca,'Fontname','Helvetica');

hxlabel = xlabel('Frequency Hz','Interpreter', 'latex');
set(hxlabel,'Fontsize',24);
set(hxlabel,'FontWeight','bold');
% set(hxlabel,'Fontname','Helvetica');

hylabel = ylabel(['Spectral density  fT/cm- $$\sqrt{Hz}$$'],'Interpreter', 'latex');
set(hylabel,'Fontsize',24);
set(hylabel,'FontWeight','bold');
% set(hylabel,'Fontname','Helvetica');
% export_fig( gcf, 'PowPlots','-transparent', ...
%         '-painters','-pdf', '-r250' ); 

%%    

cfg_sel = [];
cfg_sel.trials = [240];
trlData{1,1} = ft_selectdata(cfg_sel, dat_fT{1,1}.dataArtRej);
trlData{1,2} = ft_selectdata(cfg_sel, dat_fT{1,2}.dataArtRej);

%%
cfg                 = [];
cfg.layout          = 'neuromag306all.lay';
cfg.channel         = 'MEGGRAD';
cfg.hlim            = [-0.1 0.2];
cfg.vlim            = [-2000 2000];
figure,
ft_multiplotER(cfg,trlData{1,1},trlData{1,2});


hlegend = legend({'tSSS', 'no tSSS'});
set(hlegend,'Fontsize',18);
set(hlegend,'Fontname','Helvetica');
set(hlegend,'Location','Northeast');
legend boxoff
box off

export_fig( gcf, 'multiER-full-LH','-transparent', ...
        '-painters','-pdf', '-r250' ); 

%% Timelock analysis

cfgtlck = [];
cfgtlck.channel = 'MEGGRAD';
cfgtlck.removemean = 'yes';

for loop = 1:2
    tlck{loop}     = ft_timelockanalysis(cfgtlck,dat_fT{loop}.dataArtRej);
    tlckCmb{loop}  = ft_combineplanar([], tlck{loop});
end















