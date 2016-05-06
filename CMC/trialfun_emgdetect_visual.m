function [trl] = trialfun_emgdetect_visual(cfg, varargin);
 
    % read the header and determine the channel number corresponding with the EMG

    emgDat = varargin{1};
    hdr=emgDat.hdr;
    
    tok=tokenize(emgDat.cfg.previous.previous.datafile,'.');
    disp('#############')
    disp(tok{1})
    disp('#############')
%     thresh = cfg.thresh;
    disp(['#####  Using thresh value: ', num2str(cfg.thresh)]);
%     hdr         = ft_read_header(cfg.headerfile);

%     chanindx    = match_str(hdr.label,(ft_channelselection('EMG*', hdr.label)));%strmatch('EMG*', hdr.label);
%     if length(chanindx)>1
%       error('only one EMG channel supported');
%     end

    % read all data of the EMG channel, assume continuous file format
%     emg = ft_read_data(cfg.datafile, 'header', hdr, ...
%                     'begsample', 1, 'endsample', hdr.nSamples*hdr.nTrials, ...
%                     'chanindx', chanindx, 'checkboundary', false);

    % apply filtering, hilbert transformation and boxcar convolution (for smoothing)
    emgflt      = ft_preproc_highpassfilter(emgDat.trial{1}, hdr.Fs, 10); % highpassfilter
    emghlb      = abs(hilbert(emgflt')');                     % hilbert transform
    %     emgcnv      = conv2([1], ones(1,hdr.Fs/2), emghlb, 'same'); % smooth using convolution
    emgstd      = ft_preproc_standardize(emghlb, 2);          % z-transform, i.e. mean=0 and stdev=1

    % plotting
    
    [peaks, locs] = findpeaks(emgstd(hdr.Fs*cfg.shift+1:(length(emgstd)-hdr.Fs*cfg.shift)),...
        'MinPeakHeight', cfg.thresh, 'MinPeakDistance',hdr.Fs*cfg.shift*(1-cfg.overlap));
    locs=locs+499;

        disp('##################')
        disp(['Epochs found: ', num2str(length(locs))]);
        disp('##################')

    time_vector = 0:1./hdr.Fs:(length(emgstd)-1)/hdr.Fs;
    hFig=figure(1);
    set(hFig, 'Position', [100 120 1656 968]);
    subplot 211
    plot(time_vector,emgstd, 'b'), hold on, plot(locs./hdr.Fs, peaks, '*r');
    plot([0 time_vector(end)], [cfg.thresh cfg.thresh], '--m'), 
    xlabel('Time(s)'), ylabel('z-trans(hilb(EMG))')
    legend('EMG-hilbert-std', 'Peaks', 'Threshold')
    title('Full data')

    subplot 212
    plot(time_vector,emgstd, 'b'), hold on, plot(locs./hdr.Fs, peaks, '*r');
    plot([0 time_vector(end)], [cfg.thresh cfg.thresh], '--m')
    xlim([10 14]), xlabel('Time(s)'), ylabel('z-trans(hilb(EMG))')
    legend('EMG-hilbert-std', 'Peaks', 'Threshold')
    title('Zoomed data')
    suptitle(['Found ', num2str(length(locs)), ' epochs'])
% 
    if strcmp(cfg.savepdf,'yes')
        export_fig( gcf,['J:\MEG_Research\CMC\raw_ArtRej\EMGthresh\',...
            num2str(randn()),'-EMGthresh'] ,...
            '-transparent', '-painters','-pdf', '-r250' ); 
    end
%     close(hFig)

%%
% 
%     emgtrl      = diff(emgtrl, [], 2);
% 
%     emgon       = find(emgtrl(:)== 1);
%     emgoff      = find(emgtrl(:)==-1);
    
    trl(:,1) = locs - round(hdr.Fs*cfg.shift)./2; %emgon (:) + hdr.Fs*0.5;  % as a consequence of the convolution with a one-second boxcar
    trl(:,2) = locs + round(hdr.Fs*cfg.shift)./2; %emgoff(:) - hdr.Fs*0.5;  % as a consequence of the convolution with a one-second boxcar
    trl(:,3) = 0;
       
return