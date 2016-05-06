function trl = trialfun_EMG(cfg)

% Find total length of time and cut into 1 s

time_tot    = floor(proc.data_import.time{1}(end));
fs          = proc.data_import.fsample;
sample_tot  = time_tot*fs;
trl         = [];
trl         = [[0:fs:sample_tot-fs]' [fs:fs:sample_tot]'];
trl(:,3)    = [fs.*size(trlok,1)]';

end