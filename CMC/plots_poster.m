%% grand avg CMC

cfg = []; 
cfg.parameter = 'powspctrm'

for loop = 1:7
   
    PD_CMC_grndAvg{loop} = ft_freqgrandaverage(cfg, PDcoh_cmb{5:10,loop});
    Ctrl_CMC_grndAvg{loop} = ft_freqgrandaverage(cfg, Ctrlcoh_cmb{:,loop});
    
end

%% CMC 
% Left
L_par_ch = {'MEG0412+0413', 'MEG0422+0423', 'MEG0432+0433', 'MEG0442+0443',...
    'MEG0632+0633', 'MEG0712+0713', 'MEG0742+0743', 'MEG1812+1813', 'MEG1822+1823'};

% Right
R_par_ch = {'MEG0722+0723', 'MEG0732+0733', 'MEG1042+1043', 'MEG1112+1113',...
    'MEG1122+1123', 'MEG1132+1133', 'MEG1142+1143', 'MEG2212+2213', 'MEG2222+2223'};


