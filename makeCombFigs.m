params.savepath = 'C:\Users\kpb8927\Documents\data\cuneate\Lando_COactpas_20170917_TD_001\Figures\NewChain\';
params.filetype = '.pdf';
params.savefig = 0;
params.muscles = {'brachioradialis','brachialis','tricep_lat','tricep_lon',...
    'tricep_sho','pectoralis_sup','pectoralis_inf','lat_dorsi_sup',...
    'flex_carpi_ulnaris','flex_carpi_radialis','flex_digit_superficialis',...
    'deltoid_ant','deltoid_med','deltoid_pos','teres_major','infraspinatus'...
    'ext_carp_rad_brevis','ext_carpi_ulnaris','ext_digitorum','bicep_sh',...
    'bicep_lh'};
params.emgToPlot =8;
params.trialType = 'bum';



params.musIdx = 3;


% plotMS(trial_data,out_AGDyn_LowStc,params)
% plotMS(trial_data,out_AGStc_LowDyn,params)
% plotMS(trial_data,out_AGDyn_AGStc,params)
% plotMS(trial_data,out_HighDyn_HighStc,params)
% plotMS(trial_data,out_AGDyn_HighStc,params)
% plotMS(trial_data,out_HighDyn_AGStc,params)
plotCombinedSignals(trial_data,spindleOut_constG(:,params.musIdx),params)


