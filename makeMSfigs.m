params.savepath = '~/LimbLab/Projects/EMGanalysis/Figures/Han/Han_20171207_COactpas/';
params.filetype = '.pdf';
params.savefig = 0;
params.musIdx = 9;




plotMS(trial_data,out_AGDyn_LowStc,params)
plotMS(trial_data,out_AGStc_LowDyn,params)
plotMS(trial_data,out_AGDyn_AGStc,params)
plotMS(trial_data,out_HighDyn_HighStc,params)
plotMS(trial_data,out_AGDyn_HighStc,params)
plotMS(trial_data,out_HighDyn_AGStc,params)

