params.savepath = '~/LimbLab/Projects/EMGanalysis/Figures/Han/Han_20171207_COactpas/';
params.filetype = '.pdf';
params.savefig = 0;
params.musIdx = 14;
params.trialType = 'bum';
for a = [1:11]
    params.emgToPlot = a;
    plotEMG(trial_data,params)
end