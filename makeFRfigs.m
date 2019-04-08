params.savepath = '~/LimbLab/Projects/EMGanalysis/Figures/Han/Han_20171207_COactpas/';
params.filetype = '.pdf';
params.savefig = 0;
params.musIdx = 6;
params.trialType = 'bum';


for a = 4
    params.frToPlot = a;
    plotFR(trial_data,params)
end