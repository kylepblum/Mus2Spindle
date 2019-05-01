params.savepath = '~/LimbLab/Projects/EMGanalysis/Figures/Han/Han_20171207_COactpas/';
params.filetype = '.pdf';
params.savefig = 0;
params.musIdx = 6;
params.trialType = 'bum';
params.spikeArray = 'cuneate_spikes';


for a = 1:19
    params.frToPlot = a;
    plotFR(td,params)
end