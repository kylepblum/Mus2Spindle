function bumpTrialIdx = getBumpTrials(trial_data,params)

if numel(params.bumpDir) == 1
    bumpTrialIdx = trial_data(1).trialID(trial_data(1).bumpDir==params.bumpDir);
else
    bumpTrialIdx = [];
    for i = 1:numel(params.bumpDir)
        bumpDir = params.bumpDir(i);
        bumpTrialIdx = [bumpTrialIdx trial_data(1).trialID(trial_data(1).bumpDir==bumpDir)];
    end
end
end