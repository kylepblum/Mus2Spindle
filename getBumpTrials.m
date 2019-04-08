function bumpTrialIdx = getBumpTrials(trial_data,params)

if numel(params.bumpDir) == 1
    try
        bumpTrialIdx = trial_data(1).trialID(trial_data(1).bumpDir==params.bumpDir);
    catch
        bumpTrialIdx = trial_data(1).trial_id(trial_data(1).bumpDir==params.bumpDir);
    end
else
    bumpTrialIdx = [];
    for i = 1:numel(params.bumpDir)
        bumpDir = params.bumpDir(i);
        if isfield(trial_data,'trialID')
            bumpTrialIdx = [bumpTrialIdx trial_data(1).trialID(trial_data(1).bumpDir==bumpDir)];
        end
    end
end
end