function actTrialIdx = getActTrials(trial_data,params)

if numel(params.bumpDir) == 1
    try
        actTrialIdx = trial_data(1).trialID(trial_data(1).target_direction==params.targDir...
            & isnan(trial_data(1).bumpDir));
    catch
        actTrialIdx = trial_data(1).trial_id(trial_data(1).target_direction==params.targDir...
            & isnan(trial_data(1).bumpDir));
    end
else
    actTrialIdx = [];
    for i = 1:numel(params.bumpDir)
        targDir = params.targDir(i);
        if isfield(trial_data,'trialID')
            actTrialIdx = [actTrialIdx trial_data(1).trialID(trial_data(1).target_direction==targDir...
                & isnan(trial_data(1).bumpDir))];
        end
    end
end
end