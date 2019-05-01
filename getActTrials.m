function actTrialIdx = getActTrials(trial_data,params)

if numel(params.targDir) == 1
    if isfield(trial_data,'target_direction')
        
        try
            actTrialIdx = trial_data(1).trialID(trial_data(1).target_direction==params.targDir...
                & isnan(trial_data(1).bumpDir));
        catch
            actTrialIdx = trial_data(1).trial_id(trial_data(1).target_direction==params.targDir...
                & isnan(trial_data(1).bumpDir));
        end
        
    elseif isfield(trial_data,'tgtDir')
        
        try
            actTrialIdx = trial_data(1).trialID(trial_data(1).tgtDir==params.targDir...
                & isnan(trial_data(1).bumpDir));
        catch
            actTrialIdx = trial_data(1).trial_id(trial_data(1).tgtDir==params.targDir...
                & isnan(trial_data(1).bumpDir));
        end
        
    end
else
    actTrialIdx = [];
    for i = 1:numel(params.bumpDir)
        targDir = params.targDir(i);
        if isfield(trial_data,'trialID')
            if isfield(trial_data,'target_direction')
                actTrialIdx = [actTrialIdx trial_data(1).trialID(trial_data(1).target_direction==targDir...
                    & isnan(trial_data(1).bumpDir))];
            elseif isfield(trial_data,'tgtDir')
                actTrialIdx = [actTrialIdx trial_data(1).trialID(trial_data(1).tgtDir==targDir...
                    & isnan(trial_data(1).bumpDir))];
                
            end
        end
    end
end
end