% bumpParams.bumpDir = [0 90 180 270];
% spindleParams.trialInd = getBumpTrials(trial_data,bumpParams);
% spindleParams.trialInd = spindleParams.trialInd(:);

spindleParams.trialInd = 1:numel(trial_data);



spindleParams.emgName = 'DeltMid';
spindleParams.musName = 'deltoid_med';
spindleParams.bufferSize = 1/trial_data(1).bin_size;

spindleParams.startIdx = {'idx_bumpTime',-100}; %reference index and relative idx
spindleParams.endIdx = {'idx_bumpTime',500};  
spindleParams.dataStore = 'lean';


out = getAffPotFromMusState(trial_data,spindleParams);

