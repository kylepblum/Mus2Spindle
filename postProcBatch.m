% This is to get td files that weren't run through my procBatch code into a
% similar state for analysis

%% Normalize muscle lengths
% This is probably best trial-to-trial anyway (normalize to first few ms)

% Chris's opensim data is in its own structure: which to normalize?
lenParams.opensimChris = 1;
lenParams.idx_opensimLen = 15:53;
lenParams.idx_opensimVel = 54:92;
lenParams.L0 = 'trial_init';

% trial_data(1).musLenRel = [];
% trial_data(1).musVelRel = [];
% trial_data(1).muscle_names = {};

for trial = 1:numel(td)
    temp{trial} = getRelMusLen(td(trial),lenParams);
    trial_data(trial) = temp{trial};
end

%% Normalize EMGs
%This is probably best as one big dataset (normalizing to trial = bad)

trial_data(1).emgNorm = [];
tempEMG = struct('emg',[]); %need to do this to mimic trial_data structure

% Get all emg data into a single array
for trial = 1:numel(trial_data)
    idxTrialStart(trial) = length(tempEMG.emg) + 1;
    tempEMG.emg = [tempEMG.emg; trial_data(trial).emg];
    idxTrialEnd(trial) = length(tempEMG.emg);
end

tempEMG = getNormEMG(tempEMG,[]);

for trial = 1:numel(trial_data)
    trial_data(trial).emgNorm = tempEMG.emg(idxTrialStart(trial):idxTrialEnd(trial),:);
end