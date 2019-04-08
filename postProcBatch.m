% This is to get td files that weren't run through my procBatch code into a
% similar state for analysis

%% Get rid of non-reward trials
% clearvars -except td

for trial = 1:numel(td)
    if strcmpi(td(trial).result,'R')
        tempTD(trial) = td(trial);
    end
end
trial_data = tempTD;
%% Normalize muscle lengths
% 

% Chris's opensim data is in its own structure: which to normalize?
lenParams.opensimChris = 1;
lenParams.idx_opensimLen = 15:53;
lenParams.idx_opensimVel = 54:92;
lenParams.L0 = 'session_mean';
lenParams.L0idx = 150;

trial_data(1).musLenRel = [];
trial_data(1).musVelRel = [];
trial_data(1).musNames = [];


% trial_data(1).musLenRel = [];
% trial_data(1).musVelRel = [];
% trial_data(1).muscle_names = {};
if strcmpi(lenParams.L0,'session_mean')
    tempKin = struct('opensim',[],'opensim_names',[]);
    for trial = 1:numel(td)
        idxTrialStart(trial) = length(tempKin.opensim) + 1;
        tempKin.opensim = [tempKin.opensim; td(trial).opensim];
        idxTrialEnd(trial) = length(tempKin.opensim);
    end
    tempKin.opensim_names = td(1).opensim_names;
    tempKin = getRelMusLen(tempKin,lenParams);
    for trial = 1:numel(trial_data)
        trial_data(trial).musLenRel = tempKin.musLenRel(idxTrialStart(trial):idxTrialEnd(trial),:);
        trial_data(trial).musVelRel = tempKin.musVelRel(idxTrialStart(trial):idxTrialEnd(trial),:);
        trial_data(trial).musNames = tempKin.musNames;
    end
else
    for trial = 1:numel(td)
        temp{trial} = getRelMusLen(trial_data(trial),lenParams);
        trial_data(trial) = temp{trial};
    end
end
%% Normalize EMGs
%This is probably best as one big dataset (normalizing to trial = bad)

trial_data(1).emgNorm = [];
tempEMG = struct('emg',[]); %need to do this to mimic trial_data structure

% Get all emg data into a single array
for trial = 1:numel(td)
    idxTrialStart(trial) = length(tempEMG.emg) + 1;
    tempEMG.emg = [tempEMG.emg; td(trial).emg];
    idxTrialEnd(trial) = length(tempEMG.emg);
end

tempEMG = getNormEMG(tempEMG,[]);

for trial = 1:numel(trial_data)
    trial_data(trial).emgNorm = tempEMG.emgNorm(idxTrialStart(trial):idxTrialEnd(trial),:);
end

%% Misc (mostly collect certain variables together for every trial)
% Get 'trialID' for every trial into each trial struct
if ~isfield(trial_data,'trialID')

        trialID = 1:numel(trial_data);

    for trial = 1:numel(trial_data)
        trial_data(trial).trialID = trialID;
    end
    trial_data = rmfield(trial_data,'trial_id');
end

% Get 'bumpDir' and 'target_direction' for every trial into each trial struct
if numel(trial_data(1).bumpDir) == 1
    bumpDir = [];
    for trial = 1:numel(trial_data)
        bumpDir = [bumpDir; trial_data(trial).bumpDir];
    end
    for trial = 1:numel(trial_data)
        trial_data(trial).bumpDir = bumpDir;
    end
end
if numel(trial_data(1).target_direction) == 1
    target_direction = [];
    for trial = 1:numel(trial_data)
        target_direction = [target_direction; trial_data(trial).target_direction];
    end
    for trial = 1:numel(trial_data)
        trial_data(trial).target_direction = target_direction;
    end
end   
    
%% Get smooth firing rates
% smoothParams.signals = 'cuneate_spikes';
smoothParams.signals = 'cuneate_spikes';
smoothParams.kernel_SD = 0.05;
smoothParams.calc_rate = true;
trial_data = smoothSignals(trial_data,smoothParams);
    