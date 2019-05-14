% bumpParams.bumpDir = [0 90 180 270];
% spindleParams.trialInd = getBumpTrials(trial_data,bumpParams);
% spindleParams.trialInd = spindleParams.trialInd(:);




emgNames = trial_data.emg_names;
emgNames(8) = []; %get rid of trapezius
musNamesFull = trial_data.muscle_names;

musNames = musNamesFull([4, 18, 19, 20, 8, 9, 10, 24, 23, 35, 38, 37, ... 
    39, 6, 13, 14, 15, 28, 29, 5, 3]);

%% Set up params structure array


for i = 1:numel(emgNames) %Not including Trapezius
    spindleParams(i).trialInd = 1:numel(trial_data);
    spindleParams(i).bufferSize = 0; %This was taken care of in TD procoessing
    spindleParams(i).time_step = 0.005;

    spindleParams(i).startIdx = {'idx_startTime',0}; %reference index and relative idx
    spindleParams(i).endIdx = {'idx_endTime',0};  
    spindleParams(i).dataStore = 'lean';
    spindleParams(i).emgName = emgNames{i};
    spindleParams(i).musName = musNames{i};

end

%%
spindleData = struct('SimOut',[],'metaParams',[]);
for mus = 1:length(spindleParams)

    disp([spindleParams(mus).emgName])

    
    spindleData(mus).SimOut = getAffPotFromMusState(trial_data,spindleParams(mus));

    spindleData(mus).metaParams = spindleParams(mus);

end
