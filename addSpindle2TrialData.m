% bumpParams.bumpDir = [0 90 180 270];
% spindleParams.trialInd = getBumpTrials(trial_data,bumpParams);
% spindleParams.trialInd = spindleParams.trialInd(:);

filename = 'C:\Users\kpb8927\data\td-library\Han_20171106_TRT_5ms.mat';

load(filename)


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

%% Run spindle model on muscle kinematics data
tic;
spindleData = struct('SimOut',[],'metaParams',[]);
parfor mus = 1:length(spindleParams)

    disp([spindleParams(mus).emgName])

    
    spindleData(mus).SimOut = getAffPotFromMusState(trial_data,spindleParams(mus));

    spindleData(mus).metaParams = spindleParams(mus);
end

toc;

%% Rearrange spindle model output and put into TD
temp = trial_data;
% temp.spindle = [];

parfor trial = 1:numel(trial_data)
    for mus = 1:numel(spindleData)
        temp(trial).spindle(:,mus) = spindleData(mus).SimOut(trial).r;
    end
end

trial_data = temp;

%% Save file

save(filename,trial_data)



