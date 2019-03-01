function out = getAffPotFromMusState(trial_data,params)

% params
%   .trialInd
%   .emgName
%   .musName
%   .bufferSize

if ~isfield(params,'trialInd')
    params.trialInd = 1:numel(trial_data);
end

if ~isfield(params,'emgName') || ~isfield(params,'musName')
    params.musName = 'deltoid_ant';
    params.emgName = 'DeltAnt';
end

if ~isfield(params,'bufferSize')
    params.bufferSize = 100;
end

if ~isfield(params,'startIdx') || ~isfield(params,'endIdx')
    params.startIdx = {'idx_startTime',0};
    params.endIdx = {'idx_endTime',0};
end



emgInd = find(strcmp(params.emgName,trial_data(1).emg_names));
musInd = find(strcmp(params.musName,trial_data(1).muscle_names));




bs = params.bufferSize;
bin_size = trial_data(1).bin_size;




%%%%%%%%%%%%%%%%


tic
parfor a = 1:numel(params.trialInd)
    disp([params.emgName ': Beginning trial ' num2str(params.trialInd(a))])

% Take out relevant data from trial a and pad so model can initialize

    thisTrial = params.trialInd(a);
    startIdx = trial_data(thisTrial).(params.startIdx{1}) + params.startIdx{2};
    endIdx = trial_data(thisTrial).(params.endIdx{1}) + params.endIdx{2};
    if isnan(startIdx) || isnan(endIdx)
        startIdx = 1;
        endIdx = 601;
    end
    
    timeIdx = startIdx:endIdx;

    
    
    t = (-bs*bin_size:bin_size:(endIdx-startIdx)*bin_size);
    

    %%% Process emg into gamma signal
    gamma = trial_data(thisTrial).emgNorm(timeIdx,emgInd);
    gamma = smooth(gamma,5);
    gamma = 0.5*zeros(size(gamma));
    
    gammaD = gamma;
    gammaD = [ones(bs,1)*gammaD(1); gammaD];
    delta_gammaD = diff(gammaD);
    delta_gammaD = [0.5; delta_gammaD];
    
    gammaS = gamma;
    gammaS = [ones(bs,1)*gammaS(1); gammaS];
    delta_gammaS = diff(gammaS);
    delta_gammaS = [0.5; delta_gammaS];
    
    % GET L0 as a struct field?
    % Initialize hs with real initial lengths!
    % develop driver to take bagParams/chainParams to initialize hs models
    
    
    cdl = trial_data(thisTrial).musLenRel(timeIdx,musInd);
    delta_cdl = diff(cdl)*1300;
    delta_cdl(end+1) = delta_cdl(end);
%     delta_cdl = trial_data(thisTrial).musVelRel(timeIdx,musInd)*bin_size*1300; %transform velocity from L0/s into model units (nm/dt): 0.01 s/dt, 1300 nm/L0
%     delta_cdl = trial_data(thisTrial).musVelRel(timeIdx,musInd)*bin_size*1300; %transform velocity from L0/s into model units (nm/dt): 0.01 s/dt, 1300 nm/L0
    delta_cdl = [zeros(bs,1)*delta_cdl(1); delta_cdl];
%     delta_cdl = smooth(delta_cdl,10);
%     delta_cdl = delta_cdl;

    if sum(isnan(delta_cdl)) || sum(isnan(gamma))
       delta_cdl = zeros(size(delta_cdl));
       delta_gammaS = zeros(size(delta_cdl));
       delta_gammaD = zeros(size(delta_cdl));
       out(a).nanflag = 1;
    else
       out(a).nanflag = 0;
    end
       

% Use MATLAB's built-in parallel computing to run simulations and store
% data: 
        [hsB(a),dataB(a),hsC(a),dataC(a),flags(a)] = sarcSimDriver(t,delta_gammaD',delta_gammaS',delta_cdl');

        if flags(a).error == 1
            warning('Problem with simulation. Null data will be stored.')
            out(a).r = [];
            out(a).dataB = struct;
            out(a).dataC = struct;
            out(a).hsB = struct;
            out(a).hsC = struct;
            out(a).delta_cdl = [];
            out(a).trialInd = thisTrial;
            out(a).errorflag = 1;
        else 
            [r,~,~] = sarc2spindle(dataB(a),dataC(a),1,1,0.03,1,0.0);
            out(a).r = r(bs+1:end);
            out(a).dataB = removeBufferFromStruct(dataB(a),bs);
            out(a).dataC = removeBufferFromStruct(dataC(a),bs);
            out(a).hsB = hsB(a);
            out(a).hsC = hsC(a);
            out(a).delta_cdl = delta_cdl(bs+1:end);
            out(a).trialInd = thisTrial;

        end
%     trial_data(a).spindle = out(a);  %Weird indexing because of parfor
    
%     out(a).gamma = gammaD(bs+1:end);
end
toc
end
