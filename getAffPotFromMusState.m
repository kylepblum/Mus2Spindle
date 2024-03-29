function out = getAffPotFromMusState(in_data,params)

% params
%   .trialInd
%   .emgName
%   .musName
%   .bufferSize

if ~isfield(params,'trialInd')
    params.trialInd = 1:numel(in_data);
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

if ~isfield(params,'time_step')
    params.time_step = in_data(1).bin_size;
end

if ~isfield(params,'bumps')
    params.bumps = false;
end

try %Cheap way to deal with different nomenclature for these variables
    emgInd = find(strcmp(params.emgName,in_data(1).emg_names));
    musInd = find(strcmp(params.musName,in_data(1).muscle_names));
catch
    emgInd = find(strcmp(params.emgName,in_data(1).emg_names));
    musInd = find(strcmp(params.musName,in_data(1).musNames));
end



bs = params.bufferSize;
bin_size = in_data(1).bin_size;

% Do we need to interpolate?
if bin_size ~= params.time_step
    interpolate = true;
    time_step = params.time_step;
else
    interpolate = false;
    time_step = params.time_step;
end



%%%%%%%%%%%%%%%%


% spindlePool = parpool;
for a = 1

% Take out relevant data from trial a and pad so model can initialize

    thisTrial = a;
    if params.bumps
        if isnan(in_data(thisTrial).bumpDir(thisTrial))
            startIdx = in_data(thisTrial).idx_goCueTime - 10;
            endIdx = in_data(thisTrial).idx_goCueTime + 100;
        else
            startIdx = in_data(thisTrial).idx_bumpTime - 10;
            endIdx = in_data(thisTrial).idx_bumpTime + 100;
        end
        if isnan(startIdx) || isnan(endIdx)
            startIdx = 1;
            endIdx = 601;
        end
    else
        if numel(in_data)==1
            startIdx = 1;
            endIdx = size(in_data.musLenRel,1);
        else 
            startIdx = 1;
            endIdx = size(in_data(a).musLenRel,1);
        end
    end
    
    timeIdx = startIdx:endIdx;

    
    
    t = (-bs*bin_size:bin_size:(endIdx-startIdx)*bin_size);
    

    %%% Get gamma and cdl variables, interpolate if necessary
    if strcmpi(params.gamma,'agc')
        gammaD = in_data(thisTrial).emgNorm(timeIdx,emgInd);
        gammaS = in_data(thisTrial).emgNorm(timeIdx,emgInd);
    elseif strcmpi(params.gamma,'pas')
        gammaD = 0.3*ones(size(in_data(thisTrial).emgNorm(timeIdx,emgInd)));
        gammaS = 0.3*ones(size(in_data(thisTrial).emgNorm(timeIdx,emgInd)));
    elseif strcmpi(params.gamma,'static')
        gammaD = 0.3*ones(size(in_data(thisTrial).emgNorm(timeIdx,emgInd)));
        gammaS = in_data(thisTrial).emgNorm(timeIdx,emgInd);    
    end
    %     gamma = [ones(bs,1)*gamma(1); gamma];
%     cdl = in_data(thisTrial).musLenRel(timeIdx,musInd)*1300;
%     cdl = [ones(bs,1)*cdl(1); cdl]';
    
    delta_cdl = -in_data(thisTrial).musVelRel(timeIdx,musInd)*1300*bin_size;
%     delta_cdl(1) = 1300; 

    % If we need to interpolate
    if interpolate == 1
        t_interp = (-bs*bin_size:time_step:(endIdx-startIdx)*bin_size);
        gammaS = interp1(t,gammaS,t_interp);
        gammaD = interp1(t,gammaD,t_interp);
%         cdl = interp1(t,cdl,t_interp);
        delta_cdl = interp1(t,delta_cdl,t_interp);
        t = t_interp;
    end
    
    %%% Process emg into gamma signal
%     gammaS = movmean(gammaS,0.03/(t(2)-t(1))); %30 ms smoothing filter
%     gammaD = movmean(gammaD,0.03/(t(2)-t(1))); %30 ms smoothing filter
    %if we want constant gamma, use initial value   
    
%     gammaD = gamma;
    gammaDinit = gammaD(1);
    delta_gammaD = diff(gammaD);
    delta_gammaD = [gammaDinit+0.1; delta_gammaD];
    
    delta_gammaD(gammaD>=1) = 0; %saturate gamma at 1
    
%     gammaS = gamma;
    gammaSinit = gammaS(1);
    delta_gammaS = diff(gammaS);
    delta_gammaS = [gammaSinit+0.1; delta_gammaS];
    
    delta_gammaS(gammaS>=1) = 0; %saturate gamma at 1
    
    % GET L0 as a struct field?
    % Initialize hs with real initial lengths!
    % develop driver to take bagParams/chainParams to initialize hs models
%     cdlInit = cdl(1)-1300;
%     delta_cdl = diff(cdl);
%     delta_cdl(end+1) = delta_cdl(end);
%     delta_cdl = in_data(thisTrial).musVelRel(timeIdx,musInd)*bin_size*1300; %transform velocity from L0/s into model units (nm/dt): 0.01 s/dt, 1300 nm/L0
%     delta_cdl = in_data(thisTrial).musVelRel(timeIdx,musInd)*bin_size*1300; %transform velocity from L0/s into model units (nm/dt): 0.01 s/dt, 1300 nm/L0
%     delta_cdl(1) = cdlInit;
%     delta_cdl = movmean(delta_cdl,10);
%     delta_cdl = delta_cdl;

    if sum(isnan(delta_cdl)) || sum(isnan(gammaD)) || sum(isnan(gammaS))
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
%             [r,rs,rd] = sarc2spindle(dataB(a),dataC(a),1,2,0.03,0,0.0);
            if strcmpi(params.dataStore,'lean')
%                 out(a).r = r(bs+1:end)';
%                 out(a).rs = rs(bs+1:end)';
%                 out(a).rd = rd(bs+1:end)';
                out(a).bin_size = time_step; % Rename to be consistent with TD
%                 out(a).dataB.f_activated = dataB(a).f_activated(bs+1:end)';
%                 out(a).dataC.f_activated = dataC(a).f_activated(bs+1:end)';
                out(a).dataB.hs_force = dataB(a).hs_force(bs+1:end)';
                out(a).dataC.hs_force = dataC(a).hs_force(bs+1:end)';
                out(a).dataB.cmd_length = dataB(a).cmd_length(bs+1:end)';
                out(a).dataC.cmd_length = dataC(a).cmd_length(bs+1:end)';
                out(a).dataB.time_step = time_step;
                out(a).dataC.time_step = time_step;
%                 out(a).dataB.t
%                 out(a).dataC.t = t(bs+1:end)';
            else
                out(a).r = r(bs+1:end);
                out(a).dataB = removeBufferFromStruct(dataB(a),bs);
                out(a).dataC = removeBufferFromStruct(dataC(a),bs);
                out(a).hsB = hsB(a);
                out(a).hsC = hsC(a);
                out(a).delta_cdl = delta_cdl(bs+1:end);
                out(a).trialInd = thisTrial;
                
            end
        end
%     in_data(a).spindle = out(a);  %Weird indexing because of parfor
    
%     out(a).gamma = gammaD(bs+1:end);
end
% delete(spindlePool);
end
