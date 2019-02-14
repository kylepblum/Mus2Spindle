function trial_data = getAffPotFromMusState(trial_data,params)

% params
%   .trialInd
%   .emgName
%   .musName
%   .bufferSize

if ~isfield(params,'trialInd')
    params.trialInd = length(trial_data);

elseif ~isfield(params,'.emgName')
    params.musName = 'deltoid_ant';
    params.emgName = 'DeltAnt';

elseif ~isfield(params,'.musName')
    warning('No muscle names were given in params struct.')
    params.musName = 'deltoid_ant';
    params.emgName = 'DeltAnt';

end




% TODOs
% EMG normalization?
% Muscle length normalization?
% 

emgInd = find(strcmp(params.emgName,td_s(1).emg_names));
musInd = find(strcmp(params.musName,td_s(1).muscle_names));
bs = params.bufferSize;





%%%%%%%%%%%%%%%%
numSims = numel(trial_data);       % Number of simulations to run in parallel



parfor a = 1:numSims
% Take out relevant data from trial a and pad so model can initialize
    t = (0:trial_data(a).bin_size:trial_data(a).idx_endTime)';
    delta_gammaD = trial_data(a).emg(:,emgInd);
    delta_gammaD = [ones(bs,1)*delta_gammaD(1);...
                    delta_gammaD;
                    ones(bs,1)*delta_gammaD(end)];
    
    delta_gammaS = trial_data(a).emg(:,emgInd);
    delta_gammaS = [ones(bs,1)*delta_gammaS(1);...
                    delta_gammaS;
                    ones(bs,1)*delta_gammaS(end)];
    
    delta_cdl = trial_data(a).muscle_len(musInd);


    
    
    
% Use MATLAB's built-in parallel computing to run simulations and store
% data: 
    [hsB(a),dataB(a),hsC(a),dataC(a)] = sarcSimDriver(t,delta_gammaD,delta_gammaS,delta_cdl(a,:));
    [r(a),~,~] = sarc2spindle(dataB(a),dataC(a),1,1,0.03,1,0.09);
    disp(['Done with simulation number ' num2str(a)])
end
