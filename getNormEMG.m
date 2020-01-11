function trial_data = getNormEMG(trial_data,params)
%
%
emgNorm = normalize(trial_data.emg,'range');
% Deal with NaNs
emgNorm = fillmissing(emgNorm,'linear',1,'EndValues','nearest');
% Add to TrialData structure
trial_data.emgNorm = emgNorm;

end