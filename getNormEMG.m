function trial_data = getNormEMG(trial_data,params)
%
%

emgNorm = normalize(trial_data.emg,'zscore','robust');

trial_data.emgNorm = emgNorm;


end