function trial_data = getNormEMG(trial_data,params)
%
%

emgNorm = normalize(trial_data.emg,'range');

% Deal with NaNs
t = (1:size(emgNorm,1))';

for mus = 1:size(emgNorm,2)
    nanEMG = isnan(emgNorm(:,mus));
    emgNorm(nanEMG,mus) = interp1(t(~nanEMG), emgNorm(~nanEMG,mus), t(nanEMG));
    %Special consideration for NaNs at the margins
    if nanEMG(1)
       firstNonNaN = find(~isnan(emgNorm(:,mus)),1);
       emgNorm(1:firstNonNaN-1,mus) = emgNorm(firstNonNaN,mus);
    end
    if nanEMG(end)
       finalNonNaN = find(~isnan(emgNorm(:,mus)),1,'last');
       emgNorm(end-finalNonNaN:end,mus) = emgNorm(finalNonNaN,mus);      
    end
end
  
% Add to TrialData structure
trial_data.emgNorm = emgNorm;


end