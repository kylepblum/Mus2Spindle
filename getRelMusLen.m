function trial_data = getRelMusLen(trial_data,params)
%
% TO-DOs: 
% - Pennation angles?
% - 

L0 = nanmean(trial_data.muscle_len);
musLenRel = trial_data.muscle_len./L0;
musVelRel = trial_data.muscle_vel./L0;

trial_data.musLenRel = musLenRel;
trial_data.musVelRel = musVelRel;

end