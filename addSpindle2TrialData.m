%% Load and preprocess trial_data

clear, clc

td_path = '/Users/kylepblum/Koofr/Limblab/Data/S1-adaptation/';
biomch_path = fullfile(td_path,'biomech');
spdl_path = fullfile(td_path,'spindle');


td_name_pre = 'Duncan_20190911*5ms*';
% td_name_pre = 'Han_20200113*5ms*';
filenames = dir(fullfile(td_path,td_name_pre));
for i = 1:numel(filenames)
    temp(i) = load(fullfile(td_path,filenames(i).name));
end

% Note -- they were loaded in alphabetically 
temp(1).trial_data.epoch = 'AD';
temp(2).trial_data.epoch = 'BL';
temp(3).trial_data.epoch = 'WO';

splitParams.split_idx_name = 'idx_startTime';
splitParams.linked_fields = {'result','trialID','bumpDir','tgtDir'};


td = [];

for i = [2, 1, 3] %files are loaded in alphabetical order -- this rearragnes them
    td = [td splitTD(temp(i).trial_data,splitParams)];
end

td = td(strfind([td.result],'R'));

td = td(1:end-1); %Last trial is wonky sometimes... 


clearvars -except td td_path biomch_path spdl_path td_name_pre
    
%% Setup biomech simulation
clc

td_biomech = td;
% td_biomech = td(7);

% for iTrial = 1:length(td_biomech)
%     td_biomech(iTrial).force = movmean(td_biomech(iTrial).force,30);
% end

body_weight = 10; % body mass in kg 11.7
M1 = 34.4*body_weight/1000; % from Scott 2000
M2 = 25.2*body_weight/1000;
L1 = 0.24;
L2 = 0.22; 
dt = td_biomech(1).bin_size;
left_arm = 0;

biomechParams = struct('M1',M1,'M2',M2,'L1',L1,'L2',L2,'left_arm',left_arm,'dt',dt);
biomechParams.muscle_names = {'sh_flx','sh_ext','el_flx','el_ext','bi_flx','bi_ext'};
% clearvars -except td sim_data td_biomech biomechParams td_path biomch_path td_name_pre

%% Run biomech simulation
sim_data = biomechPlanarArm(td_biomech,biomechParams);

% close all
% for i = 1
% figure; plot(sim_data(i).muscle_L);
% end

% Save biomech simulation
save_sim = 1;

if save_sim
    save(fullfile(biomch_path,[td_name_pre(1:end-5) '_sim_data_FV_5ms']),'sim_data')
end
% clearvars -except td sim_data biomechParams td_path biomch_path td_name_pre
%% Plot the simulation

% close all;
%
%     params.type = 'time_signals';
biomechParams.type = 'freeze_video';
biomechParams.resolution = 5;
biomechParams.signals = {'vel','torques','muscles','muscle_neurons'};
biomechParams.kin_model = 'real';
biomechParams.comp_blocks = [1,2,3];
biomechParams.pos_origin = [6 30];

td_plot = td(:);
trial_idx = 1:length(td_plot);
bumpTrials = trial_idx(~isnan([td_plot.bumpDir]));

% idx = bumpTrials(200)

idx = 1;


% biomechParams.idx_start = td_plot(idx).idx_bumpTime;
% biomechParams.idx_end = td_plot(idx).idx_goCueTime;
% plot_arm_sim(sim_data(idx),biomechParams);

biomechParams.idx_start = td_plot(idx).idx_goCueTime;
biomechParams.idx_end = td_plot(idx).idx_goCueTime+70;
plot_arm_sim(sim_data(idx),biomechParams);

td_plot(idx).bumpDir
td_plot(idx).tgtDir

cols = linspecer(6);
figure; hold on
subplot(2,1,1); hold on
title('flexors')
for mus = [1,3,5]
    plot(sim_data(idx).muscles(:,mus),'Color',cols(mus,:))
    plot(sim_data(idx).acts(:,mus),'Color',cols(mus,:),'LineWidth',2)
end
% plot([td_plot(idx).idx_bumpTime;td_plot(idx).idx_bumpTime],[0 25],'color','k')
plot([td_plot(idx).idx_goCueTime;td_plot(idx).idx_goCueTime],[0 25],'color','k')

subplot(2,1,2); hold on
for mus = [1,3,5]
    plot(sim_data(idx).muscle_L(:,mus),'Color',cols(mus,:))
end
% plot([td_plot(idx).idx_bumpTime;td_plot(idx).idx_bumpTime],[0.5 1.5],'color','k')
plot([td_plot(idx).idx_goCueTime;td_plot(idx).idx_goCueTime],[0.5 1.5],'color','k')


figure; hold on
subplot(2,1,1); hold on
title('extensors')
for mus = [2,4,6]
    plot(sim_data(idx).muscles(:,mus),'Color',cols(mus,:))
    plot(sim_data(idx).acts(:,mus),'Color',cols(mus,:),'LineWidth',2)
end
% plot([td_plot(idx).idx_bumpTime;td_plot(idx).idx_bumpTime],[0 25],'color','k')
plot([td_plot(idx).idx_goCueTime;td_plot(idx).idx_goCueTime],[0 25],'color','k')

subplot(2,1,2); hold on
for mus = [2,4,6]
    plot(sim_data(idx).muscle_L(:,mus),'Color',cols(mus,:))
end
% plot([td_plot(idx).idx_bumpTime;td_plot(idx).idx_bumpTime],[0.5 1.5],'color','k')
plot([td_plot(idx).idx_goCueTime;td_plot(idx).idx_goCueTime],[0.5 1.5],'color','k')

% figure; hold on
% plot(td_plot(idx).force(:,1),td_plot(idx).force(:,2))

clearvars -except td sim_data biomechParams td_path biomch_path td_name_pre
%% Set up spindle params structure array

% Normalize the muscle activations
if ~exist('sim_data','var')
    load(fullfile(biomch_path,[td_name_pre(1:end-5) '_sim_data_FV_5ms']))
end

for iTrial = 1:length(sim_data)
    sim_data(iTrial).acts = movmean(sim_data(iTrial).acts,10);
end

catActs = cat(1,sim_data.acts);
normCatActs = normalize(catActs,'range',[0.0 1.5]);
actFactors = catActs./normCatActs;
actFactors = actFactors(1,:);

% catLens = cat(1,sim_data.muscle_L);



for iTrial = 1:length(sim_data)
    sim_data(iTrial).normActs = sim_data(iTrial).acts./actFactors;
end

td_spindle = td;
emgNames = biomechParams.muscle_names;



for iTrial = 1:length(td_spindle)
    td_spindle(iTrial).muscle_names = emgNames;
    td_spindle(iTrial).emg_names = emgNames;
    td_spindle(iTrial).musLenRel = sim_data(iTrial).muscle_L;
    td_spindle(iTrial).musVelRel = sim_data(iTrial).muscle_v;
    td_spindle(iTrial).emgNorm = sim_data(iTrial).normActs;
    td_spindle(iTrial).musForces = sim_data(iTrial).muscles;
end

for i = 1:numel(emgNames) 
    spindleParams(i).trialInd = 1;
    spindleParams(i).bufferSize = 0; %This was taken care of in TD procoessing
    spindleParams(i).time_step = 0.005;
    spindleParams(i).gamma = 'agc';

    spindleParams(i).startIdx = {'idx_startTime',0}; %reference index and relative idx
    spindleParams(i).endIdx = {'idx_endTime',0};  
    spindleParams(i).dataStore = 'lean';
    spindleParams(i).emgName = biomechParams.muscle_names{i};
    spindleParams(i).musName = biomechParams.muscle_names{i};

end


% clearvars -except td td_spindle sim_data biomechParams td_path biomch_path td_name_pre spindleParams

%% Run spindle model on muscle kinematics data
tic;
% spindleData = struct('SimOut',struct,'metaParams',struct);
clear spindleData
parfor iTrial = 1:length(td_spindle)
    for mus = 1:numel(spindleParams)
        spindleData(iTrial).SimOut(mus) = getAffPotFromMusState(td_spindle(iTrial),spindleParams(mus));
        spindleData(iTrial).metaParams(mus) = spindleParams(mus);
       
        thisMus = spindleParams(mus).musName;
        disp(['Now Simulating Trial: ' thisMus ' ' num2str(iTrial)])
        
    end
end

toc;

%% Get firing rate from spindle forces

% if ~exist('spindleData','var')
%     load(fullfile(spdl_path,[td_name_pre(1:end-5) '_spindleSym2']));
% end

for iTrial = 1:length(td_spindle)
    for mus = 1:numel(spindleParams)
        dataB = spindleData(iTrial).SimOut(mus).dataB;
        dataC = spindleData(iTrial).SimOut(mus).dataC;
        [r,rs,rd] = sarc2spindle(dataB,dataC,2,1,0.05,0.0,0.0);

        spindleData(iTrial).SimOut(mus).r = r;
    end
end

% close all
trial = 600;

figure; hold on;
plot(cat(2,spindleData(trial).SimOut.r))

figure; hold on
plot(-td_spindle(trial).musVelRel)

figure; hold on
plot(td_spindle(trial).musLenRel)

figure; hold on;
plot(td_spindle(trial).emgNorm)
%% Save spindle data

save_sim = 1;
spdl_path = fullfile(td_path,'spindle');

if save_sim
    save(fullfile(spdl_path,[td_name_pre(1:end-5) '_spindleASymFV']),'spindleData','spindleParams')
end
% clearvars -except td sim_data biomechParams td_path biomch_path td_name_pre



%% Rearrange spindle model output and put into TD
temp = td_spindle;
% temp.spindle = [];

interp_factor = temp(1).bin_size/spindleParams(1).time_step;


for trial = 1:numel(temp)
    for mus = 1:numel(spindleData(1).metaParams)
        if ~mod(interp_factor,1)
            temp(trial).spindle(:,mus) = spindleData(trial).SimOut(mus).r(1:interp_factor:end);
        else
            warning('Spindle data was not generated with easy interpolation factor. Storing complete data...')
            temp(trial).spindle(:,mus) = spindleData(trial).SimOut(mus).r;
        end   
    end
end

trial_data = temp;



%% Save td

tdname = fullfile(td_path,[td_name_pre(1:end-5) '_adaptation_spindleASymFV.mat']);

save(tdname,'trial_data')



