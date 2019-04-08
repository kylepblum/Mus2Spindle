function plotSpindleData(spindle_data)

numTrials = numel(spindle_data);
close all
    


for i = 1:numTrials
   
    h(i) = figure;
    if i <= 5
        h(i).Position = h(i).Position + [4*(i-1) 0 0 0];
    elseif i <= 10
        h(i).Position = h(i).Position + [4*(i-6) -8 0 0];
    elseif i <= 15
    end
    subplot(4,1,1); hold on; 
    line(spindle_data(i).dataB.t,spindle_data(i).r)
    set(gca,'xticklabel',[],...
        'fontsize',8,...
        'fontname','Helvetica',...
        'Clipping','off'); 
    axis([0 4 0 0.3]); ylabel('Ia firing rate (au)')
    
    subplot(4,1,2); hold on;
    line(spindle_data(i).dataB.t,spindle_data(i).dataB.cmd_length/1300)
    set(gca,'xticklabel',[],...
        'fontsize',8,...
        'fontname','Helvetica',...
        'Clipping','off'); 
    axis([0 4 0.9 1.1]); ylabel('Muscle Length (L0)')
    
    subplot(4,1,3); hold on;
    line(spindle_data(i).dataB.t,spindle_data(i).delta_cdl/(0.01*1300))
    set(gca,'xticklabel',[],...
        'fontsize',8,...
        'fontname','Helvetica',...
        'Clipping','off'); 
    axis([0 4 -0.2 0.2]); ylabel('Muscle Vel (L0/s)')
    
    subplot(4,1,4); hold on;
    line(spindle_data(i).dataB.t,spindle_data(i).dataC.f_activated)
    set(gca,'fontsize',8,...
        'fontname','Helvetica',...
        'Clipping','off'); 
    axis([0 4 0.0 0.5]); ylabel('Gamma (au)')
    
    align_Ylabels(h(i));
end
end