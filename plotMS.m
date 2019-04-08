function plotMS(trial_data,spindle_data,params)


msd = spindle_data;
tds = trial_data;
numCond = 4;

if strcmpi(params.trialType,'bump')
    conds = [0 90 180 270];
else
    conds = [0 pi/2 pi 3*pi/2]; %Chris's data mixes radians and degrees use
end

MStemp = [];
MStemp2 = [];
lenTemp = [];
gamTemp = [];
gamTemp2 = [];


hfig = figure; hold on;
set(hfig,'units','normalized','outerposition',[0 0.0 0.4 1],'PaperSize',[8.5 11],...
    'Renderer','Painters','Color','white');
htitle = sgtitle(params.muscles{params.musIdx},'FontName','Helvetica','FontSize',16,'Interpreter','None');


for i = 1:numCond
    bump_params.bumpDir = conds(i);
    bump_params.targDir = conds(i);
%     trialsToPlot = getBumpTrials(trial_data,bump_params);
    trialsToPlot = getActTrials(trial_data,bump_params);
    
    if strcmpi(params.trialType,'bump')
        trialsToPlot = getBumpTrials(tds,bump_params);
    else
        trialsToPlot = getActTrials(tds,bump_params);
    end

    ms_idx = trialsToPlot;

    topAxesPosition = [0.1+0.21*(i-1) 0.7 0.18 0.23];
    midAxesPosition = [0.1+0.21*(i-1) 0.4 0.18 0.23];
    botAxesPosition = [0.1+0.21*(i-1) 0.1 0.18 0.23];
    htop(i) = subplot(3,4,i); hold on; axis([0.0 1.2 0.0 0.5]);
    set(htop(i),'Position',topAxesPosition,'xticklabel',[],...
        'FontName','Helvetica','FontSize',12)
    text(0.3, 0.4,[num2str(conds(i)) '(n=' num2str(numel(trialsToPlot)) ')']);
    
    hmid(i) = subplot(3,4,i+4); hold on; axis([0.0 1.2 0.1 0.9]);
    set(hmid(i),'Position',midAxesPosition,'xticklabel',[],...
        'FontName','Helvetica','FontSize',12)

    hbot(i) = subplot(3,4,i+8); hold on; axis([0.0 1.2 0.9 1.1]);
    set(hbot(i),'Position',botAxesPosition,...
        'FontName','Helvetica','FontSize',12)
    xlabel('time (s)')

    
    if i > 1
        set(htop(i),'yticklabel',[])
        set(hmid(i),'yticklabel',[])
        set(hbot(i),'yticklabel',[])
    else
        htop(i).YLabel.String = 'Firing Rate (au)';
        htop(i).YLabel.FontName = 'Helvetica';
        hmid(i).YLabel.String = 'Gamma (au)';
        hmid(i).YLabel.FontName = 'Helvetica';        
        hbot(i).YLabel.String = 'len (L0)';
        hbot(i).YLabel.FontName = 'Helvetica';
    end
    
    for trial = 1:numel(trialsToPlot)
        
        thisTrial = trialsToPlot(trial);
        %         bumpIdx = (msd(thisTrial).idx_bumpTime-100):(msd(thisTrial).idx_bumpTime+500);
        %         if ~isnan(bumpIdx)
        
        if strcmpi(params.trialType,'bump')
            timeIdx = 101:221;
        else
            timeIdx = 132:282;
        end
        
        MSsignal = msd(thisTrial).r(timeIdx);
%         EMGsignal = smooth(EMGsignal,50);
%         motorOn = tds(thisTrial).motor_control(:,:)>5;
%         POSsignalx = (tds(thisTrial).pos(1:200,1) - tds(1).pos(1,1))/tds(1).pos(1,1);
% %         POSsignaly = tds(thisTrial).pos(bumpIdx,2) - tds(1).pos(1,2);
        %
        POSsignalMus = msd(thisTrial).dataB.cmd_length(timeIdx)/1300;
%         MSsignal = MSsignal + (POSsignalMus - 1)*5;
        
        lenTemp(:,end+1) = POSsignalMus;
        
        GAMsignalDyn = msd(thisTrial).dataB.f_activated(timeIdx);
        GAMsignalStc = msd(thisTrial).dataC.f_activated(timeIdx);
        gamTemp(:,end+1) = GAMsignalDyn;
        gamTemp2(:,end+1) = GAMsignalStc;
        
        time = (0:numel(MSsignal)-1)'*0.005;
        line(time,MSsignal,'Parent',htop(i),'LineWidth',0.5,'Color',[0.7 0.7 0.8])
        MStemp(:,end+1) = MSsignal;
        
%         line(time,POSsignalx,'Parent',hbot(i))
%         line(time,POSsignaly,'Parent',hbot(i),'Color',[1 0 0])
        line(time,POSsignalMus,'Parent',hbot(i),'LineWidth',0.5,'Color',[0.7 0.7 0.8])
        line(time,GAMsignalDyn,'Parent',hmid(i),'LineWidth',0.5,'Color',[0.7 0.7 0.8])
        line(time,GAMsignalStc,'Parent',hmid(i),'LineWidth',0.5,'Color',[0.8 0.7 0.7])

        
        
        %         end
    end
    meanSignal(:,i) = mean(MStemp,2);
    semSignal(:,i) = std(MStemp,[],2)./sqrt(size(MStemp,2));
    meanLen(:,i) = mean(lenTemp,2);
    meanGam(:,i) = mean(gamTemp,2);
    meanGam2(:,i) = mean(gamTemp2,2);
    line(time,meanLen(:,i),'linewidth',5,'Color',[0 0 1],'Parent',hbot(i))
    line(time,meanSignal(:,i),'lineWidth',5,'Color',[0 0 1],'Parent',htop(i))
    line(time,meanSignal(:,i)-semSignal(:,i),'lineWidth',2,'Color',[0 0 0.5],'Parent',htop(i))
    line(time,meanSignal(:,i)+semSignal(:,i),'lineWidth',2,'Color',[0 0 0.5],'Parent',htop(i))
    line(time,meanGam(:,i),'linewidth',5,'Color',[0 0 1],'Parent',hmid(i))
    line(time,meanGam2(:,i),'linewidth',5,'Color',[1 0 0],'Parent',hmid(i))

%     line(time(motorOn(:,1)),0.15*ones(size(time(motorOn(:,1)))),'lineWidth',5,'color',[0.3 0.8 0.3],'Parent',htop(i))
    
    MStemp = [];
    lenTemp = [];
    gamTemp = [];
    gamTemp2 = [];
    
    
    

end
if params.savefig == 1
    saveas(hfig,[params.savepath params.muscles{params.musIdx} '_' params.trialType params.filetype])
    close(hfig)
end

end