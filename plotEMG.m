function plotEMG(trial_data,params)


tds = trial_data;
emgToPlot = params.emgToPlot;
numCond = 4;
conds = [0 90 180 270];
EMGtemp = [];
EMGtemp2 = [];
lenTemp = [];


hfig = fig; hold on;
set(hfig,'units','normalized','outerposition',[0 0 0.5 1],'PaperSize',[5 7],...
    'Renderer','Painters');
htitle = sgtitle(tds(1).emg_names{emgToPlot},'FontName','Helvetica','FontSize',16);


for i = 1:numCond
    bump_params.bumpDir = conds(i);
    trialsToPlot = getBumpTrials(tds,bump_params);
    emg_idx = emgToPlot;

    topAxesPosition = [0.1+0.21*(i-1) 0.5 0.18 0.35];
    botAxesPosition = [0.1+0.21*(i-1) 0.1 0.18 0.35];
    htop(i) = subplot(2,4,i); hold on; axis([-0.1 0.1 0 0.4]);
    set(htop(i),'Position',topAxesPosition,'xticklabel',[],...
        'FontName','Helvetica','FontSize',12)
    title([num2str(conds(i)) '(n=' num2str(numel(trialsToPlot)) ')']);

    hbot(i) = subplot(2,4,i+4); hold on; axis([-0.1 0.1 0.95 1.05]);
    set(hbot(i),'Position',botAxesPosition,...
        'FontName','Helvetica','FontSize',12)
    xlabel('time (s)')

    
    if i > 1
        set(htop(i),'yticklabel',[])
        set(hbot(i),'yticklabel',[])
    else
        htop(i).YLabel.String = 'EMG (au)';
        htop(i).YLabel.FontName = 'Helvetica';
        hbot(i).YLabel.String = 'len (L0)';
        hbot(i).YLabel.FontName = 'Helvetica';
    end
    
    for trial = 1:numel(trialsToPlot)
        
        thisTrial = trialsToPlot(trial);
        bumpIdx = (tds(thisTrial).idx_bumpTime-100):(tds(thisTrial).idx_bumpTime+100);
        if ~isnan(bumpIdx)
            
            EMGsignal = tds(thisTrial).emgNorm(bumpIdx,emg_idx);
%             EMGsignal = smooth(EMGsignal,50);
%             motorOn = tds(thisTrial).motor_control(bumpIdx,:)>100;
              motorOn = 101:200;
%             POSsignalx = tds(thisTrial).pos(bumpIdx,1) - tds(1).pos(1,1);
%             POSsignaly = tds(thisTrial).pos(bumpIdx,2) - tds(1).pos(1,2);
%             
            POSsignalMus = tds(thisTrial).musLenRel(bumpIdx,params.musIdx);
            lenTemp(:,end+1) = POSsignalMus;
            
            time = (-100:numel(EMGsignal)-101)'*0.001;
            line(time,EMGsignal,'Parent',htop(i))
            EMGtemp(:,end+1) = EMGsignal;
%             line(time,POSsignalx,'Parent',hbot(i))
%             line(time,POSsignaly,'Parent',hbot(i),'Color',[1 0 0])
            line(time,POSsignalMus,'Parent',hbot(i))
            
            
        end 
    end
    meanSignal(:,i) = mean(EMGtemp,2);
    semSignal(:,i) = std(EMGtemp,[],2)./sqrt(size(EMGtemp,2));
    meanLen(:,i) = nanmean(lenTemp,2);
    line(time,meanLen(:,i),'linewidth',5,'Color',[0 0 1],'Parent',hbot(i))
    line(time,meanSignal(:,i),'lineWidth',5,'Color',[0 0 1],'Parent',htop(i))
    line(time,meanSignal(:,i)-semSignal(:,i),'lineWidth',2,'Color',[0 0 0.5],'Parent',htop(i))
    line(time,meanSignal(:,i)+semSignal(:,i),'lineWidth',2,'Color',[0 0 0.5],'Parent',htop(i))
    line(time(motorOn(:)),0.15*ones(size(time(motorOn(:)))),'lineWidth',5,'color',[0.3 0.8 0.3],'Parent',htop(i))
    
    EMGtemp = [];
    EMGtemp2 = [];
    lenTemp = [];
    
    

end
if params.savefig == 1
    saveas(hfig,[params.savepath tds(1).emg_names{emgToPlot} params.filetype])
    close(hfig)
end

end