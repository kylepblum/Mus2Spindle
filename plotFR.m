function plotFR(trial_data,params)


td = trial_data;
frToPlot = params.frToPlot;
numCond = 4;
if strcmpi(params.trialType,'bump')
    conds = 0:45:315;
else
%     conds = [0 90 180 270];
    conds = [0 pi/2 pi 3*pi/2]; %Chris's data mixes radians and degrees use
end
FRtemp = [];
FRtemp2 = [];
lenTemp1 = [];
lenTemp2 = [];


hfig = figure; hold on;
set(hfig,'units','normalized','outerposition',[0 0 0.5 0.4],'PaperSize',[5 7],...
    'Renderer','Painters','Color','White');
% htitle = sgtitle(tds(1).emg_names{frToPlot},'FontName','Helvetica','FontSize',16,'Interpreter','None');


for i = 1:numCond
    bump_params.bumpDir = conds(i);
    bump_params.targDir = conds(i);
%     bump_params.target_direction = conds(i);
    if strcmpi(params.trialType,'bump')
        trialsToPlot = getBumpTrials(td,bump_params);
    else
        trialsToPlot = getActTrials(td,bump_params);
    end
    fr_idx = frToPlot;

    topAxesPosition = [0.1+0.21*(i-1) 0.5 0.18 0.35];
    botAxesPosition = [0.1+0.21*(i-1) 0.1 0.18 0.35];
    htop(i) = subplot(2,4,i); hold on; axis([-0.1 1 0 0.1]);
    set(htop(i),'Position',topAxesPosition,'xticklabel',[],...
        'FontName','Helvetica','FontSize',12)
    title([num2str(conds(i)) '(n=' num2str(numel(trialsToPlot)) ')']);
    
    hbot(i) = subplot(2,4,i+4); hold on; axis([-0.1 1 -20 20]);
    set(hbot(i),'Position',botAxesPosition,...
        'FontName','Helvetica','FontSize',12)
    xlabel('time (s)')

    
    if i > 1
        set(htop(i),'yticklabel',[])
        set(hbot(i),'yticklabel',[])
    else
        htop(i).YLabel.String = 'Firing Rate (spikes/s)';
        htop(i).YLabel.FontName = 'Helvetica';
        hbot(i).YLabel.String = 'len (L0)';
        hbot(i).YLabel.FontName = 'Helvetica';
    end
    
    for trial = 1:numel(trialsToPlot)
        
        thisTrial = trialsToPlot(trial);
        if strcmpi(params.trialType,'bump')
            if isfield(td,'idx_bumpTime')
                bumpIdx = (td(thisTrial).idx_bumpTime-100):(td(thisTrial).idx_bumpTime+1000);
            else 
                bumpIdx = (td(thisTrial).idx_endTime-1000):(td(thisTrial).idx_endTime);
            end
                
        else
            bumpIdx = (td(thisTrial).idx_goCueTime+50):td(thisTrial).idx_goCueTime+1000;
        end
        if ~isnan(bumpIdx)
            
            FRsignal = td(thisTrial).(params.spikeArray)(bumpIdx,fr_idx);
%             EMGsignal = smooth(EMGsignal,50);
%             motorOn = tds(thisTrial).motor_control(bumpIdx,:)>100;
%               motorOn = 101:200;
            POSsignalx = td(thisTrial).pos(bumpIdx,1) - td(1).pos(1,1);
            POSsignaly = td(thisTrial).pos(bumpIdx,2) - td(1).pos(1,2);
%             
%             POSsignalMus = tds(thisTrial).musLenRel(bumpIdx,params.musIdx);
            lenTemp1(:,end+1) = POSsignalx;
            lenTemp2(:,end+1) = POSsignaly;
            
%             time = (-100:numel(EMGsignal)-101)'*0.01;
            time = (-10:numel(FRsignal)-11)'*0.001;
            line(time,FRsignal,'Parent',htop(i))
            FRtemp(:,end+1) = FRsignal;
%             line(time,POSsignalx,'Parent',hbot(i))
%             line(time,POSsignaly,'Parent',hbot(i),'Color',[1 0 0])
            line(time,POSsignalx,'Parent',hbot(i),'Color',[0.8 0.8 1])
            line(time,POSsignaly,'Parent',hbot(i),'Color',[1 0.8 0.8])
                        
            
            
        end 
    end
    meanSignal(:,i) = nanmean(FRtemp,2);
    semSignal(:,i) = std(FRtemp,[],2)./sqrt(size(FRtemp,2));
    meanLen1(:,i) = nanmean(lenTemp1,2);
    meanLen2(:,i) = nanmean(lenTemp2,2);
    line(time,meanLen1(:,i),'linewidth',5,'Color',[0 0 1],'Parent',hbot(i))
    line(time,meanLen2(:,i),'linewidth',5,'Color',[1 0 0],'Parent',hbot(i))

    line(time,meanSignal(:,i),'lineWidth',5,'Color',[0 0 1],'Parent',htop(i))
    line(time,meanSignal(:,i)-semSignal(:,i),'lineWidth',2,'Color',[0 0 0.5],'Parent',htop(i))
    line(time,meanSignal(:,i)+semSignal(:,i),'lineWidth',2,'Color',[0 0 0.5],'Parent',htop(i))
%     line(time(motorOn(:)),0.15*ones(size(time(motorOn(:)))),'lineWidth',5,'color',[0.3 0.8 0.3],'Parent',htop(i))
    
    FRtemp = [];
    FRtemp2 = [];
    lenTemp1 = [];
    lenTemp2 = [];
    
    

end
% if params.savefig == 1
%     saveas(hfig,[params.savepath tds(1).musNames{frToPlot} params.filetype])
%     close(hfig)
% end

end