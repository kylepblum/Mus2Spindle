function plotCombinedSignals(trial_data,spindle_data,params)


msd = spindle_data;
tds = trial_data;
numCond = 4;

if strcmpi(params.trialType,'bump')
    conds = [0 90 180 270];
else
    conds = [0 pi/2 pi 3*pi/2]; %Chris's data mixes radians and degrees use
end

Sig1Temp = [];
Sig2Temp = [];
Sig3Temp = [];
FRTemp = [];
combTemp = [];


hfig = figure; hold on;
set(hfig,'units','normalized','outerposition',[0 0.0 0.5 0.4],'PaperSize',[8.5 11],...
    'Renderer','Painters','Color','white');
htitle = sgtitle(params.muscles{params.musIdx},'FontName','Helvetica','FontSize',16,'Interpreter','None');


for i = 1:numCond
    bump_params.bumpDir = conds(i);
    bump_params.targDir = conds(i);
    
    if strcmpi(params.trialType,'bump')
        trialsToPlot = getBumpTrials(tds,bump_params);
    else
        trialsToPlot = getActTrials(tds,bump_params);
    end


    topAxesPosition = [0.1+0.21*(i-1) 0.7 0.18 0.23];
    midAxesPosition = [0.1+0.21*(i-1) 0.4 0.18 0.23];
    botAxesPosition = [0.1+0.21*(i-1) 0.1 0.18 0.23];
    
    htop(i) = subplot(3,4,i); hold on; axis([0.0 1.2 0.0 1.0]);
    set(htop(i),'Position',topAxesPosition,'xticklabel',[],...
        'FontName','Helvetica','FontSize',12)
    text(0.3, 0.4,[num2str(conds(i)) '(n=' num2str(numel(trialsToPlot)) ')']);
    
    hmid(i) = subplot(3,4,i+4); hold on; axis([0.0 1.2 0.0 1.0]);
    set(hmid(i),'Position',midAxesPosition,'xticklabel',[],...
        'FontName','Helvetica','FontSize',12)

    hbot(i) = subplot(3,4,i+8); hold on; axis([0.0 1.2 0.0 1.0]);
    set(hbot(i),'Position',botAxesPosition,...
        'FontName','Helvetica','FontSize',12)
    xlabel('time (s)')

    
    if i > 1
%         set(htop(i),'yticklabel',[])
%         set(hmid(i),'yticklabel',[])
%         set(hbot(i),'yticklabel',[])
    else
        htop(i).YLabel.String = 'Combined Signal (au)';
        htop(i).YLabel.FontName = 'Helvetica';
        hmid(i).YLabel.String = 'Signal 1 (au)';
        hmid(i).YLabel.FontName = 'Helvetica';        
        hbot(i).YLabel.String = 'Signal 2 (au)';
        hbot(i).YLabel.FontName = 'Helvetica';
    end
    
    for trial = 1:numel(trialsToPlot)
        
        thisTrial = trialsToPlot(trial);
        
        if strcmpi(params.trialType,'bump')
            emgIdx = (tds(thisTrial).idx_bumpTime):(tds(thisTrial).idx_bumpTime+100);
            msIdx = 100:2:300; %MSdata is twice the precision and has 100 sample buffer
        else
            emgIdx = (tds(thisTrial).idx_goCueTime):(tds(thisTrial).idx_goCueTime+100);
            msIdx = 100:2:300;
        end

        try
            Sig1Temp(:,end+1) = msd(thisTrial).r(msIdx);
%             Sig1Temp(:,end+1) = tds(thisTrial).musVelRel(emgIdx,params.musIdx);
        catch
            Temp = msd(thisTrial).r(msIdx(1):2:end);
            sizeDiff = numel(msIdx) - numel(Temp);
            Temp(end:end+sizeDiff, end) = NaN;
            Sig1Temp(:,end+1) = Temp;
        end
        Sig2Temp(:,end+1) = tds(thisTrial).emgNorm(emgIdx,params.emgToPlot);
        
%         Sig3Temp(:,end+1) = msd(thisTrial).r(
        
        FRTemp(:,end+1) = tds(thisTrial).cuneate_spikes(emgIdx,4)/100;
        
%         Sig3Temp(:,end+1) = tds(thisTrial).emgNorm(emgIdx,20);
        
        
        combTemp(:,end+1) = ((Sig1Temp(:,end) + Sig2Temp(:,end)));
        
        time = (0:numel(Sig1Temp(:,end))-1)'*0.01;
%         line(time,combTemp(:,end),'Parent',htop(i),'LineWidth',0.5,'Color',[0.7 0.7 0.8])      
        line(time,Sig1Temp(:,end),'Parent',hmid(i),'LineWidth',0.5,'Color',[0.7 0.7 0.8])
        line(time,Sig2Temp(:,end),'Parent',hbot(i),'LineWidth',0.5,'Color',[0.7 0.7 0.8])

        
        
        %         end
    end
    meanSignal1(:,i) = mean(Sig1Temp,2);
    meanSignal1(:,i) =  meanSignal1(:,i);
    meanSignal1(meanSignal1<=0) = 0;
%     meanSignal1(:,i) = zeros(size(meanSignal1(:,i)));
    semSignal1(:,i) = std(Sig1Temp,[],2)./sqrt(size(Sig1Temp,2));
    
    meanSignal2(:,i) = (mean(Sig2Temp,2));
    meanSignal2(:,i) = meanSignal2(:,i)+0.5;
    meanSignal2(meanSignal2<=0) = 0;
    meanSignal2 = meanSignal2.^0.5;
%     meanSignal2(:,i) = zeros(size(meanSignal2(:,i)));
    semSignal2(:,i) = std(Sig2Temp,[],2)./sqrt(size(Sig2Temp,2));
    
%     meanSignal3(:,i) = (mean(Sig3Temp,2));
    
    meanFR(:,i) = mean(FRTemp,2);
    X = [meanSignal1(:,i) meanSignal2(:,i)];
    mdl = fitglm(X, meanFR(:,i),'linear','Distribution', 'poisson');
    
    B = mdl.Coefficients.Estimate;
    B = [-1.5 -1 1]';
    yhat = glmval(B,X,'log');
    
    R = corrcoef(yhat,meanFR(:,i));
    R2 = R(2)^2;
    text(0.5,0.95,['R^2 = ' num2str(R2)],'Parent',htop(i))
    
    meanComb(:,i) = yhat;
    semComb(:,i) = nan*zeros(size(meanComb(:,i)));
%     meanComb(:,i) = mean(combTemp,2);
%     semComb(:,i) = std(combTemp,[],2)./sqrt(size(combTemp,2));

    line(time,meanSignal1(:,i),'lineWidth',5,'Color',[0 0 1],'Parent',hmid(i))
    line(time,meanSignal1(:,i)-semSignal1(:,i),'lineWidth',2,'Color',[0 0 0.5],'Parent',hmid(i))
    line(time,meanSignal1(:,i)+semSignal1(:,i),'lineWidth',2,'Color',[0 0 0.5],'Parent',hmid(i))
   
    line(time,meanSignal2(:,i),'lineWidth',5,'Color',[0 0 1],'Parent',hbot(i))
    line(time,meanSignal2(:,i)-semSignal2(:,i),'lineWidth',2,'Color',[0 0 0.5],'Parent',hbot(i))
    line(time,meanSignal2(:,i)+semSignal2(:,i),'lineWidth',2,'Color',[0 0 0.5],'Parent',hbot(i))

    line(time,meanComb(:,i),'lineWidth',5,'Color',[0 0 1],'Parent',htop(i))
    line(time,meanComb(:,i)-semComb(:,i),'lineWidth',2,'Color',[0 0 0.5],'Parent',htop(i))
    line(time,meanComb(:,i)+semComb(:,i),'lineWidth',2,'Color',[0 0 0.5],'Parent',htop(i))
    line(time,meanFR(:,i),'lineWidth',2,'Color',[1 0 0],'Parent',htop(i));
%     line(time(motorOn(:,1)),0.15*ones(size(time(motorOn(:,1)))),'lineWidth',5,'color',[0.3 0.8 0.3],'Parent',htop(i))
    
    Sig1Temp = [];
    Sig2Temp = [];
    combTemp = [];
    
    
    

end



if params.savefig == 1
    saveas(hfig,[params.savepath params.muscles{params.musIdx} '_' params.trialType params.filetype])
    close(hfig)
end

end