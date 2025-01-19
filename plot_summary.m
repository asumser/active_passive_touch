function [p_summary_vs_baseline,p_summary_across_conditions,mean_across_conditions,sem_across_conditions]=plot_summary(Vals_to_plot,indiv_sig,RateLabels,CellTypes,DeflectionTypeLabels,plotMeanStd,YlabelName,varargin)

if ~isempty(varargin) && strcmp(varargin{1},'dontplot')
    do_plot=false;
else
    do_plot=true;
end
p_summary_vs_baseline=nan(1,numel(DeflectionTypeLabels));
p_summary_across_conditions=nan(2,1);
mean_across_conditions=nan(numel(DeflectionTypeLabels),2);
sem_across_conditions=nan(numel(DeflectionTypeLabels),2);
Vals_to_plot_in=Vals_to_plot;

maxVals_to_plot=max(Vals_to_plot(CellTypes,:,:),[],'all','omitnan');
if all(Vals_to_plot(CellTypes,:,:)>0,'all')
    minVals_to_plot=-.1*maxVals_to_plot;
else
    minVals_to_plot=min(Vals_to_plot(CellTypes,:,:),[],'all','omitnan');
end


off=0;
for DT=1:numel(DeflectionTypeLabels)
    selRates=Vals_to_plot(CellTypes,:,DT);
    selRatesSig=indiv_sig(CellTypes,DT);
    if any(~any(isnan(selRates),2))
        [p]=signrank(selRates(~any(isnan(selRates),2),1),selRates(~any(isnan(selRates),2),2));
    else
        p=1;
    end
    if do_plot
        if any(selRatesSig>0)
            plot(off+[1:2],selRates(selRatesSig>0,:)','k'); hold on
        end
        if any(selRatesSig==0)
            plot(off+[1:2],selRates(~selRatesSig,:)','Color',[.75 .75 .75]); hold on
        end
        if any(selRatesSig<0)
            plot(off+[1:2],selRates(selRatesSig<0,:)','Color',[.35 .35 .35]); hold on
        end
        switch plotMeanStd
            case 'meanstd'
                errorbar(off+[.85 2.15],mean(selRates,1,'omitnan'),std(selRates,0,1,'omitnan'),'r')
            case 'mean'
                line(off+[.85 2.15],mean(selRates,1,'omitnan'),'Color','r')
        end
    end
    p_summary_vs_baseline(DT)=p;
    mean_across_conditions(DT,:,:)=mean(selRates,1,'omitnan');
    sem_across_conditions(DT,:,:)=std(selRates,0,1,'omitnan')./sqrt(sum(~isnan(selRates)));
    if do_plot
        errorbar(off+1.5,maxVals_to_plot*1.05,.65,'horizontal','r')
        text(off+1.5,maxVals_to_plot*1.05,['p=' num2str(p)],'HorizontalAlignment','center','VerticalAlignment','bottom')

        text(off+1.5,maxVals_to_plot*.9,DeflectionTypeLabels{DT},'FontSize',14,'FontWeight','bold','HorizontalAlignment','center','VerticalAlignment','middle')
    end
    off=off+2;

    %
end
%
off=0;
if numel(DeflectionTypeLabels)>1
    for c=1:2
        if any(~any(isnan(Vals_to_plot_in(CellTypes,c,:)),3))
            [p]=signrank(Vals_to_plot_in(CellTypes,c,1),Vals_to_plot_in(CellTypes,c,2));
        else
            p=1;
        end
        p_summary_across_conditions(c)=p;
        if do_plot
        errorbar(off+1.85,maxVals_to_plot*1.15+(off)*maxVals_to_plot*0.05,1,'horizontal','r')
        text(off+1.85,maxVals_to_plot*1.15+(off)*maxVals_to_plot*0.05,['p=' num2str(p)],'HorizontalAlignment','center','VerticalAlignment','bottom')
        end
        off=off+1.3;
    end
end
if do_plot
yl=[minVals_to_plot-5 maxVals_to_plot*1.3];
if ~any(isnan(yl))
    ylim(yl)
    xlim([.75 2*numel(DeflectionTypeLabels)+.25])
end
set(gca,'Xtick',1:2*numel(DeflectionTypeLabels))
%
set(gca,'Xticklabel',RateLabels(repmat(1:2,1,2)))
box off
ylabel(YlabelName)
end
%%

