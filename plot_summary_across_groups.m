function [mean_across_groups,sem_across_groups,p_comp_conds,p_comp_groups,p_comp_0_groups_conds]=plot_summary_across_groups(Rates,indiv_sig,RateLabels,CellTypes,CellTypeLabels,plotMeanStd,YlabelName)
%%
%figure
p_comp_conds=nan(numel(CellTypeLabels),1);
p_comp_groups=nan(2,1);
p_comp_0_groups_conds=nan(2,2);
mean_across_groups=nan(numel(CellTypeLabels),2);
sem_across_groups=nan(numel(CellTypeLabels),2);
Rates_in=Rates;
if size(Rates_in,3)>1
    Rates=Rates(:,:,2);
end
maxRate=max(Rates(cat(1,CellTypes{:}),:),[],'all','omitnan');
if all(Rates(cat(1,CellTypes{:}),:)>0,'all')
minRate=-.1*maxRate;
else
minRate=min(Rates(cat(1,CellTypes{:}),:),[],'all','omitnan');
end


off=0;
for CT=1:numel(CellTypeLabels)
    selRates=Rates(CellTypes{CT},:);
    selRatesSig=indiv_sig(CellTypes{CT});
    
    [p]=signrank(selRates(~any(isnan(selRates),2),1),selRates(~any(isnan(selRates),2),2));
    p_comp_0_groups_conds(CT,1)=signrank(selRates(~any(isnan(selRates),2),1));
    p_comp_0_groups_conds(CT,2)=signrank(selRates(~any(isnan(selRates),2),2));
    if sum(selRatesSig)>0
    plot(off+[1:2],selRates(selRatesSig,:)','k'); hold on
    end
    if any(~selRatesSig)
    plot(off+[1:2],selRates(~selRatesSig,:)','Color',[.75 .75 .75]); hold on
    end
    switch plotMeanStd
        case 'meanstd'
        errorbar(off+[.85 2.15],mean(selRates,1,'omitnan'),std(selRates,0,1,'omitnan'),'r')
        case 'mean'
            line(off+[.85 2.15],mean(selRates,1,'omitnan'),'Color','r')
    end
    p_comp_conds(CT)=p;
    mean_across_groups(CT,:)=mean(selRates,1,'omitnan');
    sem_across_groups(CT,:)=std(selRates,0,1,'omitnan')./sqrt(sum(~isnan(selRates)));
   
    errorbar(off+1.5,maxRate*1.05,.65,'horizontal','r')
    text(off+1.5,maxRate*1.05,['p=' num2str(p)],'HorizontalAlignment','center','VerticalAlignment','bottom')
    
    text(off+1.5,maxRate*.9,CellTypeLabels{CT},'FontSize',14,'FontWeight','bold','HorizontalAlignment','center','VerticalAlignment','middle')
    
    off=off+2;
    
    %
end
%
off=0;
if numel(CellTypeLabels)>1
    for c=1:2
        [p]=ranksum(Rates_in(CellTypes{1},c,1),Rates_in(CellTypes{2},c,1));
        p_comp_groups(c)=p;
        errorbar(off+1.85,maxRate*1.15+(off)*maxRate*0.05,1,'horizontal','r')
        text(off+1.85,maxRate*1.15+(off)*maxRate*0.05,['p=' num2str(p)],'HorizontalAlignment','center','VerticalAlignment','bottom')
        off=off+1.3;
    end
end

ylim([1.3*minRate maxRate*1.3])
xlim([.75 2*numel(CellTypeLabels)+.25])
set(gca,'Xtick',1:2*numel(CellTypeLabels))
%
set(gca,'Xticklabel',RateLabels(repmat(1:2,1,2)))
box off
ylabel(YlabelName)
%%

