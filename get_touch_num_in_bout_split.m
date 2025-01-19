function [mean_sem_p,H_TrigGSBaseResp,TrigRateIncreasePGS,NTrigsGS]=get_touch_num_in_bout_split(CellSelectIn,isiCutoff,minTgroup,trigs,DiscreteData,ppms,params_rate_time_before_trig,params_rate_time_after_trig,params_rate_binsize,params_rate_baselinebinNum)
%
% isiCutoff=500*20;
% minTgroup=10;
% CellSelectIn={vpm_PR_withT;pom_PR_withT};
Bouts=cell(size(trigs));
trig_InBoutInd=Bouts;
groups=[0 1.5 inf];
trig_groupsplit=cell(numel(trigs),numel(groups)-1);
for X=1:size(trigs)
    if ~isempty(trigs{X})
        [Bouts{X}]=returnBursts(trigs{X},isiCutoff);
        for b=1:numel(Bouts{X})
            Bouts{X}{b}(:,2)=1:numel(Bouts{X}{b});
        end
        trig_InBoutInd{X}=cat(1,Bouts{X}{:});
        for g=1:numel(groups)-1
            trig_groupsplit{X,g}=trig_InBoutInd{X}(trig_InBoutInd{X}(:,2)>groups(g) & trig_InBoutInd{X}(:,2)<groups(g+1),1);
        end
    end
end

[H_TrigGS,TrigRateIncreasePGS,NTrigsGS]=get_PSTHs(trig_groupsplit,DiscreteData,ppms,params_rate_time_before_trig,params_rate_time_after_trig,params_rate_binsize,params_rate_baselinebinNum,'left');

H_TrigGSBase=cellfun(@(x) x(1),H_TrigGS);
H_TrigGSResp=cellfun(@(x) x(2),H_TrigGS);
H_TrigGSBaseResp=cat(2,H_TrigGSBase(:,1),H_TrigGSResp);


CellSelect={CellSelectIn{1}(all(NTrigsGS(CellSelectIn{1},:)>=minTgroup,2));CellSelectIn{2}(all(NTrigsGS(CellSelectIn{2},:)>=minTgroup,2))};
CellSelectName={'VPM';'POm'};
tt=tiledlayout(1,2);
title(tt,'N Touch in Burst resp',sprintf('%u ms cutoff, min %u touches each',isiCutoff/ppms,minTgroup))
mean_sem_p=nan(2,3,3);
for ns=1:2
    nt(ns)=nexttile(ns);
    plotx=categorical({'base','1st','2nd+'},{'base','1st','2nd+'},{'base','1st','2nd+'});
    ploty=H_TrigGSBaseResp(CellSelect{ns},:);
    plot(plotx,ploty,'Color','k');hold on
    plot(plotx,mean(ploty,1),'Color','r','LineWidth',3);
    maxRate=max(ploty,[],'all');
    p(1)=signrank(ploty(:,1),ploty(:,2));
    p(2)=signrank(ploty(:,1),ploty(:,3));
    p(3)=signrank(ploty(:,2),ploty(:,3));
    errorbar(1.5,maxRate*1.1,0.5,'horizontal','r')
    text(1.5,maxRate*1.1,['p=' num2str(p(1))],'HorizontalAlignment','center','VerticalAlignment','bottom')
    errorbar(2,maxRate*1.3,1,'horizontal','r')
    text(2,maxRate*1.3,['p=' num2str(p(2))],'HorizontalAlignment','center','VerticalAlignment','bottom')
    errorbar(2.5,maxRate*1.5,0.5,'horizontal','r')
    text(2.5,maxRate*1.5,['p=' num2str(p(3))],'HorizontalAlignment','center','VerticalAlignment','bottom')
    box off
    title(CellSelectName{ns},sprintf('N=%u',numel(CellSelect{ns})))
    ylabel('Rate [Hz]')
    mean_sem_p(ns,:,:)=cat(3,mean(ploty,1),std(ploty,0,1)./sqrt(size(ploty,1)),p);
end
nt(1).YTick=[0:10:100];
nt(2).YTick=[0:10:100];

