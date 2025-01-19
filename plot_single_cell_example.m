function [fig]=plot_single_cell_example(cell_example,trigsP,trigsT,popPSTHs_Puff,popPSTHs_Touch,MeanPtrigPuffState,MeanPTrigAngle,MeanTtrigTouchState,MeanTTrigAngle,TrialPTrigAngle,TrialTTrigAngle,plotNtrials,DiscreteData,ppms,bins,Name,fig)

spikes=DiscreteData(cell_example).Spikes;
trigs_e=[trigsP(cell_example);trigsT(cell_example)];
trig_names={'Puff';'Touch'};
psth_e={popPSTHs_Puff(cell_example,:);popPSTHs_Touch(cell_example,:)};
trace_e=cat(3,cat(1,MeanPtrigPuffState{cell_example},MeanPTrigAngle{cell_example}),cat(1,MeanTtrigTouchState{cell_example},MeanTTrigAngle{cell_example}));
sem_trace_e=cat(3,cat(1,nan(size(MeanPtrigPuffState{cell_example})),std(TrialPTrigAngle{cell_example},0,1)./sqrt(size(TrialPTrigAngle{cell_example},1))),cat(1,nan(size(MeanTtrigTouchState{cell_example})),std(TrialTTrigAngle{cell_example},0,1)./sqrt(size(TrialTTrigAngle{cell_example},1))));
trace_name={'trig','whisker angle'};
pbins=bins(1:end-1)+diff(bins);
psth_e_smoothed=cell(size(psth_e));
spike_raster=cell(size(trigs_e));
for s=1:numel(trigs_e)
    [spike_raster{s}]=get_raster_single(trigs_e{s}, spikes,ppms,-bins(1),bins(end));
    psth_e_smoothed{s}=smoothdata(psth_e{s},2,'gaussian',milliseconds(3),'SamplePoints',milliseconds(pbins));
end
m=max(cat(1,psth_e_smoothed{:}),[],'all');

tt=tiledlayout(fig,5,2);
title(tt,Name)
for s=1:2
    nexttile(1+(s-1),[2 1]);
    scatter(spike_raster{s}(spike_raster{s}(:,2)<=plotNtrials,1),spike_raster{s}(spike_raster{s}(:,2)<=plotNtrials,2),10,'k.');box off;
    xlim(bins([1 end]));ylim([1 plotNtrials]);
    ylabel('trials')
    title(trig_names{s})

    nexttile(5+(s-1));
    yyaxis left
   if ~any(isnan(sem_trace_e(1,:,s)))
        semplot = cat(2,trace_e(1,:,s) + sem_trace_e(1,:,s),fliplr(trace_e(1,:,s) - sem_trace_e(1,:,s)));
        fill(cat(2,pbins,fliplr(pbins)), semplot,[.7 .7 .7]);
        hold on;
    end
     plot(pbins,trace_e(1,:,s)); hold off;
    ylabel(trace_name{1})
    yyaxis right
    if ~any(isnan(sem_trace_e(2,:,s)))
        semplot = cat(2,trace_e(2,:,s) + sem_trace_e(2,:,s),fliplr(trace_e(2,:,s) - sem_trace_e(2,:,s)));
        fill(cat(2,pbins,fliplr(pbins)), semplot,[.7 .7 .7]);
        hold on;
    end
     plot(pbins,trace_e(2,:,s)); hold off;
    ylabel(trace_name{2})
    xlim(pbins([1 end]));

    nexttile(7+(s-1),[2 1]);
    ax(s)=bar(pbins,psth_e_smoothed{s},1,'FaceColor',ones(1,3)/4,'EdgeColor','none');box off;
    xlim(pbins([1 end]));
    ylabel('Rate (Hz)')
    xlabel('Time (ms)')
    ylim([0 1.2*m])
end