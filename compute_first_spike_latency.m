function [mean_fsl_cellselect,sem_fsl_cellselect,p_comp_fsl_acrosscond_cellselect,p_comp_fsl_acrosscellselect,FSLavg,FSLsig]=compute_first_spike_latency(trigs,DiscreteData,max_first_spike_lat,ppms,CellSelect)

Raster=cell(size(DiscreteData,2),numel(trigs));
first_spike_times=cell(size(Raster));
FSLsig=nan(size(Raster,1),1);
for cellN=1:size(DiscreteData,2)
    for trigtype=1:numel(trigs)
        first_spike_times{cellN,trigtype}=nan(size(trigs{trigtype}{cellN},1),1);
        [Raster{cellN,trigtype}]=get_raster_single(trigs{trigtype}{cellN}, DiscreteData(cellN).Spikes,ppms,0,max_first_spike_lat);
        if ~isempty(Raster{cellN,trigtype})
        [a,ia,~]=unique(Raster{cellN,trigtype}(:,2));
        first_spike_times{cellN,trigtype}(a)=Raster{cellN,trigtype}(ia,1);
        end
    end
    if all(cellfun(@(x) any(~isnan(x)),first_spike_times(cellN,:)))
        [FSLsig(cellN)]=ranksum(first_spike_times{cellN,1},first_spike_times{cellN,2});
    end
end

FSLavg=cellfun(@(x) mean(x,'omitnan'),first_spike_times);

figure('Position',[ 910   176   500   588]);
[mean_fsl_cellselect,sem_fsl_cellselect,p_comp_fsl_acrosscond_cellselect,p_comp_fsl_acrosscellselect]=plot_summary_across_groups(FSLavg,true(size(FSLavg,1),1),{'1','2'},CellSelect,{'VPM';'POm'},'mean','fsl');title('first spike latency')
