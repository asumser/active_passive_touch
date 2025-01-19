function [NumSD,WModRates,bincenters,corrs_with_wparams,PhaseMod]=comp_whisking_parameter_modulation(wbins,WP,WPbasenames,WPnames,min_t,DiscreteData,exclude,include,safety_marginEx,safety_marginIn,ppms,amp_thresh)

FN=fieldnames(DiscreteData);
AllSpikes=squeeze(struct2cell(rmfield(DiscreteData,FN(2:end))));
fprintf('purging spikes with exclude and include parameters\n')
SpikesSelected=select_trigs_by_state(DiscreteData,exclude,include,safety_marginEx*ppms,safety_marginIn*ppms,AllSpikes);
fprintf('determining timepoints with exclude and include parameters\n')
[spont_idx]=get_all_spont_idx(DiscreteData,exclude,safety_marginEx*ppms,safety_marginIn*ppms,include);


phase_fittype =  fittype('mean_r + amp*cos(phase-pref_phase)', 'dependent',{'lam'},'independent',{'phase'},'coefficients',{'mean_r','amp','pref_phase'});
phase_fitoptions=fitoptions(phase_fittype);

PhaseMod.fitresult=cell(size(DiscreteData));
PhaseMod.ftgof=PhaseMod.fitresult;
PhaseMod.mean_r=zeros(size(DiscreteData));
PhaseMod.amp=PhaseMod.mean_r;
PhaseMod.pref_phase=PhaseMod.mean_r;
PhaseMod.MD=PhaseMod.mean_r;
PhaseMod.SNR=PhaseMod.mean_r;
PhaseMod.p=PhaseMod.mean_r;
WModRates=cell(size(wbins));
NumSD=cell(size(wbins));
bincenters=cell(size(wbins));
corrs_with_wparams=nan(size(DiscreteData,2),numel(wbins),2);
for wParamInd=1:numel(WPnames)
    WModRates{wParamInd}=nan(size(DiscreteData,2),numel(wbins{wParamInd})-1);
    NumSD{wParamInd}=nan(size(DiscreteData,2),numel(wbins{wParamInd})-1);

end


%
w=waitbar(0);
for RecInd=1:size(DiscreteData,2)

    if ~isempty(WP{1}{RecInd})
        timelineBehav=linspace(0,DiscreteData(RecInd).LengthTime,numel(WP{1}{RecInd}));
        timelineEphys=linspace(0,DiscreteData(RecInd).LengthTime,DiscreteData(RecInd).LengthInd);
        selectedSpikesIndBehav=interp1(timelineBehav,1:numel(timelineBehav),timelineEphys(SpikesSelected{RecInd}),'nearest');%convert spike indices to behavior timeline indices
        includeStateIndBehav=unique(interp1(timelineBehav,1:numel(timelineBehav),timelineEphys(spont_idx{RecInd}),'nearest'));
        if numel(includeStateIndBehav)>0
            for wParamInd=1:numel(WPnames)
                if strcmp('Phase',WPnames{wParamInd})
                    BehavBaseInd=find(strcmp(WPbasenames,WPnames{wParamInd}));
                    Phase=WP{BehavBaseInd}{RecInd};
                    PhaseSelected=Phase(includeStateIndBehav);
                    Amp=WP{strcmp(WPbasenames,'Amplitude')}{RecInd};
                    NumSD{wParamInd}(RecInd,:)=histcounts(PhaseSelected(Amp(includeStateIndBehav)>amp_thresh),wbins{wParamInd});
                    whiskingSpikes=selectedSpikesIndBehav(Amp(selectedSpikesIndBehav)>amp_thresh);
                    spikesInPhaseBin=histcounts(Phase(whiskingSpikes),wbins{wParamInd});
                    spikesInPhaseBin(NumSD{wParamInd}(RecInd,:)   <min_t/10)=nan;
                    PhaseRate=1000*spikesInPhaseBin./NumSD{wParamInd}(RecInd,:);
                    WModRates{wParamInd}(RecInd,:)=PhaseRate;
                    try
                        pval= circ_kuipertest(Phase(whiskingSpikes), PhaseSelected(Amp(includeStateIndBehav)>amp_thresh), 200, false);
                    catch
                        pval=nan;
                    end
                    bincenters{wParamInd}=wbins{wParamInd}(1:end-1)+mean(diff(wbins{wParamInd}))/2;
                    if sum(~isnan(PhaseRate))>=numel(bincenters{wParamInd})/2
                        phase_fitoptions.StartPoint=[mean(PhaseRate,'omitnan') range(PhaseRate)/2 0];
                        phase_fitoptions.Lower=[0 0 -pi];
                        phase_fitoptions.Upper=[inf inf pi];

                        [PhaseMod.fitresult{RecInd}, PhaseMod.ftgof{RecInd}] = fit( reshape(bincenters{wParamInd}(~isnan(PhaseRate)),[],1), reshape(PhaseRate(~isnan(PhaseRate)),[],1), phase_fittype, phase_fitoptions);
                        PhaseMod.mean_r(RecInd)=PhaseMod.fitresult{RecInd}.mean_r;
                        PhaseMod.amp(RecInd)=PhaseMod.fitresult{RecInd}.amp;
                        PhaseMod.pref_phase(RecInd)=PhaseMod.fitresult{RecInd}.pref_phase;
                        PhaseMod.MD(RecInd)=2*PhaseMod.amp(RecInd)/PhaseMod.mean_r(RecInd);
                        PhaseMod.SNR(RecInd)= PhaseMod.MD(RecInd).*sqrt(PhaseMod.mean_r(RecInd) * .111);
                        PhaseMod.p(RecInd)=pval;
                    end

                else
                    is_abs=strcmp(WPnames{wParamInd}([1 end]),'||');
                    if is_abs
                        BehavBaseInd=find(strcmp(WPbasenames,WPnames{wParamInd}(2:end-1)));
                    else
                        BehavBaseInd=find(strcmp(WPbasenames,WPnames{wParamInd}));
                    end
                    behav=WP{BehavBaseInd}{RecInd};
                    if any(strcmp(WPbasenames{BehavBaseInd},{'Angle','Setpoint'}))
                        behav=behav-median(behav,'omitnan');
                    end
                    if is_abs
                        behav=abs(behav);
                    end
                    timeInBin=histcounts(behav(includeStateIndBehav),wbins{wParamInd});
                    NumSD{wParamInd}(RecInd,:)    =timeInBin;
                    timeInBin(timeInBin<min_t)=nan;
                    [spikesInBin]=histcounts(behav(selectedSpikesIndBehav),wbins{wParamInd});
                    WModRates{wParamInd}(RecInd,:)=1000*spikesInBin./timeInBin;
                    bincenters{wParamInd}=wbins{wParamInd}(1:end-1)+1;
                    temp=cat(1,WModRates{wParamInd}(RecInd,:),bincenters{wParamInd})';
                    temp=temp(~any(isnan(temp),2),:);
                    if size(temp,1)>3
                    [corrs_with_wparams(RecInd,wParamInd,1),corrs_with_wparams(RecInd,wParamInd,2)]=corr(temp(:,1),temp(:,2));
                    end
                end
            end
            waitbar(RecInd/size(DiscreteData,2),w)
        end
    end
end
    close(w)


