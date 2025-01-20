%% Settings
%script to run analyses and generate figures for Sumser et al. "Differential
%representation of active and passive touch in mouse somatosensory
%thalamus". Before running for the first time, "whisker_decomposition.m" has
%to be run.

DataPath='D:\active_passive_touch_data\data\';
write_figs=true;
ppms=20;%points per ms (Ephys)
% PSTHs
params_psth_time_before_trig=25;%ms
params_psth_time_after_trig=75;%ms
params_psth_binsize=1;%ms
params_gausssmooth_PSTH_win=5;%ms
params_safety_margin=params_psth_time_after_trig;%ms
%Rates
params_rate_time_before_trig=50;%ms
params_rate_time_after_trig=params_rate_time_before_trig;%ms
params_rate_binsize=params_rate_time_before_trig;%ms
params_rate_baselinebinNum=1;
params_sig=.05; %significance level
params_amp_thresh=3; %minimum amplitude to be assigned as whisking
params_ppms_behav=1; %samples per ms in behavior data
params_puff_triglength_limit=[28 32];%ms filter to select only air puffs with 30 ms length
params_touch_triglength_limit=[0 inf];%ms
params_whisking_state_time_before_trig=500;%ms
params_minTouchEvents=10;%minimum touch events for computing
params_minFillQ=4;%minimum fill quality to include in localized (1=best, 4=soma visible)
params_shuffN=1e4;%number of random spike interval shuffles
params_state_exclude_safety=250;%ms exclude transition times for state computations
params_state_include_safety=25;%ms 
params_max_first_spike_lat=75;%ms maximum first spike latency to be considered
% whisking modulation
params_safety_marginEx_wmod=[100 200];%safety margin to spontaneous ms
params_safety_marginIn_wmod=[];
params_safety_marginEx_wmodL=[5 75];% ms exclude transition times for state computations
params_safety_marginIn_wmodL=[10 0];%ms include light

params_min_t=500;%minimal time in bin for modulation
params_phasebinN=16;
params_phasebins=linspace(-pi,pi,params_phasebinN+1);%phase
params_ampbins=[0:5:30 inf];%amplitude
params_abssetpbins=[0:5:30 inf];%|setpoint|
params_absangbins=[0:5:40 inf];%|angle|
params_setpbins=[-90 -25:5:30 inf];%setpoint
params_angbins=[-90 -30:5:45 inf];%angle

% kinematic splits
params_split_behav_time_before_trig=25;%ms
params_split_behav_after_trig=25;%ms
params_split_binsize=1;%ms
params_smoothwidth=milliseconds(10);%ms
params_splitsmoothmeth='gaussian';
params_split_allR=true;%split by global quantile (false: each cell quantile)
params_splitname='Tertile';
params_splits=linspace(0,100,4);
% load & prepare

load([DataPath 'EphysData.mat'])
load([DataPath 'WhiskerBehavData.mat'])
figdir=[DataPath 'Figs_' datestr(now,'YYmmDD') '\'];
if write_figs;mkdir(figdir);end
set(0,'defaultfigurecolor',[1 1 1], 'DefaultAxesBox', 'off','defaultAxesFontName', 'Arial','defaultTextFontName', 'Arial');
rate_bins=-params_rate_time_before_trig:params_rate_binsize:params_rate_time_after_trig;
psth_bins=-params_psth_time_before_trig:params_psth_binsize:params_psth_time_after_trig;
psth_relativeTime=psth_bins(1:end-1)+diff(psth_bins);

% %% Get events
%Puff
parameter='Puff';exclude={'Pole';'Light';'Exclude';'Grooming'};include=[];
[trigsP,trigsP_length]=get_event_idx(parameter, DiscreteData,params_puff_triglength_limit,ppms,exclude,params_safety_margin,'trigstart',include);
[~,meanAmpPTrigTrial]=get_triggered_avg(trigsP,-params_whisking_state_time_before_trig:0,Amp,ppms/params_ppms_behav,DiscreteData);
trigsPWhisking=cellfun(@(x,y) x(y>params_amp_thresh),trigsP,meanAmpPTrigTrial,'UniformOutput',0);
trigsPNonwhisking=cellfun(@(x,y) x(y<=params_amp_thresh),trigsP,meanAmpPTrigTrial,'UniformOutput',0);
%Puff Light
parameter='Puff';exclude={'Pole';'Exclude';'Grooming'};include={{'Light'};[-5 35]};
[trigsPL,trigsPL_length]=get_event_idx(parameter, DiscreteData,params_puff_triglength_limit,ppms,exclude,params_safety_margin,'trigstart',include);
% Touch
parameter='Touch';exclude={'Puff';'Light';'Exclude';'Grooming'};include=[];
[trigsT,trigsT_length]=get_event_idx(parameter, DiscreteData,params_touch_triglength_limit,ppms,exclude,params_safety_margin,'trigstart',include);

% %% get triggered behavior
[MeanPTrigAngle,~,TrialPTrigAngle]=get_triggered_avg(trigsP,psth_relativeTime,Angle,ppms/params_ppms_behav,DiscreteData);
SEM_PTrigAngle=cellfun(@(x) std(x,0,1)./sqrt(size(x,1)),TrialPTrigAngle,'UniformOutput',false);
[MeanPTrigAmp]=get_triggered_avg(trigsP,psth_relativeTime,Amp,ppms/params_ppms_behav,DiscreteData);
[MeanPWTrigAngle]=get_triggered_avg(trigsPWhisking,psth_relativeTime,Angle,20,DiscreteData);
[MeanPnWTrigAngle]=get_triggered_avg(trigsPNonwhisking,psth_relativeTime,Angle,20,DiscreteData);
[MeanPWTrigAmp]=get_triggered_avg(trigsPWhisking,psth_relativeTime,Amp,20,DiscreteData);
[MeanPnWTrigAmp]=get_triggered_avg(trigsPNonwhisking,psth_relativeTime,Amp,20,DiscreteData);
[MeanPLTrigAngle,~,TrialPLTrigAngle]=get_triggered_avg(trigsPL,psth_relativeTime,Angle,20,DiscreteData);
SEM_PLTrigAngle=cellfun(@(x) std(x,0,1)./sqrt(size(x,1)),TrialPLTrigAngle,'UniformOutput',false);
[MeanPLTrigAmp]=get_triggered_avg(trigsPL,psth_relativeTime,Amp,20,DiscreteData);
[MeanTTrigAngle,~,TrialTTrigAngle]=get_triggered_avg(trigsT,psth_relativeTime,Angle,20,DiscreteData);
SEM_TTrigAngle=cellfun(@(x) std(x,0,1)./sqrt(size(x,1)),TrialTTrigAngle,'UniformOutput',false);
[MeanTTrigAmp]=get_triggered_avg(trigsT,psth_relativeTime,Amp,20,DiscreteData);

% get state vectors 
[PuffFullStateVector]=mk_state_vector(DiscreteData,'Puff');
[MeanPtrigPuffState]=get_triggered_avg(trigsP,psth_relativeTime,PuffFullStateVector,1,DiscreteData);
[MeanPWtrigPuffState]=get_triggered_avg(trigsPWhisking,psth_relativeTime,PuffFullStateVector,1,DiscreteData);
[MeanPnWtrigPuffState]=get_triggered_avg(trigsPNonwhisking,psth_relativeTime,PuffFullStateVector,1,DiscreteData);
[MeanPLtrigPuffState]=get_triggered_avg(trigsPL,psth_relativeTime,PuffFullStateVector,1,DiscreteData);
[MeanPLtrigLightState]=get_triggered_avg(trigsPL,psth_relativeTime,'Light',1,DiscreteData);
[MeanTtrigTouchState]=get_triggered_avg(trigsT,psth_relativeTime,'Touch',1,DiscreteData);

% %% Compute Rates
%puff rates
fprintf('computing puff rates...\n')
[H_Puff,PuffTrigRateIncreaseP,NPTrigs]=get_PSTHs(trigsP,DiscreteData,ppms,params_rate_time_before_trig,params_rate_time_after_trig,params_rate_binsize,params_rate_baselinebinNum,'left');
PuffTrigRateIncreaseP=PuffTrigRateIncreaseP{1}(:,2);
PuffTrigRateSigIncrease=PuffTrigRateIncreaseP<params_sig;
PuffTrigRatePrePost=cat(2,H_Puff{:})';
PuffTrigRateSigIncreaseSign=double(PuffTrigRateSigIncrease).*sign(diff(PuffTrigRatePrePost,1,2));

%touch rates
fprintf('computing touch rates...\n')
[H_Touch,TouchTrigRateIncreaseP,NTTrigs]=get_PSTHs(trigsT,DiscreteData,ppms,params_rate_time_before_trig,params_rate_time_after_trig,params_rate_binsize,params_rate_baselinebinNum,'both');
TouchTrigRateIncreaseP=TouchTrigRateIncreaseP{1}(:,2);
TouchTrigRateSigIncrease=TouchTrigRateIncreaseP<params_sig;
TouchTrigRatePrePost=cat(2,H_Touch{:})';
TouchTrigRateSigIncreaseSign=double(TouchTrigRateSigIncrease).*sign(diff(TouchTrigRatePrePost,1,2));


%qpuff rates
fprintf('computing Qpuff rates...\n')
[H_QPuff,PuffTrigRateIncreaseQP,NQPTrigs]=get_PSTHs(trigsPNonwhisking,DiscreteData,ppms,params_rate_time_before_trig,params_rate_time_after_trig,params_rate_binsize,params_rate_baselinebinNum,'left');
PuffTrigRateIncreaseQP=PuffTrigRateIncreaseQP{1}(:,2);
QPuffTrigRateSigIncrease=PuffTrigRateIncreaseQP<params_sig;
QPuffTrigRatePrePost=cat(2,H_QPuff{:})';

%wpuffrates
fprintf('computing Wpuff rates...\n')
[H_WPuff,PuffTrigRateIncreaseWP,NWPTrigs]=get_PSTHs(trigsPWhisking,DiscreteData,ppms,params_rate_time_before_trig,params_rate_time_after_trig,params_rate_binsize,params_rate_baselinebinNum,'left');
PuffTrigRateIncreaseWP=PuffTrigRateIncreaseWP{1}(:,2);
WPuffTrigRateSigIncrease=PuffTrigRateIncreaseWP<params_sig;
WPuffTrigRatePrePost=cat(2,H_WPuff{:})';

%lpuffrates
fprintf('computing Lpuff rates...\n')
[H_LPuff,PuffTrigRateIncreaseLP,NLPTrigs]=get_PSTHs(trigsPL,DiscreteData,ppms,params_rate_time_before_trig,params_rate_time_after_trig,params_rate_binsize,params_rate_baselinebinNum,'left');
PuffTrigRateIncreaseLP=PuffTrigRateIncreaseLP{1}(:,2);
LPuffTrigRateSigIncrease=PuffTrigRateIncreaseLP<params_sig;
LPuffTrigRatePrePost=cat(2,H_LPuff{:})';

% whisking / quiescence rates
fprintf('computing whisking state rates & stats...\n')
include={'Whisking'};exclude={'Puff';'Touch';'Light';'Exclude';'Grooming'};
[RatesW,RatesWp]=state_rates(DiscreteData,exclude,include, params_state_exclude_safety,params_state_include_safety,ppms,params_shuffN);
RatesWp(RatesWp>.5)=RatesWp(RatesWp>.5)-1;
fprintf('computing nonwhisking state rates & stats...\n')
include={'NoWhisking'};exclude={'Puff';'Touch';'Light';'Exclude';'Grooming'};
[RatesQ,RatesQp]=state_rates(DiscreteData,exclude,include, params_state_exclude_safety,params_state_include_safety,ppms,params_shuffN);
RatesQp(RatesQp>.5)=RatesQp(RatesQp>.5)-1;

% whisking light rates
fprintf('computing whisking light rates...\n')
include={'Light'};exclude={'NoWhisking';'Puff';'Touch';'Exclude';'Grooming'};
state_exclude_safetyL=200;
[RatesWL,RatesWLp]=state_rates(DiscreteData,exclude,include, state_exclude_safetyL,params_state_include_safety,ppms,params_shuffN);
RatesWLp(RatesWLp>.5)=RatesWLp(RatesWLp>.5)-1;

% quiescence light rates

include={'Light'};exclude={'Whisking';'Puff';'Touch';'Exclude';'Grooming'};
[RatesQL,RatesQLp]=state_rates(DiscreteData,exclude,include, state_exclude_safetyL,params_state_include_safety,ppms,params_shuffN);
RatesQLp(RatesQLp>.5)=RatesQLp(RatesQLp>.5)-1;

% whisking parameter modulation
fprintf('computing whisking in air modulation...\n')
exclude={'Pole';'Light';'Exclude';'Grooming';'Puff'};
include=[];
wbins={params_ampbins;params_setpbins;params_abssetpbins;params_angbins;params_absangbins;params_phasebins};
WPnames={'Amplitude';'Setpoint';'|Setpoint|';'Angle';'|Angle|';'Phase'};
WP={Angle,Amp,Setp,Phase};
WPbasenames={'Angle','Amplitude','Setpoint','Phase'};
[NumSD,WModRates,WModBincenters,corrs_with_wparams,PhaseMod]=comp_whisking_parameter_modulation(wbins,WP,WPbasenames,WPnames,params_min_t,DiscreteData,exclude,include,params_safety_marginEx_wmod,params_safety_marginIn_wmod,ppms,params_amp_thresh);
corrwithwparamssig=corrs_with_wparams(:,:,2)<params_sig;

% whisking light parameter modulation
fprintf('computing whisking in air modulation with light...\n')
exclude={'Pole';'Exclude';'Grooming';'Puff'};
include={'Light'};
wbins={params_ampbins;params_setpbins;params_abssetpbins;params_angbins;params_absangbins;params_phasebins};
WPnames={'Amplitude';'Setpoint';'|Setpoint|';'Angle';'|Angle|';'Phase'};
WP={Angle,Amp,Setp,Phase};
WPbasenames={'Angle','Amplitude','Setpoint','Phase'};
[NumSD_L,WModRates_L,WModBincenters_L,corrs_with_wparams_L,PhaseMod_L]=comp_whisking_parameter_modulation(wbins,WP,WPbasenames,WPnames,params_min_t,DiscreteData,exclude,include,params_safety_marginEx_wmodL,params_safety_marginIn_wmodL,ppms,params_amp_thresh);
corrwithwparamssig_L=corrs_with_wparams_L(:,:,2)<params_sig;
pref_phase_nLL=cat(2,PhaseMod.pref_phase',PhaseMod_L.pref_phase');
Phase_SNR_nLL=cat(2,PhaseMod.SNR',PhaseMod_L.SNR');
PhaseMod_p_nLL=cat(2,PhaseMod.p',PhaseMod_L.p');





% response types
is_good_unit = RecordingInfo.FillQuality<=params_minFillQ & RecordingInfo.RecordingType=='juxta';
has_Light=~cellfun(@isempty,MeanPLtrigPuffState);
has_Touch=NTTrigs>=params_minTouchEvents;
is_vpm_of_all=RecordingInfo.DefinedNucleus=='VPM';
is_pom_of_all=RecordingInfo.DefinedNucleus=='POm';
has_PuffResponse=PuffTrigRateSigIncrease;
has_TouchResponse=TouchTrigRateSigIncrease;
any_included=find(is_good_unit & (is_pom_of_all | is_vpm_of_all));
vpm=find(is_good_unit & is_vpm_of_all);
vpm_PR_withT=find(is_good_unit & has_PuffResponse & has_Touch & is_vpm_of_all);
vpm_PR_anyT=find(is_good_unit & has_PuffResponse & is_vpm_of_all);
vpm_PR_anyT_withL=find(is_good_unit & has_PuffResponse & has_Light & is_vpm_of_all);
vpm_nPR_withT=find(is_good_unit & ~has_PuffResponse & has_Touch & is_vpm_of_all);
pom=find(is_good_unit & is_pom_of_all);
pom_PR_withT=find(is_good_unit & has_PuffResponse & has_Touch & is_pom_of_all);
pom_PR_anyT=find(is_good_unit & has_PuffResponse & is_pom_of_all);
pom_PR_anyT_withL=find(is_good_unit & has_PuffResponse & has_Light & is_pom_of_all);
pom_nPR_withT=find(is_good_unit & ~has_PuffResponse & has_Touch & is_pom_of_all);


fprintf('test interval: %u ms; cell numbers: \n vpm puff responsive and active touches:%u, vpm puff responsive overall:%u , vpm puff responsive and cortical inactivation:%u \n pom puff responsive and active touches:%u, pom puff responsive overall:%u , pom puff responsive and cortical inactivation:%u \n',params_rate_time_after_trig,numel(vpm_PR_withT),numel(vpm_PR_anyT),numel(vpm_PR_anyT_withL),numel(pom_PR_withT),numel(pom_PR_anyT),numel(pom_PR_anyT_withL))

% group stats
fprintf('computing group statistics...\n')
% vpm vs pom touch&puff
p_vpmpom_PT_baseline_resp=nan(2,2);
for c=1:2
p_vpmpom_PT_baseline_resp(1,c)=ranksum(PuffTrigRatePrePost(vpm_PR_withT,c),PuffTrigRatePrePost(pom_PR_withT,c));
p_vpmpom_PT_baseline_resp(2,c)=ranksum(TouchTrigRatePrePost(vpm_PR_withT,c),TouchTrigRatePrePost(pom_PR_withT,c));
end


% vpm vs pom qpuff&wpuff
p_vpmpom_QPWP_baseline_resp=nan(2,2);
for c=1:2
p_vpmpom_QPWP_baseline_resp(1,c)=ranksum(QPuffTrigRatePrePost(vpm_PR_withT,c),QPuffTrigRatePrePost(pom_PR_withT,c));
p_vpmpom_QPWP_baseline_resp(2,c)=ranksum(WPuffTrigRatePrePost(vpm_PR_withT,c),WPuffTrigRatePrePost(pom_PR_withT,c));
end

% vpm vs pom light whisking
p_vpmpom_WLWnL_baseline_resp=nan(2,2);
HR_nL=cat(2,RatesQ,RatesW);
HR_L=cat(2,RatesQL,RatesWL);
for c=1:2
p_vpmpom_WLWnL_baseline_resp(1,c)=ranksum(HR_nL(vpm_PR_anyT_withL,c),HR_nL(pom_PR_anyT_withL,c));
p_vpmpom_WLWnL_baseline_resp(2,c)=ranksum(HR_L(vpm_PR_anyT_withL,c),HR_L(pom_PR_anyT_withL,c));
end

p_vpmpom_PLPnL_baseline_resp=nan(2,2);
for c=1:2
p_vpmpom_PLPnL_baseline_resp(1,c)=ranksum(PuffTrigRatePrePost(vpm_PR_anyT_withL,c,1),PuffTrigRatePrePost(pom_PR_anyT_withL,c,1));
p_vpmpom_PLPnL_baseline_resp(2,c)=ranksum(LPuffTrigRatePrePost(vpm_PR_anyT_withL,c,1),LPuffTrigRatePrePost(pom_PR_anyT_withL,c,1));
end


% compute first spike latencies
fprintf('computing first spike latencies...\n')
CellSelect={vpm_PR_withT;pom_PR_withT};
trigs={trigsP;trigsT};
[mean_fsl_vpmpom_PT,sem_fsl_vpmpom_PT,p_comp_fsl_acrosscond_vpmpom_PT,p_comp_fsl_acrossvpmpom_PT]=compute_first_spike_latency(trigs,DiscreteData,params_max_first_spike_lat,ppms,CellSelect);

CellSelect={vpm_PR_withT;pom_PR_withT};
trigs={trigsPNonwhisking;trigsPWhisking};
[mean_fsl_vpmpom_QPWP,sem_fsl_vpmpom_QPWP,p_comp_fsl_acrosscond_vpmpom_QPWP,p_comp_fsl_acrossvpmpom_QPWP]=compute_first_spike_latency(trigs,DiscreteData,params_max_first_spike_lat,ppms,CellSelect);

CellSelect={vpm_PR_anyT_withL;pom_PR_anyT_withL};
trigs={trigsP;trigsPL};
[mean_fsl_vpmpom_PL,sem_fsl_vpmpom_PL,p_comp_fsl_acrosscond_vpmpom_PL,p_comp_fsl_acrossvpmpom_PL]=compute_first_spike_latency(trigs,DiscreteData,params_max_first_spike_lat,ppms,CellSelect);
%% export data for clustering
Cell_numbers=[vpm_PR_withT;pom_PR_withT];
Puff_PSTH=P_PSTH_instrate(Cell_numbers,:);
Touch_PSTH=T_PSTH_instrate(Cell_numbers,:);
Bins_PSTH=-params_psth_time_before_trig:params_psth_binsize:params_psth_time_after_trig;
Bins_PSTH=Bins_PSTH(1:end-1)+params_psth_binsize/2;

is_vpm=false(size(Cell_numbers));is_vpm(1:numel(vpm_PR_withT))=true;
is_pom=false(size(Cell_numbers));is_pom(numel(vpmr)+1:end)=true;
save('vpm_pom_cluster_data.mat','Puff_PSTH','Touch_PSTH',...
    'Bins_PSTH','Cell_numbers','is_vpm','is_pom')

%% Figure 2
Fig_no='Fig2_';

%puff psths
contPlot={MeanPTrigAngle,MeanPtrigPuffState,MeanPTrigAmp};
Name='Puff';
figure('Position',[38 260 1413 525],'Color','White');
[popPSTHs_Puff]=plot_population_psth(trigsP,DiscreteData,ppms,params_psth_time_before_trig,params_psth_time_after_trig,params_psth_binsize,params_gausssmooth_PSTH_win,Name,{vpm_PR_withT;pom_PR_withT},{'VPM';'POm'},'latency',-inf,contPlot,{'meanbefore','scaleminmax','mono'},'baseline corr. whisker angle');
if write_figs;exportgraphics(gcf,[figdir Fig_no 'Pop_' parameter '_PSTH_' datestr(now,'yymmdd') '.pdf'],'BackgroundColor','none');end

%touch psths
Name='Touch';
contPlot={MeanTTrigAngle,MeanTtrigTouchState,MeanTTrigAmp};
figure('Position',[38 260 1413 525],'Color','White');
[popPSTHs_Touch]=plot_population_psth(trigsT,DiscreteData,ppms,params_psth_time_before_trig,params_psth_time_after_trig,params_psth_binsize,params_gausssmooth_PSTH_win,Name,{vpm_PR_withT;pom_PR_withT},{'VPM';'POm'},'latency',-inf,contPlot,{'meanbefore','scaleminmax','mono'},'baseline corr. whisker angle');
if write_figs;exportgraphics(gcf,[figdir Fig_no 'Pop_' parameter '_PSTH_' datestr(now,'yymmdd') '.pdf'],'BackgroundColor','none','ContentType','vector');end

% summary puff / touch
figure('Position',[265         633        1220         587])
tt=tiledlayout(1,3);
title(tt,'Puff/Touch Response')
nexttile(1);
[p_PT_vs_baseline{1},p_P_vs_T{1},mean_PT{1},sem_PT{1}]=plot_summary(cat(3,PuffTrigRatePrePost,TouchTrigRatePrePost),cat(2,PuffTrigRateSigIncrease,TouchTrigRateSigIncreaseSign),{'base';'resp';'base';'resp'},vpm_PR_withT,{'Puff','Touch'},'mean','Rate [Hz]');
title('VPM')
nexttile(2);
[p_PT_vs_baseline{2},p_P_vs_T{2},mean_PT{2},sem_PT{2}]=plot_summary(cat(3,PuffTrigRatePrePost,TouchTrigRatePrePost),cat(2,PuffTrigRateSigIncrease,TouchTrigRateSigIncreaseSign),{'base';'resp';'base';'resp'},pom_PR_withT,{'Puff','Touch'},'mean','Rate [Hz]');
title('POm')
if write_figs;exportgraphics(gcf,[figdir Fig_no 'RatesPT_POm_RespondingCellsWithTouch.pdf'], 'BackgroundColor','none','ContentType','vector');end
nexttile(3);
[mean_mod_vpmpom_PT,sem_mod_vpmpom_PT,p_comp_mod_acrossPT_vpmpom,p_comp_mod_acrossvpmpom_PT,p_comp_mod_vs_0_vpmpom_PT]=plot_summary_across_groups(cat(2,comp_modulation(PuffTrigRatePrePost),comp_modulation(TouchTrigRatePrePost)),true(size(PuffTrigRateSigIncrease)),{'Puff','Touch'},{vpm_PR_withT;pom_PR_withT},{'VPM';'POm'},'mean','modulation');title('Modulation')
if write_figs;exportgraphics (gcf,[figdir Fig_no 'Touch_Puff_Summaries.pdf'], 'BackgroundColor','none','ContentType','vector');end

celldist=cat(1, histcounts(PuffTrigRateSigIncreaseSign(vpm_PR_withT),[-1.5:1.5],'Normalization','probability'),...
                histcounts(PuffTrigRateSigIncreaseSign(pom_PR_withT),[-1.5:1.5],'Normalization','probability'),...
                histcounts(TouchTrigRateSigIncreaseSign(vpm_PR_withT),[-1.5:1.5],'Normalization','probability'),...
                histcounts(TouchTrigRateSigIncreaseSign(pom_PR_withT),[-1.5:1.5],'Normalization','probability'));
%
celldistnum=cat(1, histcounts(PuffTrigRateSigIncreaseSign(vpm_PR_withT),[-1.5:1.5]),...
                histcounts(PuffTrigRateSigIncreaseSign(pom_PR_withT),[-1.5:1.5]),...
                histcounts(TouchTrigRateSigIncreaseSign(vpm_PR_withT),[-1.5:1.5]),...
                histcounts(TouchTrigRateSigIncreaseSign(pom_PR_withT),[-1.5:1.5]));
%
figure('Position',[265         922        1224         298])
distnames={'VPM Puff','POm Puff','VPM Touch','POm Touch'};
modtypes={'neg.','nonmod.','pos'};
tiledlayout(1,size(celldist,1))
for n=1:size(celldist,1)
        nexttile;
        pie(celldist(n,:),modtypes)
        title(distnames{n})
end

if write_figs;exportgraphics (gcf,[figdir Fig_no 'Touch_Puff_mod_sign_pie.pdf'], 'BackgroundColor','none','ContentType','vector');end

%% Figure 1 
Fig_no='Fig1_';
plotNtrials=100;
fig1v=figure('Position',[10         245         574         412],'Color','White');
Name='VPM example';
[fig1v]=plot_single_cell_example(vpm_PR_withT(1),trigsP,trigsT,popPSTHs_Puff,popPSTHs_Touch,MeanPtrigPuffState,MeanPTrigAngle,MeanTtrigTouchState,MeanTTrigAngle,TrialPTrigAngle,TrialTTrigAngle,plotNtrials,DiscreteData,ppms,psth_bins,Name,fig1v);
exportgraphics(fig1v,[figdir Fig_no 'example_VPM_PuffTouch.pdf'],'BackgroundColor','none','ContentType','vector')


fig1p=figure('Position',[10         245         574         412],'Color','White');
Name='POm example';
[fig1p]=plot_single_cell_example(pom_PR_withT(9),trigsP,trigsT,popPSTHs_Puff,popPSTHs_Touch,MeanPtrigPuffState,MeanPTrigAngle,MeanTtrigTouchState,MeanTTrigAngle,TrialPTrigAngle,TrialTTrigAngle,plotNtrials,DiscreteData,ppms,psth_bins,Name,fig1p);
exportgraphics(fig1p,[figdir Fig_no 'example_POm_PuffTouch.pdf'],'BackgroundColor','none','ContentType','vector')


fig1c=figure;
coords=RecordingInfo{:,'CellAtlasCoords'};
scatter3(coords(pom,1),-coords(pom,2),-coords(pom,3),5,'red','filled');hold on;
scatter3(coords(pom_PR_withT,1),-coords(pom_PR_withT,2),-coords(pom_PR_withT,3),25,'red','filled');hold on;
scatter3(coords(vpm,1),-coords(vpm,2),-coords(vpm,3),5,'blue','filled');hold on;
scatter3(coords(vpm_PR_withT,1),-coords(vpm_PR_withT,2),-coords(vpm_PR_withT,3),25,'blue','filled');hold on;
sp=100*round(min(coords([vpm;pom],1:3).*[1 -1 -1])/100)+[-100 100 100];
quiver3(sp(1),sp(2),sp(3),200,0,0,'off','color','k');
quiver3(sp(1),sp(2),sp(3),0,0,200,'off','color','k');
quiver3(sp(1),sp(2),sp(3),0,200,0,'off','color','k');
axis equal
view(gca,23.5702,16.4257)
exportgraphics(fig1c,[figdir Fig_no 'cell_locations_all.pdf'],'BackgroundColor','none','ContentType','vector')


%% Figure 3
Fig_no='Fig3_';
figure('Position',[ 910   176   500   588])
[mean_rate_vpmpom_QW,sem_rate_vpmpom_QW,p_comp_rate_acrossQW_vpmpom,p_comp_rate_acrossvpmpom_QW]=plot_summary_across_groups(...
    cat(2,RatesQ,RatesW),abs(RatesWp)<.05,{'Q','W'},{vpm_PR_anyT;pom_PR_anyT},{'VPM';'POm'},'mean','Rate (Hz)');title('Whisking Rate')
if write_figs;exportgraphics (gcf,[figdir Fig_no 'WhiskingRates.pdf'], 'BackgroundColor','none','ContentType','vector');end

figure('Position',[ 910   176   500   588])
[mean_mod_vpmpom_QW,sem_mod_vpmpom_QW,~,p_comp_mod_acrossvpmpom_QW,p_comp_mod_vs_0_vpmpom_QW]=plot_summary_across_groups(...
    repmat(comp_modulation(RatesQ,RatesW),1,2),abs(RatesWp)<.05,{'spont Whisk','spont Whisk'},{vpm_PR_anyT;pom_PR_anyT},{'VPM';'POm'},'mean','mod');title('Whisking Mod','notsaved')


ampX=find(strcmp(WPnames,'Amplitude'));
ampcorr=corrs_with_wparams(:,ampX,1);
ampcorrsig=corrwithwparamssig(:,ampX);
cell_select={vpm_PR_anyT;pom_PR_anyT};
cell_select_names={'VPM';'POm'};
cell_colors={'r';'b'};
figure('Position',[148         198        1593         688])
tt= tiledlayout(1,2);
title(tt,'mean whisking in air correlations amp')

for c=1:numel(cell_select)
    nexttile(c);

    plot(WModBincenters{ampX},WModRates{ampX}(cell_select{c}(~ampcorrsig(cell_select{c})),:),'Color',[.75 .75 .75]);hold on
    plot(WModBincenters{ampX},WModRates{ampX}(cell_select{c}(ampcorrsig(cell_select{c})),:),'Color',[.25 .25 .25]);hold on
    mA=mean(WModRates{ampX}(cell_select{c},:),1,'omitnan');
    [cox, coxp]=corr(WModBincenters{ampX}',mA');
    plot(WModBincenters{ampX},mA,cell_colors{c});hold off;box off
    ylabel('mean Rate (Hz)');
    xlabel('Amplitude (°)');
    title(cell_select_names{c},sprintf('Amp: c=%.2f, p=%.3g',cox,coxp ))
    xlim(WModBincenters{ampX}([1 end]));
end

if write_figs;exportgraphics (gcf,[figdir Fig_no 'response_relative_to_whisking_in_air_amp.pdf'], 'BackgroundColor','none','ContentType','vector');end


figure('Position',[213         667        1509         563],'Color','white')
tt=tiledlayout(4,1);
phaseX=find(strcmp(WPnames,'Phase'));
cell_select=vpm_PR_anyT(6:7);
title(tt,'phase examples')
plotx=linspace(-pi,pi,100);
for r=1:numel(cell_select)
nexttile
bar(WModBincenters{phaseX},WModRates{phaseX}(cell_select(r),:),1,'FaceColor',[.7 .7 .7],'EdgeColor','k');hold on;
plot(plotx,PhaseMod.fitresult{cell_select(r)}(plotx),'Color','b','LineWidth',3)
legend off
xlim([-pi pi])
box off;
title('VPM')
end
cell_select=pom_PR_anyT([2 8]);
for r=1:numel(cell_select)
nexttile
bar(WModBincenters{phaseX},WModRates{phaseX}(cell_select(r),:),1,'FaceColor',[.7 .7 .7],'EdgeColor','k');hold on;
plot(plotx,PhaseMod.fitresult{cell_select(r)}(plotx),'Color','r','LineWidth',3)
legend off
xlim([-pi pi])
box off;
title('POm')
end
if write_figs;exportgraphics (gcf,[figdir Fig_no 'phase_lock_example_vpmpom.pdf'], 'BackgroundColor','none','ContentType','vector');end

cell_select={vpm_PR_anyT;pom_PR_anyT};
figure('Position',[   680         558        1060         420],'Color','white')
tiledlayout(2,3);
nexttile(1,[2 2]);
polarscatter(PhaseMod.pref_phase(cell_select{1}(PhaseMod.p(cell_select{1})>=.05)),PhaseMod.SNR(cell_select{1}(PhaseMod.p(cell_select{1})>=.05)),75,'m.');hold on
polarscatter(PhaseMod.pref_phase(cell_select{1}(PhaseMod.p(cell_select{1})<.05)),PhaseMod.SNR(cell_select{1}(PhaseMod.p(cell_select{1})<.05)),75,'r.');hold on
polarscatter(PhaseMod.pref_phase(cell_select{2}(PhaseMod.p(cell_select{2})>=.05)),PhaseMod.SNR(cell_select{2}(PhaseMod.p(cell_select{2})>=.05)),75,'c.');hold on
polarscatter(PhaseMod.pref_phase(cell_select{2}(PhaseMod.p(cell_select{2})<.05)),PhaseMod.SNR(cell_select{2}(PhaseMod.p(cell_select{2})<.05)),75,'b.');hold off
legend('vpm nonsig.','vpm sig','pom nonsig.','pom sig.')
title('phase SNR','preferred phase');
nexttile(3);
histogram(PhaseMod.SNR(cell_select{1}),0:.25:max(PhaseMod.SNR(cat(1,cell_select{1},cell_select{2})))+.2,'FaceColor','r');hold on;box off
title('phase SNR VPM');
nexttile(6);
histogram(PhaseMod.SNR(cell_select{2}),0:.25:max(PhaseMod.SNR(cat(1,cell_select{1},cell_select{2})))+.2,'FaceColor','b');hold on;box off
title('phase SNR POm');
if write_figs;exportgraphics (gcf,[figdir Fig_no 'phase_mod_whisking_in_air.pdf'], 'BackgroundColor','none','ContentType','vector');end

%
%nonwhisking puff psths
contPlot={MeanPnWTrigAngle,MeanPnWtrigPuffState,MeanPnWTrigAmp};
Name='QPuff';
figure('Position',[38 260 1413 525],'Color','White');
[popPSTHs_QPuff,ax]=plot_population_psth(trigsPNonwhisking,DiscreteData,ppms,params_psth_time_before_trig,params_psth_time_after_trig,params_psth_binsize,params_gausssmooth_PSTH_win,Name,{vpm_PR_withT;pom_PR_withT},{'VPM';'POm'},'latency',-inf,contPlot,{'meanbefore','scaleminmax','mono'},'baseline corr. whisker angle');
ax(2,1).YLim=[0 200];
ax(2,2).YLim=[0 50];
if write_figs;exportgraphics(gcf,[figdir Fig_no 'Pop_QPuff_PSTH_' datestr(now,'yymmdd') '.pdf'],'BackgroundColor','none','ContentType','vector');end

%whisking puff psths
Name='WPuff';
contPlot={MeanPWTrigAngle,MeanPWtrigPuffState,MeanPWTrigAmp};
figure('Position',[38 260 1413 525],'Color','White');
[popPSTHs_WPuff,ax2]=plot_population_psth(trigsPWhisking,DiscreteData,ppms,params_psth_time_before_trig,params_psth_time_after_trig,params_psth_binsize,params_gausssmooth_PSTH_win,Name,{vpm_PR_withT;pom_PR_withT},{'VPM';'POm'},'latency',-inf,contPlot,{'meanbefore','scaleminmax','mono'},'baseline corr. whisker angle');
ax2(2,1).YLim=[0 200];
ax2(2,2).YLim=[0 50];

if write_figs;exportgraphics(gcf,[figdir Fig_no 'Pop_WPuff_PSTH_' datestr(now,'yymmdd') '.pdf'],'BackgroundColor','none','ContentType','vector');end
%
% summary Qpuff / Wpuff
figure('Position',[265         633        1220         587])
tt=tiledlayout(1,4);
title(tt,'QPuff/WPuff Response')
wpax(1)=nexttile(1);
[p_QPWP_vs_baseline{1},p_QP_vs_WP{1},mean_QPWP{1},sem_QPWP{1}]=plot_summary(cat(3,QPuffTrigRatePrePost,WPuffTrigRatePrePost),cat(2,QPuffTrigRateSigIncrease,WPuffTrigRateSigIncrease),{'base';'resp';'base';'resp'},vpm_PR_withT,{'QPuff','WPuff'},'mean','Rate [Hz]');
title('VPM')
wpax(2)=nexttile(2);
[p_QPWP_vs_baseline{2},p_QP_vs_WP{2},mean_QPWP{2},sem_QPWP{2}]=plot_summary(cat(3,QPuffTrigRatePrePost,WPuffTrigRatePrePost),cat(2,QPuffTrigRateSigIncrease,WPuffTrigRateSigIncrease),{'base';'resp';'base';'resp'},pom_PR_withT,{'QPuff','WPuff'},'mean','Rate [Hz]');
title('POm')
wpax(1).YLim=[0 200];
wpax(1).YTick=[0:20:180];
modax(1)=nexttile(3);
[mean_mod_vpmpom_QPWP,sem_mod_vpmpom_QPWP,p_comp_mod_acrossQPWP_vpmpom,p_comp_mod_acrossvpmpom_QPWP,p_comp_mod_vs_0_vpmpom_QPWP]=plot_summary_across_groups(cat(2,comp_modulation(QPuffTrigRatePrePost),comp_modulation(WPuffTrigRatePrePost)),true(size(WPuffTrigRatePrePost,1),1),{'QPuff','WPuff'},{vpm_PR_withT;pom_PR_withT},{'VPM';'POm'},'mean','modulation');title('Modulation')
modax(2)=nexttile(4);
wpax(2).YLim=[0 120];
wpax(2).YTick=[0:10:100];
[mean_mod_vpmpom_WPT,sem_mod_vpmpom_WPT,p_comp_mod_acrossWPT_vpmpom,p_comp_mod_acrossvpmpom_WPT,p_comp_mod_vs_0_vpmpom_WPT]=plot_summary_across_groups(cat(2,comp_modulation(WPuffTrigRatePrePost),comp_modulation(TouchTrigRatePrePost)),true(size(WPuffTrigRatePrePost,1),1),{'WPuff','Touch'},{vpm_PR_withT;pom_PR_withT},{'VPM';'POm'},'mean','modulation');title('Modulation')

for mx=1:2
    modax(mx).YLim=[-.4 1.2];
    modax(mx).YTick=[-.4:.2:1.2];
end
%
if write_figs;exportgraphics (gcf,[figdir Fig_no 'QWPuff_Summaries.pdf'], 'BackgroundColor','none','ContentType','vector');end

%% Figure 4 & Supplementary Figure 3

selc={vpm_PR_withT;pom_PR_withT};
splitRateBins=[-params_rate_time_before_trig 0 params_rate_time_after_trig];
split_var_nameC={'Acceleration';'Curvature';'Interval';'Velocity';'Setp';'Amp';'Angle';'Phase'};
type_nameC={'rawbefore','rawafter','absbefore','absafter','rel','absrel'};
plotAs='modulation';
%plotAs='rawresponse';
min_trial_per_quantile=0;
[RasterPT,trial_kinematics,splitPlotValues]=split_by_kinematics(trigsP,trigsT,min_trial_per_quantile,params_split_behav_time_before_trig,params_split_behav_after_trig, params_split_binsize,params_smoothwidth,params_splitsmoothmeth,splitRateBins,selc,params_splits,params_splitname,DiscreteData,ppms,split_var_nameC,type_nameC,params_split_allR,figdir,Amp,Phase,Angle,Setp,Velocity,Acceleration,Curvature,plotAs);


splitPlotValuesNames={'kinematic param','normalization', 'tertile cutoffs' ,...
    'vpm puff means', 'vpm puff sem','vpm p(puff=0)' ,'vpm touch means', 'vpm touch sem','vpm p(touch=0)','vpm p puff vs touch','vpm puff/touch p 1 vs 3',...
'pom puff means', 'pom puff sem','pom p(puff=0)','pom touch means', 'pom touch sem','pom p(touch=0)','pom p puff vs touch','pom puff/touch p 1 vs 3'};
%
SplitTable=cell2table(splitPlotValues,'VariableNames',splitPlotValuesNames);
writetable(SplitTable,[DataPath 'kinematicsplitstats_' datestr(now,'yymmdd') '.xlsx'])




figure('Position',[1000         974         998         364])
params_min_Touch_bout_interval=500;
params_min_Touch_bout_trials_per_group=10;
[mean_sem_p_touchboutsplit,H_TouchBoutSplit,TouchBoutSplitP,TouchBoutSplitN]=get_touch_num_in_bout_split({vpm_PR_withT;pom_PR_withT},params_min_Touch_bout_interval*ppms,params_min_Touch_bout_trials_per_group,trigsT,DiscreteData,ppms,params_rate_time_before_trig,params_rate_time_after_trig,params_rate_binsize,params_rate_baselinebinNum);
exportgraphics(gcf,[figdir 'Rates_split_by_Touch_in_bout_' datestr(now,'yymmdd') '.pdf'],'BackgroundColor','none');

Ncellssplit=cat(1,sum(all(TouchBoutSplitN(vpm_PR_withT,:)>=params_min_Touch_bout_trials_per_group,2)),sum(all(TouchBoutSplitN(pom_PR_withT,:)>=params_min_Touch_bout_trials_per_group,2)));
splitTouchBoutVals=cat(2,{'VPM';'POm'},num2cell(cat(2,Ncellssplit,mean_sem_p_touchboutsplit(:,:,1),mean_sem_p_touchboutsplit(:,:,2),mean_sem_p_touchboutsplit(:,:,3))));
splitTouchBoutNames={'Nucleus','N included',...
    'mean rate base' ,'mean rate 1st','mean rate next',...
    'sem rate base' ,'sem rate 1st','sem rate next',...
    'p 1st vs base','p next vs. base','p next vs 1st.'};
%
SplitTouchTable=cell2table(splitTouchBoutVals,'VariableNames',splitTouchBoutNames);
writetable(SplitTouchTable,[DataPath 'touchsplitstats_' datestr(now,'yymmdd') '.xlsx'])

%% Figure 5
Fig_no='Fig5_';
figure('Position',[265         633        1220         587])
tt=tiledlayout(1,3);
title(tt,'Whisking Light Response')
nexttile(1);
[p_QLWL_vs_QnLWnL{1},p_QL_vs_WL{1},mean_QLWL{1},sem_QLWL{1}]=plot_summary(cat(3,cat(2,RatesQ,RatesW),cat(2,RatesQL,RatesWL)),cat(2,abs(RatesWp)<.05,abs(RatesWLp)<.05),{'Q';'W';'QL';'WL'},vpm_PR_anyT_withL,{'nL';'L'},'mean','Rate [Hz]');
title('VPM')
nexttile(2);
[p_QLWL_vs_QnLWnL{2},p_QL_vs_WL{2},mean_QLWL{2},sem_QLWL{2}]=plot_summary(cat(3,cat(2,RatesQ,RatesW),cat(2,RatesQL,RatesWL)),cat(2,abs(RatesWp)<.05,abs(RatesWLp)<.05),{'Q';'W';'QL';'WL'},pom_PR_anyT_withL,{'nL';'L'},'mean','Rate [Hz]');
title('POm')
nexttile(3);
[mean_mod_vpmpom_QLWL,sem_mod_vpmpom_QLWL,p_comp_mod_acrossQLWL_vpmpom,p_comp_mod_acrossvpmpom_QLWL,p_comp_mod_vs_0_vpmpom_QLWL]=...
    plot_summary_across_groups(cat(2,comp_modulation(RatesQ,RatesW),comp_modulation(RatesQL,RatesWL)),RatesWLp<.05,{'nL: Q vs W','L: Q vs W'},{vpm_PR_withT;pom_PR_withT},{'VPM';'POm'},'mean','modulation');title('Modulation')
if write_figs;exportgraphics (gcf,[figdir Fig_no 'QWL_Summaries.pdf'], 'BackgroundColor','none','ContentType','vector');end

cell_select_nLL={vpm_PR_anyT_withL,pom_PR_anyT_withL};
figure('Position',[ 222         460        1018         518],'Color','white')
tt=tiledlayout (1,2);
title(tt,'Preferred phase & SNR with Light')
MarkerType={'mo','ro';'cd','bd'};
Nuc={'VPM';'POm'};
Markersize=[20 20];
for NucInd=1:2
    nexttile;
    for LightStateInd=1:2
        polarscatter(pref_phase_nLL(cell_select_nLL{NucInd}(PhaseMod_p_nLL(cell_select_nLL{NucInd},LightStateInd)>.05),LightStateInd),...
            Phase_SNR_nLL(cell_select_nLL{NucInd}(PhaseMod_p_nLL(cell_select_nLL{NucInd},LightStateInd)>.05),LightStateInd),Markersize(LightStateInd),MarkerType{LightStateInd,1},'filled');hold on
        polarscatter(pref_phase_nLL(cell_select_nLL{NucInd}(PhaseMod_p_nLL(cell_select_nLL{NucInd},LightStateInd)<=.05),LightStateInd),...
            Phase_SNR_nLL(cell_select_nLL{NucInd}(PhaseMod_p_nLL(cell_select_nLL{NucInd},LightStateInd)<=.05),LightStateInd),Markersize(LightStateInd),MarkerType{LightStateInd,2},'filled');hold on
    end
    polarplot(pref_phase_nLL(cell_select_nLL{NucInd},:)',Phase_SNR_nLL(cell_select_nLL{NucInd},:)' ,'k');hold on
    set(gca,'RLim', [0 2.5000])
    title(Nuc{NucInd})
end
if write_figs;exportgraphics (gcf,[figdir Fig_no 'phase_mod_change_whisking_in_air_light.pdf'], 'BackgroundColor','none','ContentType','vector');end

figure('Position',[ 910   176   500   588])
plot_summary_across_groups(Phase_SNR_nLL,PhaseMod_p_nLL(:,2)<=.05,{'no Light','Light'},{vpm_PR_anyT_withL;pom_PR_anyT_withL},{'VPM';'POm'},'mean','SNR');title('nL/L phase SNR')
if write_figs;exportgraphics (gcf,[figdir Fig_no 'Phase_SNR_Summaries.pdf'], 'BackgroundColor','none','ContentType','vector');end

%puff nL psths
contPlot={MeanPTrigAngle,MeanPtrigPuffState,MeanPTrigAmp};
Name='Puff no Light';
figure('Position',[38 260 1413 525],'Color','White');
[popPSTHs_Puff_nL]=plot_population_psth(trigsP,DiscreteData,ppms,params_psth_time_before_trig,params_psth_time_after_trig,params_psth_binsize,params_gausssmooth_PSTH_win,Name,{vpm_PR_anyT_withL;pom_PR_anyT_withL},{'VPM';'POm'},'latency',-inf,contPlot,{'meanbefore','scaleminmax','mono'},'baseline corr. whisker angle');
if write_figs;exportgraphics(gcf,[figdir Fig_no 'Pop_' parameter '_PSTH_' datestr(now,'yymmdd') '.pdf'],'BackgroundColor','none');end

%puff L psths
Name='Puff Light';
contPlot={MeanPLTrigAngle,MeanPLtrigPuffState,MeanPLTrigAmp};
figure('Position',[38 260 1413 525],'Color','White');
[popPSTHs_Puff_L]=plot_population_psth(trigsPL,DiscreteData,ppms,params_psth_time_before_trig,params_psth_time_after_trig,params_psth_binsize,params_gausssmooth_PSTH_win,Name,{vpm_PR_anyT_withL;pom_PR_anyT_withL},{'VPM';'POm'},'latency',-inf,contPlot,{'meanbefore','scaleminmax','mono'},'baseline corr. whisker angle');
if write_figs;exportgraphics(gcf,[figdir Fig_no 'Pop_' parameter '_PSTH_' datestr(now,'yymmdd') '.pdf'],'BackgroundColor','none','ContentType','vector');end

% summary puff nL/L
figure('Position',[265         633        1220         587])
tt=tiledlayout(1,3);
title(tt,'Puff/Light Response')
nexttile(1);
[p_PL_vs_baseline{1},p_P_vs_PL{1},mean_PL{1},sem_PL{1}]=plot_summary(cat(3,PuffTrigRatePrePost,LPuffTrigRatePrePost),cat(2,PuffTrigRateSigIncrease,LPuffTrigRateSigIncrease),{'base';'resp';'base';'resp'},vpm_PR_anyT_withL,{'nL Puff','L Puff'},'mean','Rate [Hz]');
title('VPM')
nexttile(2);
[p_PL_vs_baseline{2},p_P_vs_PL{2},mean_PL{2},sem_PL{2}]=plot_summary(cat(3,PuffTrigRatePrePost,LPuffTrigRatePrePost),cat(2,PuffTrigRateSigIncrease,LPuffTrigRateSigIncrease),{'base';'resp';'base';'resp'},pom_PR_anyT_withL,{'nL Puff','L Puff'},'mean','Rate [Hz]');
title('POm')
nexttile(3);
[mean_mod_vpmpom_PL,sem_mod_vpmpom_PL,p_comp_mod_acrossPL_vpmpom,p_comp_mod_acrossvpmpom_PL,p_comp_mod_vs_0_vpmpom_PL]=...
    plot_summary_across_groups(cat(2,comp_modulation(PuffTrigRatePrePost),comp_modulation(LPuffTrigRatePrePost)),true(size(PuffTrigRateSigIncrease)),{'nL Puff','L Puff'},{vpm_PR_anyT_withL;pom_PR_anyT_withL},{'VPM';'POm'},'mean','modulation');title('Modulation')
if write_figs;exportgraphics (gcf,[figdir Fig_no 'Puff_Light_Summaries.pdf'], 'BackgroundColor','none','ContentType','vector');end


%% Supplementary Figure 1 


% summary no puff / touch
figure('Position',[265         633        1220         587])
tt=tiledlayout(1,3);
title(tt,'non-puff-responsive neurons')
nexttile(1);
plot_summary(cat(3,PuffTrigRatePrePost,TouchTrigRatePrePost),cat(2,PuffTrigRateSigIncrease,TouchTrigRateSigIncrease),{'0';'1';'0';'1'},vpm_nPR_withT,{'Puff','Touch'},'mean','Rate [Hz]');
title('Puff/Touch Response VPM')
nexttile(2);
plot_summary(cat(3,PuffTrigRatePrePost,TouchTrigRatePrePost),cat(2,PuffTrigRateSigIncrease,TouchTrigRateSigIncrease),{'0';'1';'0';'1'},pom_nPR_withT,{'Puff','Touch'},'mean','Rate [Hz]');
title('Puff/Touch Response POm')
nexttile(3);
plot_summary_across_groups(cat(2,comp_modulation(PuffTrigRatePrePost),comp_modulation(TouchTrigRatePrePost)),true(size(PuffTrigRateSigIncrease)),{'Puff','Touch'},{vpm_nPR_withT;pom_nPR_withT},{'VPM';'POm'},'mean','mod depth');
title('Touch/Puff modulation')
 if write_figs;exportgraphics (gcf,[figdir 'SuppNonPuffRespTouch.pdf'], 'BackgroundColor','none','ContentType','vector');end




%% Response Summary Table
VN={'condition 1','condition 2','VPM cell IDs','POm cell IDs', 'VPM N','VPM N Animals','POm N','POm N Animals','Total N Animals',...
    'VPM mean condition 1 base','VPM sem condition 1 base','VPM mean condition 1 resp','VPM sem condition 1 resp','VPM p-value condition 1 response vs baseline',...
    'VPM mean condition 2 base','VPM sem condition 2 base','VPM mean condition 2 resp','VPM sem condition 2 resp','VPM p-value condition 2 response vs baseline',...
    'VPM p-value baseline across conditions','VPM p-value response across conditions',...
    'POm mean condition 1 base','POm sem condition 1 base','POm mean condition 1 resp','POm sem condition 1 resp','POm p-value condition 1 response vs baseline',...
    'POm mean condition 2 base','POm sem condition 2 base','POm mean condition 2 resp','POm sem condition 2 resp','POm p-value condition 2 response vs baseline',...
    'POm p-value baseline across conditions','POm p-value response across conditions',...
    'VPM mean modulation condition 1','VPM sem modulation condition 1','VPM p-value modulation cond1 vs zero','VPM mean modulation condition 2','VPM sem modulation condition 2','VPM p-value modulation cond2 vs zero','VPM p-value modulation across conditions',...
    'POm mean modulation condition 1','POm sem modulation condition 1','POm p-value modulation cond1 vs zero','POm mean modulation condition 2','POm sem modulation condition 2','POm p-value modulation cond2 vs zero','POm p-value modulation across conditions',...
    'p-value modulation condition 1 across nuclei', 'p-value modulation condition 2 across nuclei',...
    'VPM mean fsl condition 1','VPM sem fsl condition 1','VPM mean fsl condition 2','VPM sem fsl condition 2','VPM p-value fsl across conditions',...
    'POm mean fsl condition 1','POm sem fsl condition 1','POm mean fsl condition 2','POm sem fsl condition 2','POm p-value fsl across conditions',...
    'p-value fsl condition 1 across nuclei', 'p-value fsl condition 2 across nuclei'}';


SumStats=struct();
A={'N neurons'
'N animals'
'Baseline mean'
'Baseline SEM'
'p-Value Baseline across'
'Response mean'
'Response SEM'
'p-Value Response across'
'p-Value Response  = Baseline'
'Mean modulation'
'Mean modulation SEM'
'p-Value modulation = zero'
'p-Value modulation across'
'First spike latency mean'
'First spike latency SEM'
'p-Value First spike latency across'
'N neurons'
'N animals'
'Baseline mean'
'Baseline SEM'
'p-Value Baseline across'
'Response mean'
'Response SEM'
'p-Value Response across'
'p-Value Response  = Baseline'
'Mean modulation'
'Mean modulation SEM'
'p-Value modulation = zero'
'p-Value modulation across'
'First spike latency mean'
'First spike latency SEM'
'p-Value First spike latency across'
'N animals total'
'p-Value Baseline'
'p-Value Response'
'p-Value modulation'
'p-Value First spike latency'};
A=cat(2,cat(1,repmat({'VPM'},16,1),repmat({'POm'},16,1),repmat({'VPM vs. POm'},5,1)),A,num2cell(nan(37,10)));
SumStatTable=cell2table(A,'VariableNames',{'Nucleus','ValName','Puff','Touch','Q-Puff','W-Puff','Quiescence spont','Whisking spont','nL-Puff','L-Puff','nL-Whisking','L-Whisking'});

SumStatLine=1;SumTableCol=3;
SumStats(SumStatLine).condition_1='Puff';
SumStats(SumStatLine).condition_2='Touch';
SumStats(SumStatLine).VPM_cell_IDs=vpm_PR_withT';
SumStats(SumStatLine).POm_cell_IDs=pom_PR_withT';
[N_animals]=get_N_animals(SumStats(SumStatLine).VPM_cell_IDs,SumStats(SumStatLine).POm_cell_IDs,RecordingInfo);
SumStats(SumStatLine).VPM_N=numel(vpm_PR_withT);
SumStats(SumStatLine).VPM_N_animals=N_animals(3);
SumStats(SumStatLine).POm_N=numel(pom_PR_withT);
SumStats(SumStatLine).POm_N_animals=N_animals(4);
SumStats(SumStatLine).Total_N_animals=N_animals(1);
SumStats(SumStatLine).VPM_mean_cond1_base=mean_PT{1}(1,1);
SumStats(SumStatLine).VPM_sem_cond1_base=sem_PT{1}(1,1);
SumStats(SumStatLine).VPM_mean_cond1_resp=mean_PT{1}(1,2);
SumStats(SumStatLine).VPM_sem_cond1_resp=sem_PT{1}(1,2);
SumStats(SumStatLine).VPM_p_cond1=p_PT_vs_baseline{1}(1);
SumStats(SumStatLine).VPM_mean_cond2_base=mean_PT{1}(2,1);
SumStats(SumStatLine).VPM_sem_cond2_base=sem_PT{1}(2,1);
SumStats(SumStatLine).VPM_mean_cond2_resp=mean_PT{1}(2,2);
SumStats(SumStatLine).VPM_sem_cond2_resp=sem_PT{1}(2,2);
SumStats(SumStatLine).VPM_p_cond2=p_PT_vs_baseline{1}(2);
SumStats(SumStatLine).VPM_p_base_cond12=p_P_vs_T{1}(1);
SumStats(SumStatLine).VPM_p_resp_cond12=p_P_vs_T{1}(2);
SumStats(SumStatLine).POm_mean_cond1_base=mean_PT{2}(1,1);
SumStats(SumStatLine).POm_sem_cond1_base=sem_PT{2}(1,1);
SumStats(SumStatLine).POm_mean_cond1_resp=mean_PT{2}(1,2);
SumStats(SumStatLine).POm_sem_cond1_resp=sem_PT{2}(1,2);
SumStats(SumStatLine).POm_p_cond1=p_PT_vs_baseline{2}(1);
SumStats(SumStatLine).POm_mean_cond2_base=mean_PT{2}(2,1);
SumStats(SumStatLine).POm_sem_cond2_base=sem_PT{2}(2,1);
SumStats(SumStatLine).POm_mean_cond2_resp=mean_PT{2}(2,2);
SumStats(SumStatLine).POm_sem_cond2_resp=sem_PT{2}(2,2);
SumStats(SumStatLine).POm_p_cond2=p_PT_vs_baseline{2}(2);
SumStats(SumStatLine).POm_p_base_cond12=p_P_vs_T{2}(1);
SumStats(SumStatLine).POm_p_resp_cond12=p_P_vs_T{2}(2);
SumStats(SumStatLine).VPM_mean_modulation_cond1=mean_mod_vpmpom_PT(1,1);
SumStats(SumStatLine).VPM_sem_modulation_cond1=sem_mod_vpmpom_PT(1,1);
SumStats(SumStatLine).VPM_p_modulation_cond10=p_comp_mod_vs_0_vpmpom_PT(1,1);
SumStats(SumStatLine).VPM_mean_modulation_cond2=mean_mod_vpmpom_PT(1,2);
SumStats(SumStatLine).VPM_sem_modulation_cond2=sem_mod_vpmpom_PT(1,2);
SumStats(SumStatLine).VPM_p_modulation_cond20=p_comp_mod_vs_0_vpmpom_PT(1,2);
SumStats(SumStatLine).VPM_p_modulation_cond12=p_comp_mod_acrossPT_vpmpom(1);
SumStats(SumStatLine).POm_mean_modulation_cond1=mean_mod_vpmpom_PT(2,1);
SumStats(SumStatLine).POm_sem_modulation_cond1=sem_mod_vpmpom_PT(2,1);
SumStats(SumStatLine).POm_p_modulation_cond10=p_comp_mod_vs_0_vpmpom_PT(2,1);
SumStats(SumStatLine).POm_mean_modulation_cond2=mean_mod_vpmpom_PT(2,2);
SumStats(SumStatLine).POm_sem_modulation_cond2=sem_mod_vpmpom_PT(2,2);
SumStats(SumStatLine).POm_p_modulation_cond20=p_comp_mod_vs_0_vpmpom_PT(2,2);
SumStats(SumStatLine).POm_p_modulation_cond12=p_comp_mod_acrossPT_vpmpom(2);
SumStats(SumStatLine).VPMPOm_p_modulation_cond1=p_comp_mod_acrossvpmpom_PT(1);
SumStats(SumStatLine).VPMPOm_p_modulation_cond2=p_comp_mod_acrossvpmpom_PT(2);
SumStats(SumStatLine).VPMPOm_p_cond1_base=p_vpmpom_PT_baseline_resp(1,1);
SumStats(SumStatLine).VPMPOm_p_cond2_base=p_vpmpom_PT_baseline_resp(2,1);
SumStats(SumStatLine).VPMPOm_p_cond1_resp=p_vpmpom_PT_baseline_resp(1,2);
SumStats(SumStatLine).VPMPOm_p_cond2_resp=p_vpmpom_PT_baseline_resp(2,2);
SumStats(SumStatLine).VPM_mean_fsl_cond1=mean_fsl_vpmpom_PT(1,1);
SumStats(SumStatLine).VPM_sem_fsl_cond1=sem_fsl_vpmpom_PT(1,1);
SumStats(SumStatLine).VPM_mean_fsl_cond2=mean_fsl_vpmpom_PT(1,2);
SumStats(SumStatLine).VPM_sem_fsl_cond2=sem_fsl_vpmpom_PT(1,2);
SumStats(SumStatLine).VPM_p_fsl_cond12=p_comp_fsl_acrosscond_vpmpom_PT(1);
SumStats(SumStatLine).POm_mean_fsl_cond1=mean_fsl_vpmpom_PT(2,1);
SumStats(SumStatLine).POm_sem_fsl_cond1=sem_fsl_vpmpom_PT(2,1);
SumStats(SumStatLine).POm_mean_fsl_cond2=mean_fsl_vpmpom_PT(2,2);
SumStats(SumStatLine).POm_sem_fsl_cond2=sem_fsl_vpmpom_PT(2,2);
SumStats(SumStatLine).POm_p_fsl_cond12=p_comp_fsl_acrosscond_vpmpom_PT(2);
SumStats(SumStatLine).VPMPOm_p_fsl_cond1=p_comp_fsl_acrossvpmpom_PT(1);
SumStats(SumStatLine).VPMPOm_p_fsl_cond2=p_comp_fsl_acrossvpmpom_PT(2);

SumStatTable{1,SumTableCol}=SumStats(SumStatLine).VPM_N;
SumStatTable{2,SumTableCol}=SumStats(SumStatLine).VPM_N_animals;
SumStatTable{17,SumTableCol}=SumStats(SumStatLine).POm_N;
SumStatTable{18,SumTableCol}=SumStats(SumStatLine).POm_N_animals;
SumStatTable{33,SumTableCol}=SumStats(SumStatLine).Total_N_animals;
SumStatTable{3,SumTableCol}=SumStats(SumStatLine).VPM_mean_cond1_base;
SumStatTable{4,SumTableCol}=SumStats(SumStatLine).VPM_sem_cond1_base;
SumStatTable{6,SumTableCol}=SumStats(SumStatLine).VPM_mean_cond1_resp;
SumStatTable{7,SumTableCol}=SumStats(SumStatLine).VPM_sem_cond1_resp;
SumStatTable{9,SumTableCol}=SumStats(SumStatLine).VPM_p_cond1;
SumStatTable{3,SumTableCol+1}=SumStats(SumStatLine).VPM_mean_cond2_base;
SumStatTable{4,SumTableCol+1}=SumStats(SumStatLine).VPM_sem_cond2_base;
SumStatTable{6,SumTableCol+1}=SumStats(SumStatLine).VPM_mean_cond2_resp;
SumStatTable{7,SumTableCol+1}=SumStats(SumStatLine).VPM_sem_cond2_resp;
SumStatTable{9,SumTableCol+1}=SumStats(SumStatLine).VPM_p_cond2;
SumStatTable{5,SumTableCol}=SumStats(SumStatLine).VPM_p_base_cond12;
SumStatTable{8,SumTableCol}=SumStats(SumStatLine).VPM_p_resp_cond12;
SumStatTable{16+3,SumTableCol}=SumStats(SumStatLine).POm_mean_cond1_base;
SumStatTable{16+4,SumTableCol}=SumStats(SumStatLine).POm_sem_cond1_base;
SumStatTable{16+6,SumTableCol}=SumStats(SumStatLine).POm_mean_cond1_resp;
SumStatTable{16+7,SumTableCol}=SumStats(SumStatLine).POm_sem_cond1_resp;
SumStatTable{16+9,SumTableCol}=SumStats(SumStatLine).POm_p_cond1;
SumStatTable{16+3,SumTableCol+1}=SumStats(SumStatLine).POm_mean_cond2_base;
SumStatTable{16+4,SumTableCol+1}=SumStats(SumStatLine).POm_sem_cond2_base;
SumStatTable{16+6,SumTableCol+1}=SumStats(SumStatLine).POm_mean_cond2_resp;
SumStatTable{16+7,SumTableCol+1}=SumStats(SumStatLine).POm_sem_cond2_resp;
SumStatTable{16+9,SumTableCol+1}=SumStats(SumStatLine).POm_p_cond2;
SumStatTable{16+5,SumTableCol}=SumStats(SumStatLine).POm_p_base_cond12;
SumStatTable{16+8,SumTableCol}=SumStats(SumStatLine).POm_p_resp_cond12;
SumStatTable{10,SumTableCol}=SumStats(SumStatLine).VPM_mean_modulation_cond1;
SumStatTable{11,SumTableCol}=SumStats(SumStatLine).VPM_sem_modulation_cond1;
SumStatTable{12,SumTableCol}=SumStats(SumStatLine).VPM_p_modulation_cond10;
SumStatTable{10,SumTableCol+1}=SumStats(SumStatLine).VPM_mean_modulation_cond2;
SumStatTable{11,SumTableCol+1}=SumStats(SumStatLine).VPM_sem_modulation_cond2;
SumStatTable{12,SumTableCol+1}=SumStats(SumStatLine).VPM_p_modulation_cond20;
SumStatTable{13,SumTableCol}=SumStats(SumStatLine).VPM_p_modulation_cond12;
SumStatTable{16+10,SumTableCol}=SumStats(SumStatLine).POm_mean_modulation_cond1;
SumStatTable{16+11,SumTableCol}=SumStats(SumStatLine).POm_sem_modulation_cond1;
SumStatTable{16+12,SumTableCol}=SumStats(SumStatLine).POm_p_modulation_cond10;
SumStatTable{16+10,SumTableCol+1}=SumStats(SumStatLine).POm_mean_modulation_cond2;
SumStatTable{16+11,SumTableCol+1}=SumStats(SumStatLine).POm_sem_modulation_cond2;
SumStatTable{16+12,SumTableCol+1}=SumStats(SumStatLine).POm_p_modulation_cond20;
SumStatTable{16+13,SumTableCol}=SumStats(SumStatLine).POm_p_modulation_cond12;
SumStatTable{36,SumTableCol}=SumStats(SumStatLine).VPMPOm_p_modulation_cond1;
SumStatTable{36,SumTableCol+1}=SumStats(SumStatLine).VPMPOm_p_modulation_cond2;
SumStatTable{34,SumTableCol}=SumStats(SumStatLine).VPMPOm_p_cond1_base;
SumStatTable{34,SumTableCol+1}=SumStats(SumStatLine).VPMPOm_p_cond2_base;
SumStatTable{35,SumTableCol}=SumStats(SumStatLine).VPMPOm_p_cond1_resp;
SumStatTable{35,SumTableCol+1}=SumStats(SumStatLine).VPMPOm_p_cond2_resp;
SumStatTable{14,SumTableCol}=SumStats(SumStatLine).VPM_mean_fsl_cond1;
SumStatTable{15,SumTableCol}=SumStats(SumStatLine).VPM_sem_fsl_cond1;
SumStatTable{14,SumTableCol+1}=SumStats(SumStatLine).VPM_mean_fsl_cond2;
SumStatTable{15,SumTableCol+1}=SumStats(SumStatLine).VPM_sem_fsl_cond2;
SumStatTable{16,SumTableCol}=SumStats(SumStatLine).VPM_p_fsl_cond12;
SumStatTable{16+14,SumTableCol}=SumStats(SumStatLine).POm_mean_fsl_cond1;
SumStatTable{16+15,SumTableCol}=SumStats(SumStatLine).POm_sem_fsl_cond1;
SumStatTable{16+14,SumTableCol+1}=SumStats(SumStatLine).POm_mean_fsl_cond2;
SumStatTable{16+15,SumTableCol+1}=SumStats(SumStatLine).POm_sem_fsl_cond2;
SumStatTable{16+16,SumTableCol}=SumStats(SumStatLine).POm_p_fsl_cond12;
SumStatTable{37,SumTableCol}=SumStats(SumStatLine).VPMPOm_p_fsl_cond1;
SumStatTable{37,SumTableCol+1}=SumStats(SumStatLine).VPMPOm_p_fsl_cond2;

SumStatLine=2;SumTableCol=8;
SumStats(SumStatLine).condition_1='Quiescience';
SumStats(SumStatLine).condition_2='Whisking';
SumStats(SumStatLine).VPM_cell_IDs=vpm_PR_anyT';
SumStats(SumStatLine).POm_cell_IDs=pom_PR_anyT';
[N_animals]=get_N_animals(SumStats(SumStatLine).VPM_cell_IDs,SumStats(SumStatLine).POm_cell_IDs,RecordingInfo);
SumStats(SumStatLine).VPM_N=numel(vpm_PR_anyT);
SumStats(SumStatLine).VPM_N_animals=N_animals(3);
SumStats(SumStatLine).POm_N=numel(pom_PR_anyT);
SumStats(SumStatLine).POm_N_animals=N_animals(4);
SumStats(SumStatLine).Total_N_animals=N_animals(1);
SumStats(SumStatLine).VPM_mean_cond1_resp=mean_rate_vpmpom_QW(1,1);
SumStats(SumStatLine).VPM_sem_cond1_resp=sem_rate_vpmpom_QW(1,1);
SumStats(SumStatLine).VPM_mean_cond2_resp=mean_rate_vpmpom_QW(1,2);
SumStats(SumStatLine).VPM_sem_cond2_resp=sem_rate_vpmpom_QW(1,2);
SumStats(SumStatLine).VPM_p_resp_cond12=p_comp_rate_acrossQW_vpmpom(1);
SumStats(SumStatLine).POm_mean_cond1_resp=mean_rate_vpmpom_QW(2,1);
SumStats(SumStatLine).POm_sem_cond1_resp=sem_rate_vpmpom_QW(2,1);
SumStats(SumStatLine).POm_mean_cond2_resp=mean_rate_vpmpom_QW(2,2);
SumStats(SumStatLine).POm_sem_cond2_resp=sem_rate_vpmpom_QW(2,2);
SumStats(SumStatLine).POm_p_resp_cond12=p_comp_rate_acrossQW_vpmpom(2);
SumStats(SumStatLine).VPMPOm_p_cond1_resp=p_comp_rate_acrossvpmpom_QW(1);
SumStats(SumStatLine).VPMPOm_p_cond2_resp=p_comp_rate_acrossvpmpom_QW(2);
SumStats(SumStatLine).VPM_mean_modulation_cond2=mean_mod_vpmpom_QW(1,2);
SumStats(SumStatLine).VPM_sem_modulation_cond2=sem_mod_vpmpom_QW(1,2);
SumStats(SumStatLine).VPM_p_modulation_cond20=p_comp_mod_vs_0_vpmpom_QW(1,2);
SumStats(SumStatLine).POm_mean_modulation_cond2=mean_mod_vpmpom_QW(2,2);
SumStats(SumStatLine).POm_sem_modulation_cond2=sem_mod_vpmpom_QW(2,2);
SumStats(SumStatLine).POm_p_modulation_cond20=p_comp_mod_vs_0_vpmpom_QW(2,2);
SumStats(SumStatLine).VPMPOm_p_modulation_cond2=p_comp_mod_acrossvpmpom_QW(2);

SumStatTable{1,SumTableCol}=SumStats(SumStatLine).VPM_N;
SumStatTable{2,SumTableCol}=SumStats(SumStatLine).VPM_N_animals;
SumStatTable{17,SumTableCol}=SumStats(SumStatLine).POm_N;
SumStatTable{18,SumTableCol}=SumStats(SumStatLine).POm_N_animals;
SumStatTable{33,SumTableCol}=SumStats(SumStatLine).Total_N_animals;
SumStatTable{3,SumTableCol}=SumStats(SumStatLine).VPM_mean_cond1_resp;
SumStatTable{4,SumTableCol}=SumStats(SumStatLine).VPM_sem_cond1_resp;
SumStatTable{6,SumTableCol}=SumStats(SumStatLine).VPM_mean_cond2_resp;
SumStatTable{7,SumTableCol}=SumStats(SumStatLine).VPM_sem_cond2_resp;
 SumStatTable{9,SumTableCol}=SumStats(SumStatLine).VPM_p_resp_cond12;
SumStatTable{16+3,SumTableCol}=SumStats(SumStatLine).POm_mean_cond1_resp;
SumStatTable{16+4,SumTableCol}=SumStats(SumStatLine).POm_sem_cond1_resp;
SumStatTable{16+6,SumTableCol}=SumStats(SumStatLine).POm_mean_cond2_resp;
SumStatTable{16+7,SumTableCol}=SumStats(SumStatLine).POm_sem_cond2_resp;
SumStatTable{16+9,SumTableCol}=SumStats(SumStatLine).POm_p_resp_cond12;
SumStatTable{10,SumTableCol}=SumStats(SumStatLine).VPM_mean_modulation_cond2;
SumStatTable{11,SumTableCol}=SumStats(SumStatLine).VPM_sem_modulation_cond2;
SumStatTable{12,SumTableCol}=SumStats(SumStatLine).VPM_p_modulation_cond20;
SumStatTable{16+10,SumTableCol}=SumStats(SumStatLine).POm_mean_modulation_cond2;
SumStatTable{16+11,SumTableCol}=SumStats(SumStatLine).POm_sem_modulation_cond2;
SumStatTable{16+12,SumTableCol}=SumStats(SumStatLine).POm_p_modulation_cond20;
SumStatTable{36,SumTableCol}=SumStats(SumStatLine).VPMPOm_p_modulation_cond2;
SumStatTable{34,SumTableCol}=SumStats(SumStatLine).VPMPOm_p_cond1_resp;
SumStatTable{35,SumTableCol}=SumStats(SumStatLine).VPMPOm_p_cond2_resp;

SumStatLine=3;SumTableCol=5;
SumStats(SumStatLine).condition_1='Puff nonWhisking';
SumStats(SumStatLine).condition_2='Puff Whisking';
SumStats(SumStatLine).VPM_cell_IDs=vpm_PR_withT';
SumStats(SumStatLine).POm_cell_IDs=pom_PR_withT';
[N_animals]=get_N_animals(SumStats(SumStatLine).VPM_cell_IDs,SumStats(SumStatLine).POm_cell_IDs,RecordingInfo);
SumStats(SumStatLine).VPM_N=numel(vpm_PR_withT);
SumStats(SumStatLine).VPM_N_animals=N_animals(3);
SumStats(SumStatLine).POm_N=numel(pom_PR_withT);
SumStats(SumStatLine).POm_N_animals=N_animals(4);
SumStats(SumStatLine).Total_N_animals=N_animals(1);
SumStats(SumStatLine).VPM_mean_cond1_base=mean_QPWP{1}(1,1);
SumStats(SumStatLine).VPM_sem_cond1_base=sem_QPWP{1}(1,1);
SumStats(SumStatLine).VPM_mean_cond1_resp=mean_QPWP{1}(1,2);
SumStats(SumStatLine).VPM_sem_cond1_resp=sem_QPWP{1}(1,2);
SumStats(SumStatLine).VPM_p_cond1=p_QPWP_vs_baseline{1}(1);
SumStats(SumStatLine).VPM_mean_cond2_base=mean_QPWP{1}(2,1);
SumStats(SumStatLine).VPM_sem_cond2_base=sem_QPWP{1}(2,1);
SumStats(SumStatLine).VPM_mean_cond2_resp=mean_QPWP{1}(2,2);
SumStats(SumStatLine).VPM_sem_cond2_resp=sem_QPWP{1}(2,2);
SumStats(SumStatLine).VPM_p_cond2=p_QPWP_vs_baseline{1}(2);
SumStats(SumStatLine).VPM_p_base_cond12=p_QP_vs_WP{1}(1);
SumStats(SumStatLine).VPM_p_resp_cond12=p_QP_vs_WP{1}(2);
SumStats(SumStatLine).POm_mean_cond1_base=mean_QPWP{2}(1,1);
SumStats(SumStatLine).POm_sem_cond1_base=sem_QPWP{2}(1,1);
SumStats(SumStatLine).POm_mean_cond1_resp=mean_QPWP{2}(1,2);
SumStats(SumStatLine).POm_sem_cond1_resp=sem_QPWP{2}(1,2);
SumStats(SumStatLine).POm_p_cond1=p_QPWP_vs_baseline{2}(1);
SumStats(SumStatLine).POm_mean_cond2_base=mean_QPWP{2}(2,1);
SumStats(SumStatLine).POm_sem_cond2_base=sem_QPWP{2}(2,1);
SumStats(SumStatLine).POm_mean_cond2_resp=mean_QPWP{2}(2,2);
SumStats(SumStatLine).POm_sem_cond2_resp=sem_QPWP{2}(2,2);
SumStats(SumStatLine).POm_p_cond2=p_QPWP_vs_baseline{2}(2);
SumStats(SumStatLine).POm_p_base_cond12=p_QP_vs_WP{2}(1);
SumStats(SumStatLine).POm_p_resp_cond12=p_QP_vs_WP{2}(2);
SumStats(SumStatLine).VPM_mean_modulation_cond1=mean_mod_vpmpom_QPWP(1,1);
SumStats(SumStatLine).VPM_sem_modulation_cond1=sem_mod_vpmpom_QPWP(1,1);
SumStats(SumStatLine).VPM_p_modulation_cond10=p_comp_mod_vs_0_vpmpom_QPWP(1,1);
SumStats(SumStatLine).VPM_mean_modulation_cond2=mean_mod_vpmpom_QPWP(1,2);
SumStats(SumStatLine).VPM_sem_modulation_cond2=sem_mod_vpmpom_QPWP(1,2);
SumStats(SumStatLine).VPM_p_modulation_cond20=p_comp_mod_vs_0_vpmpom_QPWP(1,2);
SumStats(SumStatLine).VPM_p_modulation_cond12=p_comp_mod_acrossQPWP_vpmpom(1);
SumStats(SumStatLine).POm_mean_modulation_cond1=mean_mod_vpmpom_QPWP(2,1);
SumStats(SumStatLine).POm_sem_modulation_cond1=sem_mod_vpmpom_QPWP(2,1);
SumStats(SumStatLine).POm_p_modulation_cond10=p_comp_mod_vs_0_vpmpom_QPWP(2,1);
SumStats(SumStatLine).POm_mean_modulation_cond2=mean_mod_vpmpom_QPWP(2,2);
SumStats(SumStatLine).POm_sem_modulation_cond2=sem_mod_vpmpom_QPWP(2,2);
SumStats(SumStatLine).POm_p_modulation_cond20=p_comp_mod_vs_0_vpmpom_QPWP(2,2);
SumStats(SumStatLine).POm_p_modulation_cond12=p_comp_mod_acrossQPWP_vpmpom(2);
SumStats(SumStatLine).VPMPOm_p_modulation_cond1=p_comp_mod_acrossvpmpom_QPWP(1);
SumStats(SumStatLine).VPMPOm_p_modulation_cond2=p_comp_mod_acrossvpmpom_QPWP(2);
SumStats(SumStatLine).VPMPOm_p_cond1_base=p_vpmpom_QPWP_baseline_resp(1,1);
SumStats(SumStatLine).VPMPOm_p_cond2_base=p_vpmpom_QPWP_baseline_resp(2,1);
SumStats(SumStatLine).VPMPOm_p_cond1_resp=p_vpmpom_QPWP_baseline_resp(1,2);
SumStats(SumStatLine).VPMPOm_p_cond2_resp=p_vpmpom_QPWP_baseline_resp(2,2);
SumStats(SumStatLine).VPM_mean_fsl_cond1=mean_fsl_vpmpom_QPWP(1,1);
SumStats(SumStatLine).VPM_sem_fsl_cond1=sem_fsl_vpmpom_QPWP(1,1);
SumStats(SumStatLine).VPM_mean_fsl_cond2=mean_fsl_vpmpom_QPWP(1,2);
SumStats(SumStatLine).VPM_sem_fsl_cond2=sem_fsl_vpmpom_QPWP(1,2);
SumStats(SumStatLine).VPM_p_fsl_cond12=p_comp_fsl_acrosscond_vpmpom_QPWP(1);
SumStats(SumStatLine).POm_mean_fsl_cond1=mean_fsl_vpmpom_QPWP(2,1);
SumStats(SumStatLine).POm_sem_fsl_cond1=sem_fsl_vpmpom_QPWP(2,1);
SumStats(SumStatLine).POm_mean_fsl_cond2=mean_fsl_vpmpom_QPWP(2,2);
SumStats(SumStatLine).POm_sem_fsl_cond2=sem_fsl_vpmpom_QPWP(2,2);
SumStats(SumStatLine).POm_p_fsl_cond12=p_comp_fsl_acrosscond_vpmpom_QPWP(2);
SumStats(SumStatLine).VPMPOm_p_fsl_cond1=p_comp_fsl_acrossvpmpom_QPWP(1);
SumStats(SumStatLine).VPMPOm_p_fsl_cond2=p_comp_fsl_acrossvpmpom_QPWP(2);

SumStatTable{1,SumTableCol}=SumStats(SumStatLine).VPM_N;
SumStatTable{2,SumTableCol}=SumStats(SumStatLine).VPM_N_animals;
SumStatTable{17,SumTableCol}=SumStats(SumStatLine).POm_N;
SumStatTable{18,SumTableCol}=SumStats(SumStatLine).POm_N_animals;
SumStatTable{33,SumTableCol}=SumStats(SumStatLine).Total_N_animals;
SumStatTable{3,SumTableCol}=SumStats(SumStatLine).VPM_mean_cond1_base;
SumStatTable{4,SumTableCol}=SumStats(SumStatLine).VPM_sem_cond1_base;
SumStatTable{6,SumTableCol}=SumStats(SumStatLine).VPM_mean_cond1_resp;
SumStatTable{7,SumTableCol}=SumStats(SumStatLine).VPM_sem_cond1_resp;
SumStatTable{9,SumTableCol}=SumStats(SumStatLine).VPM_p_cond1;
SumStatTable{3,SumTableCol+1}=SumStats(SumStatLine).VPM_mean_cond2_base;
SumStatTable{4,SumTableCol+1}=SumStats(SumStatLine).VPM_sem_cond2_base;
SumStatTable{6,SumTableCol+1}=SumStats(SumStatLine).VPM_mean_cond2_resp;
SumStatTable{7,SumTableCol+1}=SumStats(SumStatLine).VPM_sem_cond2_resp;
SumStatTable{9,SumTableCol+1}=SumStats(SumStatLine).VPM_p_cond2;
SumStatTable{5,SumTableCol}=SumStats(SumStatLine).VPM_p_base_cond12;
SumStatTable{8,SumTableCol}=SumStats(SumStatLine).VPM_p_resp_cond12;
SumStatTable{16+3,SumTableCol}=SumStats(SumStatLine).POm_mean_cond1_base;
SumStatTable{16+4,SumTableCol}=SumStats(SumStatLine).POm_sem_cond1_base;
SumStatTable{16+6,SumTableCol}=SumStats(SumStatLine).POm_mean_cond1_resp;
SumStatTable{16+7,SumTableCol}=SumStats(SumStatLine).POm_sem_cond1_resp;
SumStatTable{16+9,SumTableCol}=SumStats(SumStatLine).POm_p_cond1;
SumStatTable{16+3,SumTableCol+1}=SumStats(SumStatLine).POm_mean_cond2_base;
SumStatTable{16+4,SumTableCol+1}=SumStats(SumStatLine).POm_sem_cond2_base;
SumStatTable{16+6,SumTableCol+1}=SumStats(SumStatLine).POm_mean_cond2_resp;
SumStatTable{16+7,SumTableCol+1}=SumStats(SumStatLine).POm_sem_cond2_resp;
SumStatTable{16+9,SumTableCol+1}=SumStats(SumStatLine).POm_p_cond2;
SumStatTable{16+5,SumTableCol}=SumStats(SumStatLine).POm_p_base_cond12;
SumStatTable{16+8,SumTableCol}=SumStats(SumStatLine).POm_p_resp_cond12;
SumStatTable{10,SumTableCol}=SumStats(SumStatLine).VPM_mean_modulation_cond1;
SumStatTable{11,SumTableCol}=SumStats(SumStatLine).VPM_sem_modulation_cond1;
SumStatTable{12,SumTableCol}=SumStats(SumStatLine).VPM_p_modulation_cond10;
SumStatTable{10,SumTableCol+1}=SumStats(SumStatLine).VPM_mean_modulation_cond2;
SumStatTable{11,SumTableCol+1}=SumStats(SumStatLine).VPM_sem_modulation_cond2;
SumStatTable{12,SumTableCol+1}=SumStats(SumStatLine).VPM_p_modulation_cond20;
SumStatTable{13,SumTableCol}=SumStats(SumStatLine).VPM_p_modulation_cond12;
SumStatTable{16+10,SumTableCol}=SumStats(SumStatLine).POm_mean_modulation_cond1;
SumStatTable{16+11,SumTableCol}=SumStats(SumStatLine).POm_sem_modulation_cond1;
SumStatTable{16+12,SumTableCol}=SumStats(SumStatLine).POm_p_modulation_cond10;
SumStatTable{16+10,SumTableCol+1}=SumStats(SumStatLine).POm_mean_modulation_cond2;
SumStatTable{16+11,SumTableCol+1}=SumStats(SumStatLine).POm_sem_modulation_cond2;
SumStatTable{16+12,SumTableCol+1}=SumStats(SumStatLine).POm_p_modulation_cond20;
SumStatTable{16+13,SumTableCol}=SumStats(SumStatLine).POm_p_modulation_cond12;
SumStatTable{36,SumTableCol}=SumStats(SumStatLine).VPMPOm_p_modulation_cond1;
SumStatTable{36,SumTableCol+1}=SumStats(SumStatLine).VPMPOm_p_modulation_cond2;
SumStatTable{34,SumTableCol}=SumStats(SumStatLine).VPMPOm_p_cond1_base;
SumStatTable{34,SumTableCol+1}=SumStats(SumStatLine).VPMPOm_p_cond2_base;
SumStatTable{35,SumTableCol}=SumStats(SumStatLine).VPMPOm_p_cond1_resp;
SumStatTable{35,SumTableCol+1}=SumStats(SumStatLine).VPMPOm_p_cond2_resp;
SumStatTable{14,SumTableCol}=SumStats(SumStatLine).VPM_mean_fsl_cond1;
SumStatTable{15,SumTableCol}=SumStats(SumStatLine).VPM_sem_fsl_cond1;
SumStatTable{14,SumTableCol+1}=SumStats(SumStatLine).VPM_mean_fsl_cond2;
SumStatTable{15,SumTableCol+1}=SumStats(SumStatLine).VPM_sem_fsl_cond2;
SumStatTable{16,SumTableCol}=SumStats(SumStatLine).VPM_p_fsl_cond12;
SumStatTable{16+14,SumTableCol}=SumStats(SumStatLine).POm_mean_fsl_cond1;
SumStatTable{16+15,SumTableCol}=SumStats(SumStatLine).POm_sem_fsl_cond1;
SumStatTable{16+14,SumTableCol+1}=SumStats(SumStatLine).POm_mean_fsl_cond2;
SumStatTable{16+15,SumTableCol+1}=SumStats(SumStatLine).POm_sem_fsl_cond2;
SumStatTable{16+16,SumTableCol}=SumStats(SumStatLine).POm_p_fsl_cond12;
SumStatTable{37,SumTableCol}=SumStats(SumStatLine).VPMPOm_p_fsl_cond1;
SumStatTable{37,SumTableCol+1}=SumStats(SumStatLine).VPMPOm_p_fsl_cond2;

SumStatLine=4;SumTableCol=11;
SumStats(SumStatLine).condition_1='whisking no Light';
SumStats(SumStatLine).condition_2='whisking Light';
SumStats(SumStatLine).VPM_cell_IDs=vpm_PR_anyT_withL';
SumStats(SumStatLine).POm_cell_IDs=pom_PR_anyT_withL';
[N_animals]=get_N_animals(SumStats(SumStatLine).VPM_cell_IDs,SumStats(SumStatLine).POm_cell_IDs,RecordingInfo);
SumStats(SumStatLine).VPM_N=numel(vpm_PR_anyT_withL);
SumStats(SumStatLine).VPM_N_animals=N_animals(3);
SumStats(SumStatLine).POm_N=numel(pom_PR_anyT_withL);
SumStats(SumStatLine).POm_N_animals=N_animals(4);
SumStats(SumStatLine).Total_N_animals=N_animals(1);
SumStats(SumStatLine).VPM_mean_cond1_base=mean_QLWL{1}(1,1);
SumStats(SumStatLine).VPM_sem_cond1_base=sem_QLWL{1}(1,1);
SumStats(SumStatLine).VPM_mean_cond1_resp=mean_QLWL{1}(1,2);
SumStats(SumStatLine).VPM_sem_cond1_resp=sem_QLWL{1}(1,2);
SumStats(SumStatLine).VPM_p_cond1=p_QLWL_vs_QnLWnL{1}(1);
SumStats(SumStatLine).VPM_mean_cond2_base=mean_QLWL{1}(2,1);
SumStats(SumStatLine).VPM_sem_cond2_base=sem_QLWL{1}(2,1);
SumStats(SumStatLine).VPM_mean_cond2_resp=mean_QLWL{1}(2,2);
SumStats(SumStatLine).VPM_sem_cond2_resp=sem_QLWL{1}(2,2);
SumStats(SumStatLine).VPM_p_cond2=p_QLWL_vs_QnLWnL{1}(2);
SumStats(SumStatLine).VPM_p_base_cond12=p_QL_vs_WL{1}(1);
SumStats(SumStatLine).VPM_p_resp_cond12=p_QL_vs_WL{1}(2);
SumStats(SumStatLine).POm_mean_cond1_base=mean_QLWL{2}(1,1);
SumStats(SumStatLine).POm_sem_cond1_base=sem_QLWL{2}(1,1);
SumStats(SumStatLine).POm_mean_cond1_resp=mean_QLWL{2}(1,2);
SumStats(SumStatLine).POm_sem_cond1_resp=sem_QLWL{2}(1,2);
SumStats(SumStatLine).POm_p_cond1=p_QLWL_vs_QnLWnL{2}(1);
SumStats(SumStatLine).POm_mean_cond2_base=mean_QLWL{2}(2,1);
SumStats(SumStatLine).POm_sem_cond2_base=sem_QLWL{2}(2,1);
SumStats(SumStatLine).POm_mean_cond2_resp=mean_QLWL{2}(2,2);
SumStats(SumStatLine).POm_sem_cond2_resp=sem_QLWL{2}(2,2);
SumStats(SumStatLine).POm_p_cond2=p_QLWL_vs_QnLWnL{2}(2);
SumStats(SumStatLine).POm_p_base_cond12=p_QL_vs_WL{2}(1);
SumStats(SumStatLine).POm_p_resp_cond12=p_QL_vs_WL{2}(2);
SumStats(SumStatLine).VPM_mean_modulation_cond1=mean_mod_vpmpom_QLWL(1,1);
SumStats(SumStatLine).VPM_sem_modulation_cond1=sem_mod_vpmpom_QLWL(1,1);
SumStats(SumStatLine).VPM_p_modulation_cond10=p_comp_mod_vs_0_vpmpom_QLWL(1,1);
SumStats(SumStatLine).VPM_mean_modulation_cond2=mean_mod_vpmpom_QLWL(1,2);
SumStats(SumStatLine).VPM_sem_modulation_cond2=sem_mod_vpmpom_QLWL(1,2);
SumStats(SumStatLine).VPM_p_modulation_cond20=p_comp_mod_vs_0_vpmpom_QLWL(1,2);
SumStats(SumStatLine).VPM_p_modulation_cond12=p_comp_mod_acrossQLWL_vpmpom(1);
SumStats(SumStatLine).POm_mean_modulation_cond1=mean_mod_vpmpom_QLWL(2,1);
SumStats(SumStatLine).POm_sem_modulation_cond1=sem_mod_vpmpom_QLWL(2,1);
SumStats(SumStatLine).POm_p_modulation_cond10=p_comp_mod_vs_0_vpmpom_QLWL(2,1);
SumStats(SumStatLine).POm_mean_modulation_cond2=mean_mod_vpmpom_QLWL(2,2);
SumStats(SumStatLine).POm_sem_modulation_cond2=sem_mod_vpmpom_QLWL(2,2);
SumStats(SumStatLine).POm_p_modulation_cond20=p_comp_mod_vs_0_vpmpom_QLWL(2,2);
SumStats(SumStatLine).POm_p_modulation_cond12=p_comp_mod_acrossQLWL_vpmpom(2);
SumStats(SumStatLine).VPMPOm_p_modulation_cond1=p_comp_mod_acrossvpmpom_QLWL(1);
SumStats(SumStatLine).VPMPOm_p_modulation_cond2=p_comp_mod_acrossvpmpom_QLWL(2);
SumStats(SumStatLine).VPMPOm_p_cond1_base=p_vpmpom_WLWnL_baseline_resp(1,1);
SumStats(SumStatLine).VPMPOm_p_cond2_base=p_vpmpom_WLWnL_baseline_resp(2,1);
SumStats(SumStatLine).VPMPOm_p_cond1_resp=p_vpmpom_WLWnL_baseline_resp(1,2);
SumStats(SumStatLine).VPMPOm_p_cond2_resp=p_vpmpom_WLWnL_baseline_resp(2,2);

SumStatTable{1,SumTableCol}=SumStats(SumStatLine).VPM_N;
SumStatTable{2,SumTableCol}=SumStats(SumStatLine).VPM_N_animals;
SumStatTable{17,SumTableCol}=SumStats(SumStatLine).POm_N;
SumStatTable{18,SumTableCol}=SumStats(SumStatLine).POm_N_animals;
SumStatTable{33,SumTableCol}=SumStats(SumStatLine).Total_N_animals;
SumStatTable{3,SumTableCol}=SumStats(SumStatLine).VPM_mean_cond1_base;
SumStatTable{4,SumTableCol}=SumStats(SumStatLine).VPM_sem_cond1_base;
SumStatTable{6,SumTableCol}=SumStats(SumStatLine).VPM_mean_cond1_resp;
SumStatTable{7,SumTableCol}=SumStats(SumStatLine).VPM_sem_cond1_resp;
SumStatTable{9,SumTableCol}=SumStats(SumStatLine).VPM_p_cond1;
SumStatTable{3,SumTableCol+1}=SumStats(SumStatLine).VPM_mean_cond2_base;
SumStatTable{4,SumTableCol+1}=SumStats(SumStatLine).VPM_sem_cond2_base;
SumStatTable{6,SumTableCol+1}=SumStats(SumStatLine).VPM_mean_cond2_resp;
SumStatTable{7,SumTableCol+1}=SumStats(SumStatLine).VPM_sem_cond2_resp;
SumStatTable{9,SumTableCol+1}=SumStats(SumStatLine).VPM_p_cond2;
SumStatTable{5,SumTableCol}=SumStats(SumStatLine).VPM_p_base_cond12;
SumStatTable{8,SumTableCol}=SumStats(SumStatLine).VPM_p_resp_cond12;
SumStatTable{16+3,SumTableCol}=SumStats(SumStatLine).POm_mean_cond1_base;
SumStatTable{16+4,SumTableCol}=SumStats(SumStatLine).POm_sem_cond1_base;
SumStatTable{16+6,SumTableCol}=SumStats(SumStatLine).POm_mean_cond1_resp;
SumStatTable{16+7,SumTableCol}=SumStats(SumStatLine).POm_sem_cond1_resp;
SumStatTable{16+9,SumTableCol}=SumStats(SumStatLine).POm_p_cond1;
SumStatTable{16+3,SumTableCol+1}=SumStats(SumStatLine).POm_mean_cond2_base;
SumStatTable{16+4,SumTableCol+1}=SumStats(SumStatLine).POm_sem_cond2_base;
SumStatTable{16+6,SumTableCol+1}=SumStats(SumStatLine).POm_mean_cond2_resp;
SumStatTable{16+7,SumTableCol+1}=SumStats(SumStatLine).POm_sem_cond2_resp;
SumStatTable{16+9,SumTableCol+1}=SumStats(SumStatLine).POm_p_cond2;
SumStatTable{16+5,SumTableCol}=SumStats(SumStatLine).POm_p_base_cond12;
SumStatTable{16+8,SumTableCol}=SumStats(SumStatLine).POm_p_resp_cond12;
SumStatTable{10,SumTableCol}=SumStats(SumStatLine).VPM_mean_modulation_cond1;
SumStatTable{11,SumTableCol}=SumStats(SumStatLine).VPM_sem_modulation_cond1;
SumStatTable{12,SumTableCol}=SumStats(SumStatLine).VPM_p_modulation_cond10;
SumStatTable{10,SumTableCol+1}=SumStats(SumStatLine).VPM_mean_modulation_cond2;
SumStatTable{11,SumTableCol+1}=SumStats(SumStatLine).VPM_sem_modulation_cond2;
SumStatTable{12,SumTableCol+1}=SumStats(SumStatLine).VPM_p_modulation_cond20;
SumStatTable{13,SumTableCol}=SumStats(SumStatLine).VPM_p_modulation_cond12;
SumStatTable{16+10,SumTableCol}=SumStats(SumStatLine).POm_mean_modulation_cond1;
SumStatTable{16+11,SumTableCol}=SumStats(SumStatLine).POm_sem_modulation_cond1;
SumStatTable{16+12,SumTableCol}=SumStats(SumStatLine).POm_p_modulation_cond10;
SumStatTable{16+10,SumTableCol+1}=SumStats(SumStatLine).POm_mean_modulation_cond2;
SumStatTable{16+11,SumTableCol+1}=SumStats(SumStatLine).POm_sem_modulation_cond2;
SumStatTable{16+12,SumTableCol+1}=SumStats(SumStatLine).POm_p_modulation_cond20;
SumStatTable{16+13,SumTableCol}=SumStats(SumStatLine).POm_p_modulation_cond12;
SumStatTable{36,SumTableCol}=SumStats(SumStatLine).VPMPOm_p_modulation_cond1;
SumStatTable{36,SumTableCol+1}=SumStats(SumStatLine).VPMPOm_p_modulation_cond2;
SumStatTable{34,SumTableCol}=SumStats(SumStatLine).VPMPOm_p_cond1_base;
SumStatTable{34,SumTableCol+1}=SumStats(SumStatLine).VPMPOm_p_cond2_base;
SumStatTable{35,SumTableCol}=SumStats(SumStatLine).VPMPOm_p_cond1_resp;
SumStatTable{35,SumTableCol+1}=SumStats(SumStatLine).VPMPOm_p_cond2_resp;


SumStatLine=5;SumTableCol=9;
SumStats(SumStatLine).condition_1='Puff no Light';
SumStats(SumStatLine).condition_2='Puff Light';
SumStats(SumStatLine).VPM_cell_IDs=vpm_PR_anyT_withL';
SumStats(SumStatLine).POm_cell_IDs=pom_PR_anyT_withL';
SumStats(SumStatLine).VPM_N=numel(vpm_PR_anyT_withL);
SumStats(SumStatLine).POm_N=numel(pom_PR_anyT_withL);
[N_animals]=get_N_animals(SumStats(SumStatLine).VPM_cell_IDs,SumStats(SumStatLine).POm_cell_IDs,RecordingInfo);
SumStats(SumStatLine).VPM_N=numel(vpm_PR_anyT_withL);
SumStats(SumStatLine).VPM_N_animals=N_animals(3);
SumStats(SumStatLine).POm_N=numel(pom_PR_anyT_withL);
SumStats(SumStatLine).POm_N_animals=N_animals(4);
SumStats(SumStatLine).Total_N_animals=N_animals(1);
SumStats(SumStatLine).VPM_mean_cond1_base=mean_PL{1}(1,1);
SumStats(SumStatLine).VPM_sem_cond1_base=sem_PL{1}(1,1);
SumStats(SumStatLine).VPM_mean_cond1_resp=mean_PL{1}(1,2);
SumStats(SumStatLine).VPM_sem_cond1_resp=sem_PL{1}(1,2);
SumStats(SumStatLine).VPM_p_cond1=p_PL_vs_baseline{1}(1);
SumStats(SumStatLine).VPM_mean_cond2_base=mean_PL{1}(2,1);
SumStats(SumStatLine).VPM_sem_cond2_base=sem_PL{1}(2,1);
SumStats(SumStatLine).VPM_mean_cond2_resp=mean_PL{1}(2,2);
SumStats(SumStatLine).VPM_sem_cond2_resp=sem_PL{1}(2,2);
SumStats(SumStatLine).VPM_p_cond2=p_PL_vs_baseline{1}(2);
SumStats(SumStatLine).VPM_p_base_cond12=p_P_vs_PL{1}(1);
SumStats(SumStatLine).VPM_p_resp_cond12=p_P_vs_PL{1}(2);
SumStats(SumStatLine).POm_mean_cond1_base=mean_PL{2}(1,1);
SumStats(SumStatLine).POm_sem_cond1_base=sem_PL{2}(1,1);
SumStats(SumStatLine).POm_mean_cond1_resp=mean_PL{2}(1,2);
SumStats(SumStatLine).POm_sem_cond1_resp=sem_PL{2}(1,2);
SumStats(SumStatLine).POm_p_cond1=p_PL_vs_baseline{2}(1);
SumStats(SumStatLine).POm_mean_cond2_base=mean_PL{2}(2,1);
SumStats(SumStatLine).POm_sem_cond2_base=sem_PL{2}(2,1);
SumStats(SumStatLine).POm_mean_cond2_resp=mean_PL{2}(2,2);
SumStats(SumStatLine).POm_sem_cond2_resp=sem_PL{2}(2,2);
SumStats(SumStatLine).POm_p_cond2=p_PL_vs_baseline{2}(2);
SumStats(SumStatLine).POm_p_base_cond12=p_P_vs_PL{2}(1);
SumStats(SumStatLine).POm_p_resp_cond12=p_P_vs_PL{2}(2);
SumStats(SumStatLine).VPM_mean_modulation_cond1=mean_mod_vpmpom_PL(1,1);
SumStats(SumStatLine).VPM_sem_modulation_cond1=sem_mod_vpmpom_PL(1,1);
SumStats(SumStatLine).VPM_p_modulation_cond10=p_comp_mod_vs_0_vpmpom_PL(1,1);
SumStats(SumStatLine).VPM_mean_modulation_cond2=mean_mod_vpmpom_PL(1,2);
SumStats(SumStatLine).VPM_sem_modulation_cond2=sem_mod_vpmpom_PL(1,2);
SumStats(SumStatLine).VPM_p_modulation_cond20=p_comp_mod_vs_0_vpmpom_PL(1,2);
SumStats(SumStatLine).VPM_p_modulation_cond12=p_comp_mod_acrossPL_vpmpom(1);
SumStats(SumStatLine).POm_mean_modulation_cond1=mean_mod_vpmpom_PL(2,1);
SumStats(SumStatLine).POm_sem_modulation_cond1=sem_mod_vpmpom_PL(2,1);
SumStats(SumStatLine).POm_p_modulation_cond10=p_comp_mod_vs_0_vpmpom_PL(2,1);
SumStats(SumStatLine).POm_mean_modulation_cond2=mean_mod_vpmpom_PL(2,2);
SumStats(SumStatLine).POm_sem_modulation_cond2=sem_mod_vpmpom_PL(2,2);
SumStats(SumStatLine).POm_p_modulation_cond20=p_comp_mod_vs_0_vpmpom_PL(2,2);
SumStats(SumStatLine).POm_p_modulation_cond12=p_comp_mod_acrossPL_vpmpom(2);
SumStats(SumStatLine).VPMPOm_p_modulation_cond1=p_comp_mod_acrossvpmpom_PL(1);
SumStats(SumStatLine).VPMPOm_p_modulation_cond2=p_comp_mod_acrossvpmpom_PL(2);
SumStats(SumStatLine).VPMPOm_p_cond1_base=p_vpmpom_PLPnL_baseline_resp(1,1);
SumStats(SumStatLine).VPMPOm_p_cond2_base=p_vpmpom_PLPnL_baseline_resp(2,1);
SumStats(SumStatLine).VPMPOm_p_cond1_resp=p_vpmpom_PLPnL_baseline_resp(1,2);
SumStats(SumStatLine).VPMPOm_p_cond2_resp=p_vpmpom_PLPnL_baseline_resp(2,2);
SumStats(SumStatLine).VPM_mean_fsl_cond1=mean_fsl_vpmpom_PL(1,1);
SumStats(SumStatLine).VPM_sem_fsl_cond1=sem_fsl_vpmpom_PL(1,1);
SumStats(SumStatLine).VPM_mean_fsl_cond2=mean_fsl_vpmpom_PL(1,2);
SumStats(SumStatLine).VPM_sem_fsl_cond2=sem_fsl_vpmpom_PL(1,2);
SumStats(SumStatLine).VPM_p_fsl_cond12=p_comp_fsl_acrosscond_vpmpom_PL(1);
SumStats(SumStatLine).POm_mean_fsl_cond1=mean_fsl_vpmpom_PL(2,1);
SumStats(SumStatLine).POm_sem_fsl_cond1=sem_fsl_vpmpom_PL(2,1);
SumStats(SumStatLine).POm_mean_fsl_cond2=mean_fsl_vpmpom_PL(2,2);
SumStats(SumStatLine).POm_sem_fsl_cond2=sem_fsl_vpmpom_PL(2,2);
SumStats(SumStatLine).POm_p_fsl_cond12=p_comp_fsl_acrosscond_vpmpom_PL(2);
SumStats(SumStatLine).VPMPOm_p_fsl_cond1=p_comp_fsl_acrossvpmpom_PL(1);
SumStats(SumStatLine).VPMPOm_p_fsl_cond2=p_comp_fsl_acrossvpmpom_PL(2);

SumStatTable{1,SumTableCol}=SumStats(SumStatLine).VPM_N;
SumStatTable{2,SumTableCol}=SumStats(SumStatLine).VPM_N_animals;
SumStatTable{17,SumTableCol}=SumStats(SumStatLine).POm_N;
SumStatTable{18,SumTableCol}=SumStats(SumStatLine).POm_N_animals;
SumStatTable{33,SumTableCol}=SumStats(SumStatLine).Total_N_animals;
SumStatTable{3,SumTableCol}=SumStats(SumStatLine).VPM_mean_cond1_base;
SumStatTable{4,SumTableCol}=SumStats(SumStatLine).VPM_sem_cond1_base;
SumStatTable{6,SumTableCol}=SumStats(SumStatLine).VPM_mean_cond1_resp;
SumStatTable{7,SumTableCol}=SumStats(SumStatLine).VPM_sem_cond1_resp;
SumStatTable{9,SumTableCol}=SumStats(SumStatLine).VPM_p_cond1;
SumStatTable{3,SumTableCol+1}=SumStats(SumStatLine).VPM_mean_cond2_base;
SumStatTable{4,SumTableCol+1}=SumStats(SumStatLine).VPM_sem_cond2_base;
SumStatTable{6,SumTableCol+1}=SumStats(SumStatLine).VPM_mean_cond2_resp;
SumStatTable{7,SumTableCol+1}=SumStats(SumStatLine).VPM_sem_cond2_resp;
SumStatTable{9,SumTableCol+1}=SumStats(SumStatLine).VPM_p_cond2;
SumStatTable{5,SumTableCol}=SumStats(SumStatLine).VPM_p_base_cond12;
SumStatTable{8,SumTableCol}=SumStats(SumStatLine).VPM_p_resp_cond12;
SumStatTable{16+3,SumTableCol}=SumStats(SumStatLine).POm_mean_cond1_base;
SumStatTable{16+4,SumTableCol}=SumStats(SumStatLine).POm_sem_cond1_base;
SumStatTable{16+6,SumTableCol}=SumStats(SumStatLine).POm_mean_cond1_resp;
SumStatTable{16+7,SumTableCol}=SumStats(SumStatLine).POm_sem_cond1_resp;
SumStatTable{16+9,SumTableCol}=SumStats(SumStatLine).POm_p_cond1;
SumStatTable{16+3,SumTableCol+1}=SumStats(SumStatLine).POm_mean_cond2_base;
SumStatTable{16+4,SumTableCol+1}=SumStats(SumStatLine).POm_sem_cond2_base;
SumStatTable{16+6,SumTableCol+1}=SumStats(SumStatLine).POm_mean_cond2_resp;
SumStatTable{16+7,SumTableCol+1}=SumStats(SumStatLine).POm_sem_cond2_resp;
SumStatTable{16+9,SumTableCol+1}=SumStats(SumStatLine).POm_p_cond2;
SumStatTable{16+5,SumTableCol}=SumStats(SumStatLine).POm_p_base_cond12;
SumStatTable{16+8,SumTableCol}=SumStats(SumStatLine).POm_p_resp_cond12;
SumStatTable{10,SumTableCol}=SumStats(SumStatLine).VPM_mean_modulation_cond1;
SumStatTable{11,SumTableCol}=SumStats(SumStatLine).VPM_sem_modulation_cond1;
SumStatTable{12,SumTableCol}=SumStats(SumStatLine).VPM_p_modulation_cond10;
SumStatTable{10,SumTableCol+1}=SumStats(SumStatLine).VPM_mean_modulation_cond2;
SumStatTable{11,SumTableCol+1}=SumStats(SumStatLine).VPM_sem_modulation_cond2;
SumStatTable{12,SumTableCol+1}=SumStats(SumStatLine).VPM_p_modulation_cond20;
SumStatTable{13,SumTableCol}=SumStats(SumStatLine).VPM_p_modulation_cond12;
SumStatTable{16+10,SumTableCol}=SumStats(SumStatLine).POm_mean_modulation_cond1;
SumStatTable{16+11,SumTableCol}=SumStats(SumStatLine).POm_sem_modulation_cond1;
SumStatTable{16+12,SumTableCol}=SumStats(SumStatLine).POm_p_modulation_cond10;
SumStatTable{16+10,SumTableCol+1}=SumStats(SumStatLine).POm_mean_modulation_cond2;
SumStatTable{16+11,SumTableCol+1}=SumStats(SumStatLine).POm_sem_modulation_cond2;
SumStatTable{16+12,SumTableCol+1}=SumStats(SumStatLine).POm_p_modulation_cond20;
SumStatTable{16+13,SumTableCol}=SumStats(SumStatLine).POm_p_modulation_cond12;
SumStatTable{36,SumTableCol}=SumStats(SumStatLine).VPMPOm_p_modulation_cond1;
SumStatTable{36,SumTableCol+1}=SumStats(SumStatLine).VPMPOm_p_modulation_cond2;
SumStatTable{34,SumTableCol}=SumStats(SumStatLine).VPMPOm_p_cond1_base;
SumStatTable{34,SumTableCol+1}=SumStats(SumStatLine).VPMPOm_p_cond2_base;
SumStatTable{35,SumTableCol}=SumStats(SumStatLine).VPMPOm_p_cond1_resp;
SumStatTable{35,SumTableCol+1}=SumStats(SumStatLine).VPMPOm_p_cond2_resp;
SumStatTable{14,SumTableCol}=SumStats(SumStatLine).VPM_mean_fsl_cond1;
SumStatTable{15,SumTableCol}=SumStats(SumStatLine).VPM_sem_fsl_cond1;
SumStatTable{14,SumTableCol+1}=SumStats(SumStatLine).VPM_mean_fsl_cond2;
SumStatTable{15,SumTableCol+1}=SumStats(SumStatLine).VPM_sem_fsl_cond2;
SumStatTable{16,SumTableCol}=SumStats(SumStatLine).VPM_p_fsl_cond12;
SumStatTable{16+14,SumTableCol}=SumStats(SumStatLine).POm_mean_fsl_cond1;
SumStatTable{16+15,SumTableCol}=SumStats(SumStatLine).POm_sem_fsl_cond1;
SumStatTable{16+14,SumTableCol+1}=SumStats(SumStatLine).POm_mean_fsl_cond2;
SumStatTable{16+15,SumTableCol+1}=SumStats(SumStatLine).POm_sem_fsl_cond2;
SumStatTable{16+16,SumTableCol}=SumStats(SumStatLine).POm_p_fsl_cond12;
SumStatTable{37,SumTableCol}=SumStats(SumStatLine).VPMPOm_p_fsl_cond1;
SumStatTable{37,SumTableCol+1}=SumStats(SumStatLine).VPMPOm_p_fsl_cond2;

writetable(SumStatTable,[DataPath 'sumstats_' datestr(now,'yymmdd') '.xlsx'])
save([DataPath 'sumstats_' datestr(now,'yymmdd') '.mat'],'SumStats','params_*')
%%


