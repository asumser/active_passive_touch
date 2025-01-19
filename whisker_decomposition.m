%% compute whisker variables
DataPath='D:\active_passive_touch_data\data\';
bandpass_freqs=[4 25];%cutoff frequencies for bandpass
fs_filt=100;%sampling rate for bandpass filtering
fs_raw=1000;%sampling rate of whisker angle data
derivative_smoothing=16;%samples=ms if fs_raw=1000
slow_variable_smoothing=25;%samples=ms if fs_raw=1000
ephys_ppms=20;%1/1000 of sampling rate for ephys
amp_thresh=3;%minimum amplitude to count as whisking
min_whisk_duration=200;%samples=ms if fs_raw=1000

load([DataPath 'WhiskerBehavData.mat'],'Angle')
load([DataPath 'EphysData.mat'])
Phase=cell(size(Angle));
Velocity=cell(size(Angle));
Acceleration=cell(size(Angle));
Amp=cell(size(Angle));
Setp=cell(size(Angle));
wb=waitbar(0,'whisker parameter decomposition');
for RecN=1:numel(Angle)
    if ~isempty(Angle{RecN})
        Phase{RecN}=angle(hilbert(whisker_bpf(bandpass_freqs,Angle{RecN},fs_raw,fs_filt)));
        Velocity{RecN}=smooth(diff(smooth(Angle{RecN},derivative_smoothing)),derivative_smoothing)';
        Velocity{RecN}=cat(2,zeros(1,derivative_smoothing+1),Velocity{RecN}(derivative_smoothing+1:end-derivative_smoothing),zeros(1,derivative_smoothing));%clean edges
        Acceleration{RecN}=smooth(diff(smooth(Angle{RecN},derivative_smoothing),2),derivative_smoothing)';
        Acceleration{RecN}=cat(2,zeros(1,derivative_smoothing+2),Acceleration{RecN}(derivative_smoothing+1:end-derivative_smoothing),zeros(1,derivative_smoothing));%clean edges
        retractionIdx=find(diff(Phase{RecN})<-pi)+1;
        protractionIdx=find(diff(Phase{RecN}>0)==1);
        upper_envelope=interp1(protractionIdx,Angle{RecN}(protractionIdx),1:numel(Angle{RecN}),'linear');
        lower_envelope=interp1(retractionIdx,Angle{RecN}(retractionIdx),1:numel(Angle{RecN}),'linear');
        amp=movmedian(upper_envelope-lower_envelope,slow_variable_smoothing);
        amp(amp<0)=0;
        Amp{RecN}=amp;
        Setp{RecN}=movmedian((upper_envelope+lower_envelope)/2,slow_variable_smoothing);
    end
    waitbar(RecN/numel(Angle),wb)
end
close(wb)


DiscreteData=update_whisking_times(DiscreteData,Amp,ephys_ppms*1000/fs_raw,amp_thresh,min_whisk_duration);
DiscreteData=mk_nowhisk_sections(DiscreteData);

save([DataPath 'EphysData.mat'],"DiscreteData","MouseInfo","RecordingInfo","-append")
save([DataPath 'WhiskerBehavData.mat'],"Amp","Setp","Velocity","Acceleration","Phase","-append")

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [filtered_sig]=whisker_bpf(bp,signal,Fs_raw,Fs_filt)
signal_ds=resample(signal,Fs_filt,Fs_raw);
bp = bp * 2 / Fs_filt; % convert Hz to radians/S
[N, Wn] = buttord( bp, bp .* [.5 1.5], 3, 20);
[B,A] = butter(N,Wn);
filtered_sig_ds = filtfilt(B,A,signal_ds);
filtered_sig=interp1(1:numel(filtered_sig_ds),filtered_sig_ds,linspace(1,numel(signal_ds),numel(signal)),'spline');
end

%%

function DiscreteData=mk_nowhisk_sections(DiscreteData)

for n=1:numel(DiscreteData)
    W_start=DiscreteData(n).WhiskingStart;
    if ~isempty(W_start)
        W_end=DiscreteData(n).WhiskingStart+DiscreteData(n).WhiskingLength;
        NW_start=W_end+20;
        NW_end=W_start-20;
        if W_start(1)<=20
            NW_end(1)=[];
        else
            NW_start=[1 NW_start];
        end
        if W_end(end)>=DiscreteData(n).LengthInd-20
            NW_start(end)=[];
        else
            NW_end=[NW_end DiscreteData(n).LengthInd];
        end


        DiscreteData(n).NoWhiskingStart=NW_start;
        DiscreteData(n).NoWhiskingLength=NW_end-NW_start;
    else
        DiscreteData(n).NoWhiskingStart=[];
        DiscreteData(n).NoWhiskingLength=[];
    end

end
end
%%

function DiscreteData=update_whisking_times(DiscreteData,Amp,res,amp_thresh,min_dur)

for n=1:numel(Amp)
camp=Amp{n};
camp=camp>=amp_thresh;
wstarts=find(diff([0 camp])==1);
wends=find(diff([camp 0])==-1);
wdur=wends-wstarts;
wincl=wdur>=min_dur;
wstarts=wstarts(wincl)*res;
wdur=wdur(wincl)*res;
DiscreteData(n).WhiskingStart=wstarts;
DiscreteData(n).WhiskingLength=wdur;

end

end
