function [Rates,Rateshuf_p]=state_rates(DiscreteData,exclude,include, outsafety,insafety,ppms,shuffN)
Rates=nan(size(DiscreteData,2),1);
Rateshuf_p=Rates;
for n=1:size(DiscreteData,2)
    exclude_vec=zeros(1,DiscreteData(n).LengthInd);
    for ex=1:numel(exclude)
        switch exclude{ex}
            case {'Puff','Touch','Whisking','Light','NoWhisking'}
                starts=DiscreteData(n).(strcat(exclude{ex},'Start'));
                lengths=DiscreteData(n).(strcat(exclude{ex},'Length'));
                ends=starts+lengths;
                starts=max(1,starts-outsafety*ppms);
                ends=min(ends+outsafety*ppms,DiscreteData(n).LengthInd);
                for nex=1:numel(starts)
                    exclude_vec(starts(nex):ends(nex))=1;
                end
            otherwise
                if ~isempty(DiscreteData(n).Sections.(exclude{ex}))
                    starts=DiscreteData(n).Sections.(exclude{ex})(:,1);
                    ends=DiscreteData(n).Sections.(exclude{ex})(:,1);
                    starts=max(1,starts-outsafety*ppms);
                    ends=min(ends+outsafety*ppms,DiscreteData(n).LengthInd);
                    for nex=1:numel(starts)
                        exclude_vec(starts(nex):ends(nex))=1;
                    end
                end
        end
    end
    exclude_vec=exclude_vec>0;

    include_vec=zeros(1,DiscreteData(n).LengthInd);
    if ~isempty(include)
        for in=1:numel(include)
            switch include{in}
                case {'Puff','Touch','Whisking','Light','NoWhisking'}
                    starts=DiscreteData(n).(strcat(include{in},'Start'));
                    lengths=DiscreteData(n).(strcat(include{in},'Length'));
                    ends=starts+lengths;
                    starts=max(1,starts-insafety*ppms);
                    ends=min(ends+insafety*ppms,DiscreteData(n).LengthInd);
                    for nex=1:numel(starts)
                        include_vec(starts(nex):ends(nex))=include_vec(starts(nex):ends(nex))+1;
                    end
                otherwise
                    starts=DiscreteData(n).Sections.(include{in})(:,1);
                    ends=DiscreteData(n).Sections.(include{in})(:,1);
                    starts=max(1,starts-insafety*ppms);
                    ends=min(ends+insafety*ppms,DiscreteData(n).LengthInd);
                    for nex=1:numel(starts)
                        include_vec(starts(nex):ends(nex))=include_vec(starts(nex):ends(nex))+1;
                    end
            end
        end
    end

    include_vec= include_vec>=numel(include);

    full_vec=include_vec & ~exclude_vec;
    full_time=sum(full_vec)/(ppms*1000);
    spikeT=full_vec(DiscreteData(n).Spikes);
    N_included_spikes=sum(spikeT);
    Rates(n)=N_included_spikes/full_time;

    if shuffN>0 && nargout>1
        include_vec_crop=include_vec(~exclude_vec);
        spikex=false(size(include_vec));
        spikex(DiscreteData(n).Spikes)=true;
        spikex=spikex(~exclude_vec);
        spikes_cropped=find(spikex);
        N_shuff_incl_spikes=zeros(shuffN,1);
        spikeInt=diff([0 spikes_cropped]);

        parfor s=1:shuffN
            shuffSpikeInt=spikeInt(randperm(numel(spikeInt)));
            shuffSpikes=cumsum(shuffSpikeInt);
            spikeT=include_vec_crop(shuffSpikes);
            N_shuff_incl_spikes(s)=sum(spikeT);
        end
        p=sum(N_included_spikes<N_shuff_incl_spikes)/shuffN;
        Rateshuf_p(n)=p;
    end

end




%%
% DDcell=struct2cell(DiscreteData);
% DDnames=fieldnames(DiscreteData);
% [exclude_sections]=select_sections(DiscreteData,exclude,outsafety);
% if ~isempty(include)
% [include_sections]=select_sections(DiscreteData,include,insafety);
% else
%      include_sections=cell(size(DiscreteData,2),1);
%     for n=1:size(DiscreteData,2); include_sections{n}=[1 DiscreteData(n).LengthInd]; end
% end
% [~,time_in_sections]=split_insections(include_sections,exclude_sections);
%
% [Spikes]=cleanup_sections_exclusive(DiscreteData,exclude,outsafety*ppms,squeeze(DDcell(1,:,:)));
% if ~isempty(include);[Spikes]=cleanup_sections_inclusive(DiscreteData,include,0,Spikes);end
% if nargout>1
% SinB=nan(size(Spikes));
% for n=1:numel(Spikes);[Bursts]=returnBursts2(Spikes{n},isi*ppms);SinB(n)=mean(cellfun(@numel,Bursts));end
% end
% N_Spikes=cellfun(@numel,Spikes);
% Rates=(N_Spikes./cell2mat(time_in_sections))*(ppms*1000);

