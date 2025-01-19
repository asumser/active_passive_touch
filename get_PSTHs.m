function [H,p,ntrigs,bins]=get_PSTHs(trig,DiscreteData,ppms,time_before_trig,time_after_trig,binsize,baselinebinnum,tail)
ntrigs=cellfun(@numel,trig);
tempDD=struct2cell(DiscreteData);
spikes=squeeze(tempDD(1,1,:));clear temp*
bins=(-time_before_trig:binsize:time_after_trig);
if isempty(baselinebinnum)
    baselinebinnum=find(bins<=0,1,'last');
end
H=cell(size(trig));
p=cell(1,size(trig,2));
p(:)={nan(size(trig,1),numel(bins)-1)};
for R=1:size(trig,1)
    for S=1:size(trig,2)
        h=histcounts((reshape(spikes{R},1,[])-reshape(trig{R,S},[],1))/ppms,bins);
        H{R,S}=((1000./diff(bins)).*h./numel(trig{R,S}))';
        if strcmp(tail,'none')
            p{S}(R,:)=nan;
        else
            p{S}(R,:)=rate_stats(trig{R,S},spikes{R},bins,baselinebinnum,ppms,tail);
        end
    end
end
