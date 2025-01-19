function [Bursts,burstI]=returnBursts(sp,isiCutoff)



dim=size(sp); if dim(2)>dim(1),sp=sp';end
Bursts={};
if isempty(sp)
    Bursts={};
elseif numel(sp)==1
    Bursts{1}=sp;
else
    isis=[isiCutoff*2; diff(sp)];
    bursts=(find(isis>isiCutoff)); %indices of burst starts;
    burstI=nan(2,numel(bursts));
    for i=1:numel(bursts)
        if i<numel(bursts)
            i1=bursts(i);
            i2=bursts(i+1)-1;
            Bursts{i}=[sp(i1:i2)];
            burstI(2,i)=numel(Bursts{i});
        else
            i1=bursts(i);
            Bursts{i}=[sp(i1:end)];
            burstI(2,i)=numel(Bursts{i});
        end
    end
    burstI(1,:)=sp(bursts);
    
end


