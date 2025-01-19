function [not_regular_puffs]=exclude_regular_puff_repeats(trigP)

[Cpuff]=returnBursts(trigP,9000);
dCPuff=cellfun(@diff,Cpuff,'UniformOutput',0);
reg=false(size(dCPuff));
for n=1:numel(dCPuff)
    if numel(dCPuff{n})>3
        if ~any(dCPuff{n}(1:3)~=4400)
            reg(n)=1;
        end
    end
    dCPuff{n}=[nan ;dCPuff{n}];
end
not_regular_puffs=cell2mat(Cpuff(~reg)');

end



function [Bursts,burstI]=returnBursts(sp,intervalCutoff)



dim=size(sp); if dim(2)>dim(1),sp=sp';end
Bursts={};
if isempty(sp)
    Bursts={};
elseif numel(sp)==1
    Bursts{1}=sp;
else
    isis=[intervalCutoff*2; diff(sp)];
    bursts=(find(isis>intervalCutoff)); %indices of burst starts;
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


end
