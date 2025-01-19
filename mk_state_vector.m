function [cont]=mk_state_vector(DiscreteData,Name)
cont=cell(numel(DiscreteData),1);
for n=1:numel(DiscreteData)
    temp=zeros(1,DiscreteData(n).LengthInd);
    trigs=DiscreteData(n).([Name 'Start']);
    lengths=DiscreteData(n).([Name 'Length']);
    for t=1:numel(trigs)
    temp(trigs(t):trigs(t)+lengths(t))=1;
    end
    cont{n}=temp;
end