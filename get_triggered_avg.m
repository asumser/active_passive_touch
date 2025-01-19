function [OverallAvg,AvgTime,EachTrial]=get_triggered_avg(trig,relativeTime,Name,DynamicVariableDS,DiscreteData)

if iscell(Name)
    DynamicVariable=Name;
    dvinput=true;
else
    dvinput=false;
end
EachTrial=cell(size(trig));
OverallAvg=cell(size(trig));
AvgTime=cell(size(trig));
for n=1:numel(trig)
    if ~isempty(trig{n})
        if dvinput && isempty(DynamicVariable{n});continue;end
        eT=linspace(0,DiscreteData(n).LengthTime,DiscreteData(n).LengthInd)';
        trigT=eT(trig{n});
        trigTest=trigT(:)+(relativeTime(:)/1000)';
        if DynamicVariableDS~=1
            DynamicVariableT=downsample(eT,DynamicVariableDS);
        else
            DynamicVariableT=eT;
        end
        if ~dvinput
            temp=zeros(size(eT));
            trigs=DiscreteData(n).([Name 'Start']);
            lengths=DiscreteData(n).([Name 'Length']);
            for t=1:numel(trigs)
                temp(trigs(t):trigs(t)+lengths(t))=1;
            end
        else
            temp=DynamicVariable{n};
        end
        EachTrial{n}=interp1(DynamicVariableT,temp,trigTest,'linear');
        OverallAvg{n}=mean(EachTrial{n},1,'omitnan');
        AvgTime{n}=mean(EachTrial{n},2,'omitnan');
    end
end

