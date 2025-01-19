function [varargout]=select_trigs_by_state(DiscreteData,exclude,include,safetyEx,safetyIn,varargin)

varargout=cell(size(varargin));
DDcell=struct2cell(DiscreteData);
DDnames=fieldnames(DiscreteData);

for EventInd=1:numel(varargin)
    events=varargin{EventInd};%all events
    for RecInd=1:numel(events)
        % remove events in excluded states
        StateVector=zeros(1,DiscreteData(RecInd).LengthInd);
        for N=1:numel(exclude)
            if any(strcmp([exclude{N} 'Start'],DDnames))
                D=DDcell{strcmp([exclude{N} 'Start'],DDnames),1,RecInd};
                Dl=DDcell{strcmp([exclude{N} 'Length'],DDnames),1,RecInd};
                for t=1:numel(D)
                    StateVector(D(t):D(t)+Dl(t))=1;
                end
            else
                Stemp=DiscreteData(RecInd).Sections.(exclude{N});
                for t=1:size(Stemp,1)
                    StateVector(Stemp(t,1):Stemp(t,2))=1;
                end
            end
        end
        EventTriggeredSafetyInds=events{RecInd}(:)+[-safetyEx(1):safetyEx(end)];
        EventTriggeredSafetyInds(EventTriggeredSafetyInds<1)=1;
        EventTriggeredSafetyInds(EventTriggeredSafetyInds>DiscreteData(RecInd).LengthInd)=DiscreteData(RecInd).LengthInd;
        events{RecInd}=events{RecInd}(~any(StateVector(EventTriggeredSafetyInds)>0,2));
        
        if ~isempty(include)
            % include only events in included states
            StateVector=zeros(1,DiscreteData(RecInd).LengthInd);
            for N=1:numel(include)
                if any(strcmp([include{N} 'Start'],DDnames))
                    D=DDcell{strcmp([include{N} 'Start'],DDnames),1,RecInd};
                    Dl=DDcell{strcmp([include{N} 'Length'],DDnames),1,RecInd};
                    for t=1:numel(D)
                        StateVector(D(t):D(t)+Dl(t))=1;
                    end
                else
                    Stemp=DiscreteData(RecInd).Sections.(include{N});
                    for t=1:size(Stemp,1)
                        StateVector(Stemp(t,1):Stemp(t,2))=1;
                    end
                end
            end
            EventTriggeredSafetyInds=events{RecInd}(:)+[-safetyIn(1):safetyIn(end)];
            EventTriggeredSafetyInds(EventTriggeredSafetyInds<1)=1;
            EventTriggeredSafetyInds(EventTriggeredSafetyInds>DiscreteData(RecInd).LengthInd)=DiscreteData(RecInd).LengthInd;
            events{RecInd}=events{RecInd}(all(StateVector(EventTriggeredSafetyInds)>0,2));


        end
        varargout{EventInd}{RecInd}=events{RecInd};
    end
    varargout{EventInd}=varargout{EventInd}';
end
