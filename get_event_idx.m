function [trig,trig_length]=get_event_idx(parameter, DiscreteData,triglength_limit,ppms,excludeNames,safety_margin,trigstartend,include)

DDcell=struct2cell(DiscreteData);
DDnames=fieldnames(DiscreteData);
trig=squeeze(DDcell(strcmp([parameter 'Start'],DDnames),1,:));
trig_length=squeeze(DDcell(strcmp([parameter 'Length'],DDnames),1,:));
%convert lenth 2 ind
triglength_limit=triglength_limit*ppms;
trigL=cell(size(trig));
for R=1:size(DiscreteData,2)
    t=trig{R};
    tl=trig_length{R};
    if strcmp('trigend',trigstartend)
        t=t+tl;
    end
    t=t(tl>=triglength_limit(1) & tl<=triglength_limit(2)); % select appropriate trigs
    trig{R}=t;
    trigL{R}=tl(tl>=triglength_limit(1) & tl<=triglength_limit(2));
end

trig1=trig;
for R=1:size(DiscreteData,2)
    [trig{R}]=exclude_regular_puff_repeats(trig{R});
end

trig=cleanup_sections(DiscreteData,excludeNames,safety_margin*ppms,'exclude',trig);

if ~isempty(include)
    [trig]=cleanup_sections(DiscreteData,include{1},include{2}*ppms,'forceinclude',trig);
end

if nargout>1
    for R=1:size(DiscreteData,2)
        td=trig{R}(:)-trig1{R}(:)';
        [a,b]=find(td==0);
        trig_length{R}=trig_length{R}(b);
    end
end
end


function [varargout]=cleanup_sections(DD,selectsections,safety,exorinclude,varargin)

varargout=cell(size(varargin));
DDcell=struct2cell(DD);
DDnames=fieldnames(DD);

for a=1:numel(varargin)
    E=varargin{a};
    for b=1:numel(E)
        Ex=zeros(1,DD(b).LengthInd);
        for N=1:numel(selectsections)
            if any(strcmp([selectsections{N} 'Start'],DDnames))
                D=DDcell{strcmp([selectsections{N} 'Start'],DDnames),1,b};
                Dl=DDcell{strcmp([selectsections{N} 'Length'],DDnames),1,b};
                for t=1:numel(D)
                    Ex(D(t):D(t)+Dl(t))=1;
                end
            else
                Stemp=DD(b).Sections.(selectsections{N});
                for t=1:size(Stemp,1)
                    Ex(Stemp(t,1):Stemp(t,2))=1;
                end
            end
        end
        trigs=E{b}(:)+[-safety(1):safety(end)];
        trigs(trigs<1)=1;
        trigs(trigs>DD(b).LengthInd)=DD(b).LengthInd;
        if strcmp(exorinclude,'exclude')
            varargout{a}{b}=E{b}(~any(Ex(trigs)>0,2));
        elseif strcmp(exorinclude,'forceinclude')
            varargout{a}{b}=E{b}(all(Ex(trigs)>0,2));
        end
    end
    varargout{a}=varargout{a}';
end


end