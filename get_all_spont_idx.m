function [spont_idx]=get_all_spont_idx(DD,exclude,safetyEx,safetyIn,varargin)
spont_idx=cell(size(DD))'; %init
%%%%%%%%%%%%%%%%%%%%%
% decode exclude input
DDcell=struct2cell(DD);
DDnames=fieldnames(DD);
Snames=fieldnames(DD(1).Sections);
S=zeros(size(Snames));
D=zeros(size(DDcell,1),1);
Dl=zeros(size(DDcell,1),1);
% 
% for b=1:numel(DD)
%     
%     Ex=zeros(1,DD(b).LengthInd);
%     for N=1:numel(exclude)
%         if any(strcmp([exclude{N} 'Start'],DDnames))
%             D=DDcell{strcmp([exclude{N} 'Start'],DDnames),1,b};
%             Dl=DDcell{strcmp([exclude{N} 'Length'],DDnames),1,b};
%             for t=1:numel(D)
%                 Ex(D(t):D(t)+Dl(t))=1;
%             end
%         else
%             Stemp=DD(b).Sections.(exclude{N});
%             for t=1:size(Stemp,1)
%                 Ex(Stemp(t,1):Stemp(t,2))=1;
%             end
%         end
%     end
%     trigs=E{b}(:)+[-safety(1):safety(end)];
%     trigs(trigs<1)=1;
%     trigs(trigs>DD(b).LengthInd)=DD(b).LengthInd;
%     varargout{a}{b}=E{b}(~any(Ex(trigs)>0,2));
% end
% varargout{a}=varargout{a}';
% end

for N=1:numel(exclude)
    if any(strcmpi([exclude{N} 'Start'],DDnames)) %event like exclude parameters
        D=D+strcmpi([exclude{N} 'Start'],DDnames);
        Dl=Dl+strcmpi([exclude{N} 'Length'],DDnames);
    else
        S=S+strcmpi(exclude{N},Snames); %section exclude parameters
    end
end
if numel(varargin)>0 && ~isempty(varargin{1})
    Si=zeros(size(Snames));
    Di=zeros(size(DDcell,1),1);
    Dli=zeros(size(DDcell,1),1);
    include=varargin{1};
    for N=1:numel(include)
        if any(strcmpi([include{N} 'Start'],DDnames)) %event like exclude parameters
            Di=Di+strcmpi([include{N} 'Start'],DDnames);
            Dli=Dli+strcmpi([include{N} 'Length'],DDnames);
        else
            Si=Si+strcmpi(include{N},Snames); %section exclude parameters
        end
    end
end
%%%%%%%%%%%%%%%%%
%find random indices

for b=1:numel(DD) % all cells
    
    indexin1=[cell2mat(DDcell(find(D),1,b)')' cell2mat(DDcell(find(D),1,b)')'+cell2mat(DDcell(find(Dl),1,b)')']; %start/end indices of event parameters
    Sections=struct2cell(DD(b).Sections);
    indexin2=cell2mat(Sections(find(S)));%start/end indices of section parameters
    indexin=[indexin1; indexin2];%alltogether (double entries do not matter)
    if numel(varargin)>0 && ~isempty(varargin{1})
        indexout1=[cell2mat(DDcell(find(Di),1,b)')' cell2mat(DDcell(find(Di),1,b)')'+cell2mat(DDcell(find(Dli),1,b)')']; %start/end indices of event parameters
        Sections=struct2cell(DD(b).Sections);
        indexout2=cell2mat(Sections(find(Si)));%start/end indices of section parameters
        indexout=[indexout1; indexout2];%alltogether (double entries do not matter)
        if ~isempty(indexout)
        indexout=indexout-[-safetyIn(1) safetyIn(end)]; %expand with safety window
        end
        inidx=[];
        for i=1:size(indexout,1)
            inidx=cat(2,inidx,indexout(i,1):indexout(i,2));
        end
        inidx=unique(inidx);
    else
        inidx=1:DD(b).LengthInd;
    end
    if ~isempty(indexin)
        indexin=indexin+[-safetyEx(1) safetyEx(end)]; %expand with safety window
        
        outidx=[];
        for i=1:size(indexin,1)
         outidx=cat(2,outidx,indexin(i,1):indexin(i,2)); %remove trigs in undesirable sections
        end
         outidx=unique(outidx);
         outidx(outidx>DD(b).LengthInd)=[];
         outidx(outidx<=0)=[];
         %inidx(outidx)=[];
         spont_idx{b}=setdiff(inidx,outidx);
    else
        spont_idx{b}=inidx;%1:DD(b).LengthInd; %select all
    end
end