function md=comp_modulation(x,varargin) 

if numel(varargin)~=0
x=cat(2,x,varargin{1});
end
md=diff(x,1,2)./sum(x,2);
