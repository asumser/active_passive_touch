function p=rate_stats(trig,events,bins,baselinebinnum,ppms,tail)

   trig=trig(:);
    events=events(:)';
    T=repmat(trig,1,numel(events));
    E=repmat(events,numel(trig),1);
    M=E-T;
    M(M<=min(bins)*ppms | M>=max(bins)*ppms)=nan;
    N=discretize(M,bins*ppms);
    Nx=[];for n=1:numel(bins)-1;Nx=[Nx sum(N==n,2)];end
    Nx=1000*Nx./repmat(diff(bins),size(Nx,1),1);
    p=nan(1,numel(bins)-1);
    for n=1:numel(bins)-1
        if any(~isnan(Nx(:,n)))
            if strcmp(tail,'larger') || strcmp(tail,'smaller') || strcmp(tail,'unequal')
                [~,p(n)]=kstest2(Nx(:,baselinebinnum),Nx(:,n),'Tail',tail);
            elseif strcmp(tail,'right') || strcmp(tail,'left') || strcmp(tail,'both')
                [p(n)]=signrank(Nx(:,baselinebinnum),Nx(:,n),'Tail',tail);
            end
        end
    end
end
