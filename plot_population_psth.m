function [popPSTHs,ax]=plot_population_psth(trig,DiscreteData,ppms,time_before_trig,time_after_trig,binsize,smoothwin,Name,response_types,selnames,sortby,scalem,varargin)


[popPSTHs,~,~,bins]=get_PSTHs(trig,DiscreteData,ppms,time_before_trig,time_after_trig,binsize,[],'none');
popPSTHs=cat(2,popPSTHs{:})';
pbins=bins(1:end-1)+diff(bins);
popPSTHs_smoothed=smoothdata(popPSTHs,2,'gaussian',milliseconds(smoothwin),'SamplePoints',milliseconds(pbins));




tiledlayout(5,numel(response_types))
m=0;

for s=1:numel(response_types)
    selectedPSTHs_smoothed=popPSTHs_smoothed(response_types{s},:)';
    cellnotnan=sum(isnan(selectedPSTHs_smoothed))==0;
    
    switch sortby
        case 'none'
            sortval=cellnotnan;
        case 'latency'
            [~,peakTimeSel]=max(movmean(selectedPSTHs_smoothed,5,1),[],1);
            [~,sortval]=sort(peakTimeSel);
            cellnotnan=cellnotnan(sortval);
            sortval=sortval(cellnotnan);
    end
   
    selectedPSTHs_smoothed_log_norm=max(log2(max(selectedPSTHs_smoothed-mean(selectedPSTHs_smoothed(pbins<0,:),'omitnan'),0)),-inf);
    ax(1,s)=nexttile(s,[2 1]);
    imagesc(pbins,1:sum(cellnotnan),selectedPSTHs_smoothed_log_norm(:,sortval)',[0 10]);    colormap(gray);
    cb(s)=colorbar;cb(s).Label.String='dRate';
    cb(s).TickLabels=cellstr(num2str(cb(s).Ticks(:).^2));
    title(selnames{s});box off;set(gca,'YTick',[]);
    m=max(m,max(mean(selectedPSTHs_smoothed,2,'omitnan')+std(selectedPSTHs_smoothed,0,2,'omitnan')/sqrt(sum(~isnan(selectedPSTHs_smoothed(1,:))))));
    ax(2,s)=nexttile(s+3*numel(response_types),[2 1]);
    
    bar(pbins,mean(selectedPSTHs_smoothed,2,'omitnan'),1,'facecolor',[.5 .5 .5],'edgecolor','none'); hold on
    plot(pbins,mean(selectedPSTHs_smoothed,2,'omitnan')+std(selectedPSTHs_smoothed,0,2,'omitnan')/sqrt(sum(~isnan(selectedPSTHs_smoothed(1,:)))),'k--');
    plot(pbins,mean(selectedPSTHs_smoothed,2,'omitnan')-std(selectedPSTHs_smoothed,0,2,'omitnan')/sqrt(sum(~isnan(selectedPSTHs_smoothed(1,:)))),'k--');
    set(gca,'xminortick','on','box','off')
    xlim(pbins([1 end]))
    ylabel('rate (Hz)')
    xlabel(['time (s)'])
    if size(varargin,1)~=0
        ax(3,s)=nexttile(s+2*numel(response_types));
        a=1;
        ccont=cat(1,varargin{1}{a}{response_types{s}});
        if strcmp(varargin{2}{a},'meanbefore')
            ccont=ccont-mean(ccont(:,pbins<0),2);
        end
        ccont_mean=mean(ccont,1,'omitnan');
        
        plot(pbins,ccont_mean,'r');box off
        xlim(pbins([1 end]))
        scaley{s}=ylim;
        ylabel(varargin{3})
        if numel(varargin{1})>1
            yyaxis right
            for a=2:numel(varargin{1})
                ccont=cat(1,varargin{1}{a}{response_types{s}});
                if strcmp(varargin{2}{a},'meanbefore')
                    ccont=ccont-mean(ccont(:,pbins<0),2);
                end
                ccont_mean=mean(ccont,1,'omitnan');
                if strcmp(varargin{2}{a},'scaleminmax')
                    h(s)=nexttile(s+3*numel(response_types),[2 1]);
                    ccont_mean=ccont_mean-min(ccont_mean);
                    ccont_mean=ccont_mean./max(ccont_mean);
                    ccont_mean=normalize(ccont_mean,'range',[0 m]);
                elseif contains(varargin{2}{a},'scale')
                    w(s)=nexttile(s+2*numel(response_types));
                    scale=str2double(varargin{2}{a}(6:end));
                    ccont_mean=ccont_mean*scale;
                elseif contains(varargin{2}{a},'mono')
                    w(s)=nexttile(s+2*numel(response_types));
                    ccont_mean(:)=mean(ccont_mean);
                    fprintf('var %u mean: %.2f\n',a,ccont_mean(1));
                end
                plot(pbins,ccont_mean,'b-');box off
                xlim(pbins([1 end]))
            end
        end
        
    end
    title([selnames{s} ' Pop ' Name ' Responses (N=' num2str(sum(~isnan(selectedPSTHs_smoothed(1,:)))) ')'])
    
end
scaley2=cat(1,scaley{:});
scaley3=[min(scaley2(:,1)) max(scaley2(:,2))];
if m<scalem
    m=scalem;
end

for s=1:numel(response_types)
    if size(varargin,1)~=0
        %
        set(h(s),'YLim',[0 m]);
        yyaxis (w(s),'left')
        set(w(s),'YLim',scaley3);
        yyaxis (w(s),'right')
        temp=get(w(s),'YLim');
        temp(1)=0;
        set(w(s),'YLim',temp);
    else
        set(h(s),'YLim',[0 m]);
    end
end
