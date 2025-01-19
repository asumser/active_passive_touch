function [RasterPT,trial_kinematics,paramOut]=split_by_kinematics(trigsP,trigsT,min_trial_per_quantile,split_time_before_trig,split_time_after_trig,split_binsize,smoothwidth,meth,binsR,...
    selc,splits,splitname,DiscreteData,ppms,split_var_nameC,type_nameC,split_allR,figdir,Amp,Phase,Angle,Setp,Velocity,Acceleration,Curvature,plotRatesAs)

PTCol=[49 94 91; 214 143 0]/255;
selN={'VPM';'POm'};
selS={'Puff';'Touch'};

bins=binsR(1):split_binsize:binsR(end);
pbins=bins(1:end-1);
pbins=pbins+mean(diff(pbins)/2);

paramOut={};
paramoutline=0;
trigsAll=cat(2,trigsP,trigsT);

RasterPT=cell(size(DiscreteData,2),numel(selS));
for cellN=1:size(DiscreteData,2)
    spikes=DiscreteData(cellN).Spikes;
    for S=1:numel(selS)
        [RasterPT{cellN,S}]=get_raster_single(trigsAll{cellN,S}, spikes,ppms,-binsR(1),binsR(end));
    end
end

trial_kinematics=cell(numel(selS),1);
trial_kinematics(:)={cell(size(trigsT,1),numel(split_var_nameC),numel(type_nameC))};

for SN=1:numel(split_var_nameC)
    split_var_name=split_var_nameC{SN};
    switch split_var_name
        case 'Interval'
            for n=1:size(RasterPT,1)
                for S=1:numel(selS)
                    if ~isempty(trigsAll{n,S})
                        trial_kinematics{S}{n,SN,1}=cat(1,trigsAll{n,S}(1),diff(trigsAll{n,S}))/(ppms*1000);
                    end
                end
            end

        case {'Angle','Setp','Phase','Curvature','Acceleration','Velocity','Amp'}
            switch split_var_name
                case {'Angle'}
                    behavior_variable=Angle;
                    behavior_variable=cellfun(@(x) x-median(x,'omitnan'),behavior_variable,'UniformOutput',0);

                case {'Setp'}
                    behavior_variable=Setp;
                    behavior_variable=cellfun(@(x) x-median(x,'omitnan'),behavior_variable,'UniformOutput',0);
                case {'Amp'}
                    behavior_variable=Amp;
                case {'Acceleration'}
                    behavior_variable=Acceleration;
                case {'Velocity'}
                    behavior_variable=Velocity;
                case 'Phase'
                    behavior_variable=Phase;
                    for x=1:numel(behavior_variable)
                        behavior_variable{x}(Amp{x}<3)=nan;
                    end
                case 'Curvature'
                    behavior_variable=cellfun(@(x,y) x-median(x(y<3),'omitnan'),Curvature,Amp,'UniformOutput',0);
                    behavior_variable=cellfun(@(x) x./std(x,0,'omitnan'),behavior_variable,'UniformOutput',0);
            end

            kinematic_timerange=pbins(pbins<split_time_after_trig & pbins>-split_time_before_trig);
            % [~,temp_P]=get_trig_cont_pop(trigsP,kinematic_timerange,behavior_variable,20,DiscreteData);
            triggered_behavior_variable=cell(size(RasterPT,1),2);
            for S=1:2
                [~,~,triggered_behavior_variable(:,S)]=get_triggered_avg(trigsAll(:,S),kinematic_timerange,behavior_variable,20,DiscreteData);
            end
            removeData=cellfun(@isempty,triggered_behavior_variable);
            triggered_behavior_variable(removeData)={nan(size(kinematic_timerange))};

            for tx=1:numel(type_nameC)
                %{'rawbefore','rawafter','absbefore','absafter','rel','absrel'};
                for S=1:2
                    switch type_nameC{tx}
                        case 'rawbefore'
                            trial_kinematics{S}(:,SN,tx)=cellfun(@(x) mean(x(:,kinematic_timerange<0),2,'omitnan'),triggered_behavior_variable(:,S),'UniformOutput',0);
                        case 'rawafter'
                            trial_kinematics{S}(:,SN,tx)=cellfun(@(x) mean(x(:,kinematic_timerange>0),2,'omitnan'),triggered_behavior_variable(:,S),'UniformOutput',0);
                        case 'absbefore'
                            trial_kinematics{S}(:,SN,tx)=cellfun(@(x) abs(mean(x(:,kinematic_timerange<0),2,'omitnan')),triggered_behavior_variable(:,S),'UniformOutput',0);
                        case 'absafter'
                            trial_kinematics{S}(:,SN,tx)=cellfun(@(x) abs(mean(x(:,kinematic_timerange>0),2,'omitnan')),triggered_behavior_variable(:,S),'UniformOutput',0);
                        case 'rel'
                            trial_kinematics{S}(:,SN,tx)=cellfun(@(x) mean(x(:,kinematic_timerange>0),2,'omitnan')-mean(x(:,kinematic_timerange<0),2,'omitnan'),triggered_behavior_variable(:,S),'UniformOutput',0);
                        case 'absrel'
                            trial_kinematics{S}(:,SN,tx)=cellfun(@(x) abs(mean(x(:,kinematic_timerange>0),2,'omitnan')-mean(x(:,kinematic_timerange<0),2,'omitnan')),triggered_behavior_variable(:,S),'UniformOutput',0);
                    end
                end
            end
            for S=1:2
                trial_kinematics{S}(removeData(:,S),SN,:)={[]};
            end
    end
end


%% plot by kinematics
% split_var_nameC={'Acceleration';'Curvature';'Interval';'Velocity';'Setp';'Amp';'Angle';'Phase'};
% type_nameC={'raw','abs','rel','absrel'};
S=1;
for sa=1:numel(split_allR)
    split_all=split_allR(sa);
    if split_all

        js='joint_cells';
    else
        js='indiv_cells';
    end
    for SN=1:numel(split_var_nameC)
        split_var_name=split_var_nameC{SN};
        for TY=1:size(trial_kinematics{1},3)
            figadd=strcat(split_var_nameC{SN},'_',type_nameC{TY},'_',js);
            temp1=cat(1,cat(1,trial_kinematics{1}{:,SN,TY}),cat(1,trial_kinematics{2}{:,SN,TY}));
            temp1(isnan(temp1))=[];
            if isempty(temp1);continue;end
            paramoutline=paramoutline+1;
            paramOut{paramoutline,1}=split_var_name;
            paramOut{paramoutline,2}=type_nameC{TY};

            Ilim=prctile(temp1,splits);
            paramOut{paramoutline,3}=Ilim(2:end-1);
            ix1=discretize(temp1,Ilim);
            mean_value_bins=nan(max(ix1),1);
            for j=1:max(ix1)
                mean_value_bins(j)=mean(temp1(ix1==j));
            end
            labels=cell(size(Ilim,2)-1,1);
            for l=1:numel(Ilim)-2
                labels{l}=sprintf('v<= %.3f',Ilim(l+1));
            end
            labels{end}=sprintf('v> %.3f',Ilim(end-1));
            splitPSTH=nan(numel(trigsT),numel(Ilim)-1,2,numel(bins)-1);
            splitRates=nan(numel(trigsT),numel(Ilim)-1,2,numel(binsR)-1);
            Climx=nan(numel(trigsT),numel(Ilim));
            for n=1:numel(trigsT)
                if split_all
                    Clim=Ilim;
                else
                    Clim=prctile(cat(1,trial_kinematics{1}{n,SN,TY},trial_kinematics{2}{n,SN,TY}),splits);
                end
                Climx(n,:)=Clim;
                hx=nan(2,numel(Ilim)-1);
                for S=1:2
                    if ~isempty(RasterPT{n,S}) && ~isempty(trial_kinematics{S}{n,SN,TY})
                        ix=discretize(trial_kinematics{S}{n,SN,TY},Ilim)';
                        sb=ix(RasterPT{n,S}(:,2));
                        hx(S,:)=histcounts(sb,.5:numel(Ilim)-.5);
                    end
                end
                if all(hx>=min_trial_per_quantile,"all")
                    for S=1:2
                        if ~isempty(RasterPT{n,S}) && ~isempty(trial_kinematics{S}{n,SN,TY})
                            ix=discretize(trial_kinematics{S}{n,SN,TY},Ilim)';
                            sb=ix(RasterPT{n,S}(:,2));
                            for x=1:numel(Ilim)-1
                                if (min_trial_per_quantile<0 && hx(S,x)>=-min_trial_per_quantile) || min_trial_per_quantile>=0
                                splitPSTH(n,x,S,:)=histcounts(RasterPT{n,S}(sb==x,1),bins)./sum(ix==x);
                                splitRates(n,x,S,:)=histcounts(RasterPT{n,S}(sb==x,1),binsR)./sum(ix==x)./(diff(binsR)/1000);
                                end
                            end
                        end

                    end
                end
            end

            C=repmat(linspace(.7,0,numel(Ilim)-1)',1,3);
            fig1= figure('Position',[50   701   783   580]);
            tx=0;
            ML=0;

            tt=tiledlayout(4,6);

            title(tt,sprintf('PSTHs, split by %s of %s',splitname,split_var_name),sprintf('normalization: %s, split: %s',type_nameC{TY},js),'Interpreter','none')
            nexttile(1,[2 2]);
            switch split_var_name
                case  'Interval'
                    vbins=[linspace(Ilim(1),1.25,50) Ilim(end)];
                    pvbins=vbins(1:end-1)+median(diff(vbins))/2;
                case 'Phase'
                    vbins=linspace(-pi,pi,50);
                    pvbins=vbins(1:end-1)+median(diff(vbins))/2;
                otherwise
                    psplits=splits;
                    psplits(1)=1;psplits(end)=99;
                    Ilimp=prctile(temp1,psplits);
                    vbins=linspace(Ilimp(1),Ilimp(end),51);
                    pvbins=vbins(1:end-1)+median(diff(vbins))/2;

            end

            H1=histcounts(cat(1,trial_kinematics{1}{:,SN,TY}),vbins,'Normalization','probability');
            H2=histcounts(cat(1,trial_kinematics{2}{:,SN,TY}),vbins,'Normalization','probability');
            switch split_var_name
                case  'Phase'

                    gr=[220 220 220;144 144 144;91 91 91]/255;
                    polarhistogram('BinEdges',vbins,'BinCounts',H1,'EdgeColor','none','FaceColor',PTCol(1,:));hold on
                    polarhistogram('BinEdges',vbins,'BinCounts',H2,'EdgeColor','none','FaceColor',PTCol(2,:),'FaceAlpha',.7);
                    polarplot(linspace(-pi,Ilim(2),50),repmat(max(cat(2,H1,H2))*1.25,1,50),'LineWidth',4,'Color',gr(1,:))
                    polarplot(linspace(Ilim(2),Ilim(3),50),repmat(max(cat(2,H1,H2))*1.25,1,50),'LineWidth',4,'Color',gr(2,:))
                    polarplot(linspace(Ilim(3),pi,50),repmat(max(cat(2,H1,H2))*1.25,1,50),'LineWidth',4,'Color',gr(3,:))
                otherwise
                    bar(pvbins,H1,1,'FaceColor',PTCol(1,:),'EdgeColor','none');hold on
                    bar(pvbins,H2,1,'FaceColor',PTCol(2,:),'EdgeColor','none','FaceAlpha',.7);hold on
                    line(repmat(Ilim(2:end-1),2,1),[0 1.25*max([H1 H2])],'Color','k');
            end
            box off
            legend(selS);legend boxoff
            title(split_var_name)
            switch plotRatesAs
                case 'modulation'
                    splitRatesC=diff(splitRates,1,4)./sum(splitRates,4);
                case 'difference'
                    splitRatesC=diff(splitRates,1,4);
                case 'log2foldchange'
                    splitRatesC=log2(splitRates(:,:,:,2)./splitRates(:,:,:,1));
                case 'rawresponse'
                    splitRatesC=splitRates(:,:,:,2);
            end
            na=cat(1,selc{:});
            selectedSplitRates_all=splitRatesC(na,:,:);
            selectedSplitRates_all(isinf(selectedSplitRates_all))=nan;
            mean_selectedSplitRates_all=squeeze(mean(selectedSplitRates_all,1,'omitnan'));
            sem_selectedSplitRates_all=squeeze(std(selectedSplitRates_all,0,1,'omitnan')./sqrt(sum(~isnan(selectedSplitRates_all))));
            yla=cat(2,floor(10*min(0,min(mean_selectedSplitRates_all-sem_selectedSplitRates_all,[],'all')))/10,ceil(10*1.35*max(mean_selectedSplitRates_all+sem_selectedSplitRates_all,[],'all'))/10);
            for nx=1:2
                nt(nx)=nexttile(1+nx*2,[2 2]);
                n=selc{nx};
                selectedSplitRates=splitRatesC(n,:,:);
                selectedSplitRates(isinf(selectedSplitRates))=nan;
                mean_selectedSplitRates=squeeze(mean(selectedSplitRates,1,'omitnan'));
                sem_selectedSplitRates=squeeze(std(selectedSplitRates,0,1,'omitnan')./sqrt(sum(~isnan(selectedSplitRates))));
                p_selectedSplitRates_vs_0=nan(size(selectedSplitRates,3),size(selectedSplitRates,2));
                for pt=1:size(selectedSplitRates,3)
                    for tx=1:size(selectedSplitRates,2)
                        if any(~isnan(selectedSplitRates(:,tx,pt)))
                            p_selectedSplitRates_vs_0(pt,tx)=signrank(selectedSplitRates(:,tx,pt),[],'tail','right');
                        end
                    end
                end
                paramOut{paramoutline,4+(nx-1)*8}=mean_selectedSplitRates(:,1)';
                paramOut{paramoutline,5+(nx-1)*8}=sem_selectedSplitRates(:,1)';
                paramOut{paramoutline,6+(nx-1)*8}=p_selectedSplitRates_vs_0(1,:);
                paramOut{paramoutline,7+(nx-1)*8}=mean_selectedSplitRates(:,2)';
                paramOut{paramoutline,8+(nx-1)*8}=sem_selectedSplitRates(:,2)';
                paramOut{paramoutline,9+(nx-1)*8}=p_selectedSplitRates_vs_0(2,:);
                p_selectedSplitRates_vs_0=p_selectedSplitRates_vs_0(:);
                indicator_selectedSplitRates_vs_0=cell(size(p_selectedSplitRates_vs_0));
                for np=1:numel(p_selectedSplitRates_vs_0)
                    if p_selectedSplitRates_vs_0(np)<.001
                        indicator_selectedSplitRates_vs_0{np}='*';
                    elseif p_selectedSplitRates_vs_0(np)<.01
                        indicator_selectedSplitRates_vs_0{np}='*';
                    elseif    p_selectedSplitRates_vs_0(np)<.05
                        indicator_selectedSplitRates_vs_0{np}='*';
                    else
                        indicator_selectedSplitRates_vs_0{np}='';
                    end
                end

                p_selectedSplitRates_P_vs_T=nan(1,size(selectedSplitRates,2));

                for tx=1:size(selectedSplitRates,2)
                    if any(all(~isnan(selectedSplitRates(:,tx,[1 2])),3))
                        p_selectedSplitRates_P_vs_T(1,tx)=signrank(selectedSplitRates(:,tx,1),selectedSplitRates(:,tx,2));
                    end
                end
                paramOut{paramoutline,10+(nx-1)*8}=p_selectedSplitRates_P_vs_T(1,:);
                indicator_selectedSplitRates_P_vs_T=cell(size(p_selectedSplitRates_P_vs_T));
                for np=1:numel(p_selectedSplitRates_P_vs_T)
                    if p_selectedSplitRates_P_vs_T(np)<.001
                        indicator_selectedSplitRates_P_vs_T{np}='***';
                    elseif p_selectedSplitRates_P_vs_T(np)<.01
                        indicator_selectedSplitRates_P_vs_T{np}='**';
                    elseif    p_selectedSplitRates_P_vs_T(np)<.05
                        indicator_selectedSplitRates_P_vs_T{np}='*';
                    else
                        indicator_selectedSplitRates_P_vs_T{np}='';
                    end
                end


                b{nx}=bar(1:size(mean_selectedSplitRates,1),mean_selectedSplitRates);hold on;
                x = nan(size(mean_selectedSplitRates,2), size(mean_selectedSplitRates,1));
                for i = 1:size(mean_selectedSplitRates,2)
                    x(i,:) = b{nx}(i).XEndPoints;
                    b{nx}(i).FaceColor=PTCol(i,:);
                end

                errorbar(x',mean_selectedSplitRates,sem_selectedSplitRates,'k','linestyle','none');
                text(sort(x(:)),reshape(mean_selectedSplitRates'+sem_selectedSplitRates',[],1),indicator_selectedSplitRates_vs_0,'HorizontalAlignment','center','VerticalAlignment','bottom')

                line(x,repmat(1.15*max(mean_selectedSplitRates+sem_selectedSplitRates,[],'all'),size(x,1),size(x,2)),'color','k')
                text(mean(x,1),repmat(1.15*max(mean_selectedSplitRates+sem_selectedSplitRates,[],'all'),1,size(x,2)),indicator_selectedSplitRates_P_vs_T,'HorizontalAlignment','center','VerticalAlignment','bottom')
                yl=ylim;
                yl(2)=1.35*max(mean_selectedSplitRates+sem_selectedSplitRates,[],'all');
                if ~any(isnan(yl))
                ylim(yl);
                end
                hold off
                ylabel(plotRatesAs)
                xlabel(splitname)
                box off

            end

            tx=0;
            starttile=12;
            for s=1:2
                for nx=1:2
                    n=selc{nx};
                    tx=tx+1;
                    tind(tx)=nexttile(starttile+tx*3-2,[1 3]);
                    selectedPSTHs=permute(splitPSTH(n,:,s,:),[4 2 1 3]);
                    selectedPSTHs_mean=mean(selectedPSTHs,3,'omitnan');
                    selectedPSTHs_mean=smoothdata(selectedPSTHs_mean,1,meth,smoothwidth,'SamplePoints',milliseconds(pbins));
                    for bin=1:size(selectedPSTHs_mean,2)
                        plot(pbins',selectedPSTHs_mean(:,bin),'Color',C(bin,:));hold on;box off;
                    end
                    ML=max(ML,max(selectedPSTHs_mean,[],'all'));
                    xlabel('time around deflection (ms)')
                    ylabel({'Population mean';'spike prob. per bin'})
                    title(sprintf('%s (N=%u) %s',selN{nx},sum(any(~all(isnan(selectedPSTHs),1),2)),selS{s}));
                end
            end
            if split_all
                legend(labels);legend('boxoff');
            end
            for tx=1:4
                if  ~any(isnan(ML)) && ~any(isinf(ML)) && ML~=0
                set(tind(tx),'YLim',[0 1.1*ML])
                end
            end
            

            splitRates_1st_last=splitRatesC(:,[1 3],:);
%             fig2=figure('Position',[910         176        1000         710]);
%             tt= tiledlayout(1,2);
             label1=sprintf('%s<=%.3f',split_var_name,Ilim(2));
             labelend=sprintf('%s>%.3f',split_var_name,Ilim(end-1));
%             nexttile;
            [p_1st_vs_last_splitbin]=plot_summary(splitRates_1st_last,true(size(splitRates_1st_last)),{label1;labelend;label1;labelend},selc{1},{'Puff','Touch'},'mean',plotRatesAs,'dontplot');
%             ylim([-1.25 1.25])
%             title('VPM')
            paramOut{paramoutline,11}=p_1st_vs_last_splitbin;
            for s=1:2
                if p_1st_vs_last_splitbin(s)<0.05

                    x = b{1}(s).XEndPoints([1 end]);
                    y=repmat(nt(1).YLim(end)*.95+(s-1)*.025,1,2);
                    line(nt(1),x,y,'color','k')
                    text(nt(1),mean(x),y(1),'*','HorizontalAlignment','center','VerticalAlignment','bottom')
                end
            end
%             nexttile;
            [p_1st_vs_last_splitbin]=plot_summary(splitRates_1st_last,true(size(splitRates_1st_last)),{label1;labelend;label1;labelend},selc{2},{'Puff','Touch'},'mean',plotRatesAs,'dontplot');
%             title('POm')
%             ylim([-1.25 1.25])
%             title(tt,'1st last quantile comp',figadd,'Interpreter','none');
           % exportgraphics(fig2,[figdir 'Comp13_PT_' figadd '.pdf'],'BackgroundColor','none','ContentType','vector')
            paramOut{paramoutline,19}=p_1st_vs_last_splitbin;
            for s=1:2
                if p_1st_vs_last_splitbin(s)<0.05

                    x = b{2}(s).XEndPoints([1 end]);
                    y=repmat(nt(2).YLim(end)*.95+(s-1)*.025,1,2);
                    line(nt(2),x,y,'color','k')
                    text(nt(2),mean(x),y(1),'*','HorizontalAlignment','center','VerticalAlignment','bottom')
                end
            end


            exportgraphics(fig1,[figdir 'SplitPT_' figadd '.pdf'],'BackgroundColor','none','ContentType','vector')
        end
    end
end
