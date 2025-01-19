function [RasterHP,RasterHT,sp_tr_Pall,sp_tr_Tall,paramOut]=split_kinematics(split_time_before_trig,split_time_after_trig,split_safety_margin,Wdname,split_binsize,smoothwidth,meth,binsR,...
    selc,splits,splitname,DiscreteData,ppms,split_var_nameC,type_nameC,split_allR,figdir,TrialWp,TrialWt,Amp,Phase,Angle,Setp)


selN={'VPM';'POm'};
selS={'Puff';'Touch'};
% bins=-split_time_before_trig:split_binsize:split_time_after_trig;
bins=binsR(1):split_binsize:binsR(end);
pbins=bins(1:end-1);
pbins=pbins+mean(diff(pbins)/2);

paramOut={};
paramoutline=0;
% splits={[0 25 50 75 100];[0 100/3 200/3 100];[0 50 100]};
% splitname={'Quartile';'Tertile';'median'};

isi=10;
%Puff
parameter='Puff';
exclude={'Pole';'Light';'Exclude';'Grooming'};
triglength_limit=[28 32];%ms
tempDD=struct2cell(DiscreteData);
spikes=squeeze(tempDD(1,1,:));clear temp*
[RasterHP,P_trig]=get_raster_pop_no_self_safety(parameter,spikes, DiscreteData,triglength_limit,ppms,split_time_before_trig,split_time_after_trig,exclude,split_safety_margin,split_binsize,isi,[]);

parameter='Touch';
exclude={'Puff';'Light';'Exclude';'Grooming'};
triglength_limit=[-inf inf];%ms
[RasterHT,T_trig]=get_raster_pop_no_self_safety(parameter,spikes, DiscreteData,triglength_limit,ppms,split_time_before_trig,split_time_after_trig,exclude,split_safety_margin,split_binsize,isi,[]);



%%


%
close all
sp_tr_T=cell(size(T_trig,1),numel(split_var_nameC),numel(type_nameC));
sp_tr_P=cell(size(T_trig,1),numel(split_var_nameC),numel(type_nameC));
for SN=1:numel(split_var_nameC)
    split_var_name=split_var_nameC{SN};
    switch split_var_name
        case 'Interval'
            %             sp_tr_T=cell(size(T_trig));
            for n=1:numel(RasterHT)
                if ~isempty(T_trig{n})
                    sp_tr_T{n,SN,1}=cat(1,T_trig{n}(1),diff(T_trig{n}))/20000;
                end
            end
            sp_tr_P(:,SN,1)=cellfun(@(x) cat(1,x(1),diff(x))/20000,P_trig,'UniformOutput',0);

        case {'Angle','Setp','Phase','Curvature','Acceleration','Velocity','Amp'}
            switch split_var_name
                case {'Angle','Setp','Amp'}
                    eval(sprintf('a=%s;',split_var_name));
                    switch split_var_name
                        case {'Angle','Setp'}
                            a=cellfun(@(x) x-median(x,'omitnan'),a,'UniformOutput',0);
                    end
                case {'Acceleration','Velocity'}
                    a=load(Wdname, split_var_name);
                    a=a.(split_var_name);
                case 'Phase'
                    a=Phase;
                    for x=1:numel(a)
                        
                        a{x}(Amp{x}<3)=nan;
                    end
                case 'Curvature'
                    a=load(Wdname, split_var_name,'Amp','Curvature');
                    amp=a.Amp;
                    a=a.(split_var_name);

                    a=cellfun(@(x,y) x-median(x(y<3),'omitnan'),a,amp,'UniformOutput',0);
                    a=cellfun(@(x) x./std(x,0,'omitnan'),a,'UniformOutput',0);
            end
            pbinsx=pbins(pbins<25 & pbins>-25);
            [~,temp_P]=get_trig_cont_pop(P_trig,pbinsx,a,20,DiscreteData);
            [~,temp_T]=get_trig_cont_pop(T_trig,pbinsx,a,20,DiscreteData);
            xT=cellfun(@isempty,temp_T);
            xP=cellfun(@isempty,temp_P);
            temp_P(xP)={nan(size(pbinsx))};
            temp_T(xT)={nan(size(pbinsx))};
            for tx=1:numel(type_nameC)
                %{'rawbefore','rawafter','absbefore','absafter','rel','absrel'};
                switch type_nameC{tx}
                    case 'rawbefore'
                        sp_tr_P(:,SN,tx)=cellfun(@(x) mean(x(:,pbinsx<0),2,'omitnan'),temp_P,'UniformOutput',0);
                        sp_tr_T(:,SN,tx)=cellfun(@(x) mean(x(:,pbinsx<0),2,'omitnan'),temp_T,'UniformOutput',0);
                    case 'rawafter'
                        sp_tr_P(:,SN,tx)=cellfun(@(x) mean(x(:,pbinsx>0),2,'omitnan'),temp_P,'UniformOutput',0);
                        sp_tr_T(:,SN,tx)=cellfun(@(x) mean(x(:,pbinsx>0),2,'omitnan'),temp_T,'UniformOutput',0);
                    case 'absbefore'
                        sp_tr_P(:,SN,tx)=cellfun(@(x) abs(mean(x(:,pbinsx<0),2,'omitnan')),temp_P,'UniformOutput',0);
                        sp_tr_T(:,SN,tx)=cellfun(@(x) abs(mean(x(:,pbinsx<0),2,'omitnan')),temp_T,'UniformOutput',0);
                    case 'absafter'
                        sp_tr_P(:,SN,tx)=cellfun(@(x) abs(mean(x(:,pbinsx>0),2,'omitnan')),temp_P,'UniformOutput',0);
                        sp_tr_T(:,SN,tx)=cellfun(@(x) abs(mean(x(:,pbinsx>0),2,'omitnan')),temp_T,'UniformOutput',0);
                    case 'rel'
                        sp_tr_P(:,SN,tx)=cellfun(@(x) mean(x(:,pbinsx>0),2,'omitnan')-mean(x(:,pbinsx<0),2,'omitnan'),temp_P,'UniformOutput',0);
                        sp_tr_T(:,SN,tx)=cellfun(@(x) mean(x(:,pbinsx>0),2,'omitnan')-mean(x(:,pbinsx<0),2,'omitnan'),temp_T,'UniformOutput',0);
                    case 'absrel'
                        sp_tr_P(:,SN,tx)=cellfun(@(x) abs(mean(x(:,pbinsx>0),2,'omitnan')-mean(x(:,pbinsx<0),2,'omitnan')),temp_P,'UniformOutput',0);
                        sp_tr_T(:,SN,tx)=cellfun(@(x) abs(mean(x(:,pbinsx>0),2,'omitnan')-mean(x(:,pbinsx<0),2,'omitnan')),temp_T,'UniformOutput',0);
                end





            end
            sp_tr_P(xP,SN,:)={[]};
            sp_tr_T(xT,SN,:)={[]};
    end


end

sp_tr_T=cellfun(@(x) reshape(x,[],1),sp_tr_T,'UniformOutput',0);
sp_tr_P=cellfun(@(x) reshape(x,[],1),sp_tr_P,'UniformOutput',0);
sp_tr_Tall=cell(1,size(sp_tr_T,2),size(sp_tr_T,3));
sp_tr_Pall=sp_tr_Tall;
for SN=1:size(sp_tr_T,2)
    for ty=1:size(sp_tr_T,3)
        sp_tr_Tall{1,SN,ty}=cat(1,sp_tr_T{:,SN,ty});
        sp_tr_Pall{1,SN,ty}=cat(1,sp_tr_P{:,SN,ty});
    end
end

%%

[TrialWpN]=normalize_trig_traces(TrialWp,2,1:20);
[TrialWtN]=normalize_trig_traces(TrialWt,2,1:20);
%% split by all together
PTCol=[49 94 91; 214 143 0]/255;

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
        for TY=1:size(sp_tr_P,3)
            figadd=strcat(split_var_nameC{SN},'_',type_nameC{TY},'_',js);
            temp1=cat(1,cat(1,sp_tr_P{:,SN,TY}),cat(1,sp_tr_T{:,SN,TY}));
            temp1(isnan(temp1))=[];
            if isempty(temp1);continue;end
            %%
            paramoutline=paramoutline+1;
            paramOut{paramoutline,1}=split_var_name;
            paramOut{paramoutline,2}=type_nameC{TY};
           
            Ilim=prctile(temp1,splits);
            paramOut{paramoutline,3}=Ilim(2:end-1);
            ix1=discretize(temp1,Ilim);
            mean_value_disc=nan(max(ix1),1);
            for j=1:max(ix1)
                mean_value_disc(j)=mean(temp1(ix1==j));
            end
            labels=cell(size(Ilim,2)-1,1);
            for l=1:numel(Ilim)-2
                labels{l}=sprintf('v<= %.3f',Ilim(l+1));
            end
            labels{end}=sprintf('v> %.3f',Ilim(end-1));

            ST=cell(numel(T_trig),numel(Ilim)-1,2);
            ST2=ST;
            TTn=ST;
            TT=ST;
            STH=ST;
            STR=cell(numel(T_trig),numel(binsR)-1,2);
            STHc=ST;
            %puff trial split
            Climx=nan(numel(T_trig),numel(Ilim));
            for n=1:numel(T_trig)
                if split_all
                    Clim=Ilim;
                else
                    Clim=prctile(cat(1,sp_tr_P{n,SN,TY},sp_tr_T{n,SN,TY}),splits);
                end
                Climx(n,:)=Clim;
                if ~isempty(RasterHP{n}) && ~isempty(sp_tr_P{n,SN,TY}) && ~any(n==[80 81])
                    ix=discretize(sp_tr_P{n,SN,TY},Ilim)';
                    sb=ix(RasterHP{n}(:,2));
                    for x=1:numel(Ilim)-1
                        TTn{n,x,1}=mean(TrialWpN{n}(ix==x,:),1,'omitnan');
                        TT{n,x,1}=mean(TrialWp{n}(ix==x,:),1,'omitnan');
                        ST{n,x,1}=RasterHP{n}(sb==x,1);
                        ST2{n,x,1}=RasterHP{n}(sb==x,:);
                        STH{n,x,1}=histcounts(ST{n,x,1},bins)./sum(ix==x);
                        STR{n,x,1}=histcounts(ST{n,x,1},binsR)./sum(ix==x)./(diff(binsR)/1000);
                        if sum(ix==x) <=1
                            STH{n,x,1}(:)=nan;
                            STR{n,x,1}(:)=nan;
                        end
                        %      STHc{n,x,1}=conv(STH{n,x,1},smoothker2,'same');
                    end
                else
                    STH(n,:,1)={nan(1,numel(pbins))};
                    STR(n,:,1)={nan(1,numel(binsR)-1)};
                    %   STHc(n,:,2)={nan(1,numel(pbins))};
                end
                %touch trial split
                if ~isempty(RasterHT{n}) && ~isempty(sp_tr_T{n,SN,TY}) && ~any(n==[80 81])
                    ix=discretize(sp_tr_T{n,SN,TY},Clim);
                    sb=ix(RasterHT{n}(:,2));
                    for x=1:numel(Ilim)-1
                         TTn{n,x,2}=mean(TrialWtN{n}(ix==x,:),1,'omitnan');
                         TT{n,x,2}=mean(TrialWt{n}(ix==x,:),1,'omitnan');
                        ST{n,x,2}=RasterHT{n}(sb==x,1);
                        STH{n,x,2}=histcounts(ST{n,x,2},bins)./sum(ix==x);
                        STR{n,x,2}=histcounts(ST{n,x,2},binsR)./sum(ix==x)./(diff(binsR)/1000);
                        if sum(ix==x) <=1
                            STH{n,x,2}(:)=nan;
                            STR{n,x,2}(:)=nan;
                        end
                        %  STHc{n,x,2}=conv(STH{n,x,2},smoothker2,'same');
                    end
                else
                    STH(n,:,2)={nan(1,numel(pbins))};
                    STR(n,:,2)={nan(1,numel(binsR)-1)};
                    %   STHc(n,:,2)={nan(1,numel(pbins))};
                end
            end
            stn=cellfun(@numel,ST);
            normSN=ones(numel(T_trig),1);
            STR(80,:,2)={nan(1,numel(binsR)-1)};
            STR(81,:,2)={nan(1,numel(binsR)-1)};
            STR(80,:,1)={nan(1,numel(binsR)-1)};
            STR(81,:,1)={nan(1,numel(binsR)-1)};
            STRs=STR;

            %
            C=repmat(linspace(.7,0,numel(Ilim)-1)',1,3);
%%
            fig1= figure('Position',[50          97        1055        1184]);
            tx=0;
            ML=0;

            tt=tiledlayout(5,6);

            title(tt,sprintf('PSTHs, split by %s of %s',splitname{S},split_var_name),sprintf('normalization: %s, split: %s',type_nameC{TY},js),'Interpreter','none')
            nexttile(1,[1 2]);
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

            H1=histcounts(cat(1,sp_tr_P{:,SN,TY}),vbins,'Normalization','probability');
            H2=histcounts(cat(1,sp_tr_T{:,SN,TY}),vbins,'Normalization','probability');
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
                    %         line(Ilim([2 2]),[0 1.25*max([H1 H2])],'Color','k');
                    %         line(Ilim([end-1 end-1]),1.25*[0 max([H1 H2])],'Color','k');
                    line(repmat(Ilim(2:end-1),2,1),[0 1.25*max([H1 H2])],'Color','k');
            end
            box off
            legend(selS,'Location','eastoutside');legend boxoff
            title(split_var_name)

            STR=cellfun(@(x) (x(2)-x(1))./(x(2)+x(1)),STRs);
            %STR=cellfun(@(x) x(2)-x(1),STRs);
            %STR=cellfun(@(x) log2(x(2)./x(1)),STRs);
            %  STR=cellfun(@(x) 2.*(x(2)-x(1))./x(1),STRs);
            for nx=1:2
                nt(nx)=nexttile(1+nx*2,[1 2]);
                n=selc{nx};
                tempR=STR(n,:,:);
                tempR(isinf(tempR))=nan;

                mean_tempR=squeeze(mean(tempR,1,'omitnan'));
                sem_tempR=squeeze(std(tempR,0,1,'omitnan')./sqrt(sum(~isnan(tempR))));
               

                tempp=nan(size(tempR,3),size(tempR,2));
                for pt=1:size(tempR,3)
                    for tx=1:size(tempR,2)
                        if any(~isnan(tempR(:,tx,pt)))
                            tempp(pt,tx)=signrank(tempR(:,tx,pt),[],'tail','right');
                        end
                    end
                end
                paramOut{paramoutline,4+(nx-1)*8}=mean_tempR(:,1)';
                paramOut{paramoutline,5+(nx-1)*8}=sem_tempR(:,1)';
                paramOut{paramoutline,6+(nx-1)*8}=tempp(1,:);
                paramOut{paramoutline,7+(nx-1)*8}=mean_tempR(:,2)';
                paramOut{paramoutline,8+(nx-1)*8}=sem_tempR(:,2)';
                paramOut{paramoutline,9+(nx-1)*8}=tempp(2,:);
                tempp=tempp(:);
                tempind=cell(size(tempp));
                for np=1:numel(tempp)
                    if tempp(np)<.001
                        tempind{np}='***';
                    elseif tempp(np)<.01
                        tempind{np}='**';
                    elseif    tempp(np)<.05
                        tempind{np}='*';
                    else
                        %tempind{np}='n.s.';
                        tempind{np}='';
                    end
                end

                temppc=nan(1,size(tempR,2));

                for tx=1:size(tempR,2)
                    if any(all(~isnan(tempR(:,tx,[1 2])),3))
                        temppc(1,tx)=signrank(tempR(:,tx,1),tempR(:,tx,2));
                    end
                end
                paramOut{paramoutline,10+(nx-1)*8}=temppc(1,:);
                tempindc=cell(size(temppc));
                for np=1:numel(temppc)
                    if temppc(np)<.001
                        tempindc{np}='***';
                    elseif temppc(np)<.01
                        tempindc{np}='**';
                    elseif    temppc(np)<.05
                        tempindc{np}='*';
                    else
%                         tempindc{np}='n.s.';
                        tempindc{np}='';
                    end
                end

                
                b=bar(1:size(mean_tempR,1),mean_tempR);hold on;
                x = nan(size(mean_tempR,2), size(mean_tempR,1));
                for i = 1:size(mean_tempR,2)
                    x(i,:) = b(i).XEndPoints;
                    b(i).FaceColor=PTCol(i,:);
                end

                errorbar(x',mean_tempR,sem_tempR,'k','linestyle','none');
                text(sort(x(:)),reshape(mean_tempR'+sem_tempR',[],1),tempind,'HorizontalAlignment','center','VerticalAlignment','bottom')

                line(x,repmat(1.15*max(mean_tempR+sem_tempR,[],'all'),size(x,1),size(x,2)),'color','k')
                text(mean(x,1),repmat(1.15*max(mean_tempR+sem_tempR,[],'all'),1,size(x,2)),tempindc,'HorizontalAlignment','center','VerticalAlignment','bottom')
                yl=ylim;
                yl(2)=1.35*max(mean_tempR+sem_tempR,[],'all');
                ylim(yl);
                hold off
                ylabel('modulation')
%                 yt=get(gca,'YTick');
                xlabel('tertile')
                legend(selS,'Location','eastoutside');legend boxoff
                box off
                cR=nan(2,2);
                for pt=1:2
                    tempR2=tempR(:,:,pt);
                    tempR2=tempR2(~(sum(isnan(tempR2),2)>=size(tempR2,2)-1),:);

                    temp_mvd=repmat(mean_value_disc',size(tempR2,1),1);
                    %           temp_mvd=repmat(1:3,size(tempR2,1),1);
                    tempR2=tempR2(:);temp_mvd=temp_mvd(:);
                    re=~isnan(tempR2);
                    tempR2=tempR2(re);
                    temp_mvd=temp_mvd(re);
                    [c,p]=corrcoef(temp_mvd,tempR2);
                    cR(pt,1)=c(2,1);
                    cR(pt,2)=p(2,1);
                end
                %title(selN{nx},sprintf('c(P)=%.2f, p(P)=%.6f, c(T)=%.2f, p(T)=%.6f',cR(1,1),cR(1,2),cR(2,1),cR(2,2)))
                title(selN{nx});
            end



            tx=0;

            starttile=6;
            for s=1:2

                for nx=1:2
                    n=selc{nx};
                    tx=tx+1;
                    tind(tx)=nexttile(starttile+tx*3-2,[1 3]);
                    temp= cell2mat(STH(n,:,s))./normSN(n);
                    temp=reshape(temp,size(temp,1),numel(bins)-1,size(STH,2));
                    %         SEM=permute(std(temp,0,1,'omitnan')./sqrt(sum(~any(isnan(temp),[2 3]))),[2 3 1]);
                    %         SEM=smoothdata(SEM,1,meth,smoothwidth,'SamplePoints',milliseconds(pbins));
                    M=permute(mean(temp,1,'omitnan'),[2 3 1]);
                    M=smoothdata(M,1,meth,smoothwidth,'SamplePoints',milliseconds(pbins));
                    for b=1:size(M,2)
                        %             SEMp=cat(1,M(:,b)+SEM(:,b),flipud(M(:,b)-SEM(:,b)));
                        %             fill(fillpbins, SEMp, C(b,:),'FaceAlpha',.5,'EdgeColor','none');
                        hold on;
                        %plot(pbins',conv(M(:,b),smoothker2,'same'),'Color',C(b,:));box off;
                        plot(pbins',M(:,b),'Color',C(b,:));box off;

                    end

                    ML=max(ML,max(M,[],'all'));
                    legend(labels);legend('boxoff');
                    xlabel('time around deflection (ms)')
                    ylabel({'Population mean';'spike prob. per bin'})
                    title(sprintf('%s (N=%u) %s',selN{nx},sum(any(~all(isnan(temp),2),3)),selS{s}));
                end
            end

            %
            for tx=1:4
                set(tind(tx),'YLim',[0 1.1*ML])
            end

            tx=0;

            starttile=18;

            temp=TT(cat(1,selc{:}),:,:);
            tempn=TTn(cat(1,selc{:}),:,:);
            for x=1:size(temp,2)
                for tt=1:2
                    temp2{1,x,tt}=mean(cat(1,temp{:,x,tt}),1,'omitnan');
                    temp2{2,x,tt}=mean(cat(1,tempn{:,x,tt}),1,'omitnan');
                end
            end
            temp3{1,1}=cat(1,temp2{1,:,1});%raw puff
            temp3{1,2}=cat(1,temp2{1,:,2});%raw touch
            temp3{2,1}=cat(1,temp2{2,:,1});%norm puff
            temp3{2,2}=cat(1,temp2{2,:,2});%norm touch
            trig_trace{SN}=temp3;

            for s=1:2
                for t=1:2
                    tx=tx+1;
                    tind(tx)=nexttile(starttile+tx*3-2,[1 3]);
                    for b=1:size(temp3{s,t},1)
                        plot(pbins'+25,temp3{s,t}(b,:),'Color',C(b,:));box off;hold on;
                    end
                    xlabel('time around deflection (ms)')
                    if s==1
                        ylabel('whisker angle (raw)')
                    else
                        ylabel('whisker angle (norm)')
                    end

                    xlim([-50 50])
                end



            end
      

%%



            exportgraphics(fig1,[figdir 'SplitPT_' figadd '.pdf'],'BackgroundColor','none','ContentType','vector')

            %%

            STR=STR(:,[1 3],:);
            %%
            fig2=figure('Position',[910         176        1000         710]);
            tt= tiledlayout(1,2);
            label1=sprintf('%s<=%.3f',split_var_name,Ilim(2));
            labelend=sprintf('%s>%.3f',split_var_name,Ilim(3));
            nexttile;
            [pI,pC,mS]=plot_ratecomp3(STR,true(size(STR)),{label1;labelend;label1;labelend},selc{1},{'Puff','Touch'},'mean','dRate post-pre [Hz]');
            ylim([-1.25 1.25])
            title('VPM')
             paramOut{paramoutline,11}=pI';
            nexttile;
            [pI,pC,mS]=plot_ratecomp3(STR,true(size(STR)),{label1;labelend;label1;labelend},selc{2},{'Puff','Touch'},'mean','dRate post-pre [Hz]');
            title('POm')
             ylim([-1.25 1.25])
            title(tt,'Tert 1/3 comp',figadd,'Interpreter','none');
            exportgraphics(fig2,[figdir 'Comp13_PT_' figadd '.pdf'],'BackgroundColor','none','ContentType','vector')
             paramOut{paramoutline,19}=pI';

        end
        %%
    end


end
