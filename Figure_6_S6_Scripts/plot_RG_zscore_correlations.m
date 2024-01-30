clear all;

constructs={'POLR1F','NUCKS','NUCKSdelK'};
mytitles={'POLR1F','NUCKS','NUCKSdelK'};
posneg=[8.37 7.02 5.36];
pospos=[8.93 2.99 2.2];
kblock=[3.62 1.87 0.48];

reps=[5 5 5];


colorhex={'#90298d','#2a3b8f','#28a8e0','#008181'};
for i=1:length(colorhex)
    mycolor(i,:)=sscanf(colorhex{i}(2:end),'%2x%2x%2x',[1 3])/255;
end
mycolor2=mycolor(end:-1:1,:);

[hs, seqs]=fastaread('../Figure_6_S6_Data/nucleolar_idrs.fasta');
%bins=0:1:60;
bins=0:0.15:6;

count=0;
%% Plot all reps
for c=1:length(constructs)
    count=count+1;
    pos=find(strcmp(hs,constructs{c})==1);
    myseq=seqs{pos};

    nrgr=[];
    for r=1:reps(c)
        da=importdata(['../Figure_6_S6_Data/' constructs{c} '/340/' num2str(r) '/ana/Rg.dat']); 
        mrgr(r)=mean(da/sqrt(length(myseq)));
        clear da; 
    end

    % Plot histograms
    rgmean(c)=mean(mrgr);
    rgstd(c)=std(mrgr);
    rgmedian(c)=median(mrgr);
    rgste(c)=std(mrgr)/sqrt(reps(c));
end

figure;
subplot(1,3,1)
for c=1:3 errorbar(posneg(c),rgmean(c),rgste(c),rgste(c),'o','color','k','markeredgecolor','k','markerfacecolor',mycolor(c,:),'markersize',20); hold on; 
    xlabel('pos-neg')
    ylabel('<Rg/sqrt(N)>')
end
%legend(mytitles)
ylim([2 4])
xlim([5 9])
[fitvals,s]=polyfit(posneg,rgmean,1)
x=5:0.1:9;
y=fitvals(1)*x+fitvals(2);
plot(x,y,'-k','linewidth',4); hold on; 
[yfit, dy]=polyconf(fitvals,x,s,'predopt','curve');
%plot(x,yfit - dy,'-k');
%plot(x,yfit + dy,'-k');
%patch([x, fliplr(x)], [yfit - dy fliplr(yfit + dy)], [211 211 211]/255, 'EdgeColor','none', 'FaceAlpha',0.5)
[f gof]=fit(posneg',rgmean','poly1')
text(6,2.5,['R^2 = ' num2str(gof.rsquare)]);


subplot(1,3,2)
for c=1:3 errorbar(pospos(c),rgmean(c),rgste(c),rgste(c),'o','color','k','markeredgecolor','k','markerfacecolor',mycolor(c,:),'markersize',20); hold on; 
    xlabel('pos-pos')
    ylabel('<Rg/sqrt(N)>')
end
%legend(mytitles)
ylim([2 4])
xlim([2 9])
[fitvals,s]=polyfit(pospos,rgmean,1)
x=2:0.1:9;
y=fitvals(1)*x+fitvals(2);
plot(x,y,'-k','linewidth',4); hold on; 
[yfit, dy]=polyconf(fitvals,x,s,'predopt','curve');
%plot(x,yfit - dy,'-k');
%plot(x,yfit + dy,'-k');
%patch([x, fliplr(x)], [yfit - dy fliplr(yfit + dy)], [211 211 211]/255, 'EdgeColor','none', 'FaceAlpha',0.5)
[f gof]=fit(pospos',rgmean','poly1')
text(3,2.5,['R^2 = ' num2str(gof.rsquare)]);

subplot(1,3,3)
for c=1:3 errorbar(kblock(c),rgmean(c),rgste(c),rgste(c),'o','color','k','markeredgecolor','k','markerfacecolor',mycolor(c,:),'markersize',20); hold on; 
    xlabel('K Block')
    ylabel('<Rg/sqrt(N)>')
end
%legend(mytitles)
ylim([2 4])
xlim([0 4])
[fitvals,s]=polyfit(kblock,rgmean,1)
x=0:0.1:4;
y=fitvals(1)*x+fitvals(2);
plot(x,y,'-k','linewidth',4); hold on; 
[yfit, dy]=polyconf(fitvals,x,s,'predopt','curve');
%plot(x,yfit - dy,'-k');
%plot(x,yfit + dy,'-k');
%patch([x, fliplr(x)], [yfit - dy fliplr(yfit + dy)], [211 211 211]/255, 'EdgeColor','none', 'FaceAlpha',0.5)
[f gof]=fit(kblock',rgmean','poly1')
text(1,2.5,['R^2 = ' num2str(gof.rsquare)]);

