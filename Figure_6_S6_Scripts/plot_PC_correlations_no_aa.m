clear all;

%% Import data
constructs={'POLR1F','NUCKS','NUCKSdelK'};
mytitles={'POLR1F','NUCKS','NUCKSdelK'};
posneg=[8.37 7.02 5.36];
pospos=[8.93 2.99 2.2];
kblock=[3.62 1.87 0.48];

isilmean=[28.2 13.7 1.96];
isilste=[0.4 0.5 0.08];
isilnorm=(isilmean-min(isilmean))./(max(isilmean)-min(isilmean));




%figure;
%x=[1,2];
%bar(x,isilmean(2:end)); hold on;
%er = errorbar(x,isilmean(2:end),isilste(2:end),isilste(2:end));    
%er.Color = [0 0 0];                            
%er.LineStyle = 'none'; 
%return

ivtda=importdata('../Figure_6_S6_Data/Compiled_PC_InVitro.csv');
for i=1:3
    ivtmean(i)=nanmean(ivtda.data(:,i));
    ivtstd(i)=nanstd(ivtda.data(:,i));
    ivtmedian(i)=nanmedian(ivtda.data(:,i));
    tmp=nnz(~isnan(ivtda.data(:,i)))
    ivtste(i)=nanstd(ivtda.data(:,i))/sqrt(tmp);
end

ivvda=importdata('../Figure_6_S6_Data/Compiled_DFC.csv');
for i=1:3
    ivvmean(i)=nanmean(ivvda.data(:,i));
    ivvstd(i)=nanstd(ivvda.data(:,i));
    ivvmedian(i)=nanmedian(ivvda.data(:,i));
    tmp=nnz(~isnan(ivvda.data(:,i)))
    ivvste(i)=nanstd(ivvda.data(:,i))/sqrt(tmp);
end

colorhex={'#90298d','#2a3b8f','#28a8e0','#008181'};
for i=1:length(colorhex)
    mycolor(i,:)=sscanf(colorhex{i}(2:end),'%2x%2x%2x',[1 3])/255;
end
mycolor2=mycolor(end:-1:1,:);


%% Plot data

% Just in silico data
f=figure;
f.Position=[100 100 280 350];
for c=1:3 errorbar(c,isilmean(c),isilste(c),'o','color','k','markeredgecolor','k','markerfacecolor',mycolor(c,:),'markersize',24); hold on; 
    ylabel('in silico Partition Coefficient')
end
%legend(mytitles)
ylim([0 30])
xlim([0 4])
set(gca,'xtick',[1:1:3])
set(gca,'xticklabel',{'PolR1F IDR','WT NUCKS','delK-rich NUCKS'}); 
%return

figure;
subplot(1,4,1)
for c=1:3 errorbar(kblock(c),isilmean(c),isilste(c),'o','color','k','markeredgecolor','k','markerfacecolor',mycolor(c,:),'markersize',20); hold on; 
    xlabel('K Block z-score')
    ylabel('in silico PC')
end
%legend(mytitles)
ylim([0 30])
xlim([0 4])
[fitvals,s]=polyfit(kblock,isilmean,1)
x=0:0.1:30;
y=fitvals(1)*x+fitvals(2);
plot(x,y,'-k','linewidth',4); hold on; 
[yfit, dy]=polyconf(fitvals,x,s,'predopt','curve');
%plot(x,yfit - dy,'-k');
%plot(x,yfit + dy,'-k');
%patch([x, fliplr(x)], [yfit - dy fliplr(yfit + dy)], [211 211 211]/255, 'EdgeColor','none', 'FaceAlpha',0.5)
[f gof]=fit(kblock',isilmean','poly1')

subplot(1,4,2)
for c=1:3 errorbar(isilmean(c),ivtmean(c),ivtste(c),ivtste(c),isilste(c),isilste(c),'o','color','k','markeredgecolor','k','markerfacecolor',mycolor(c,:),'markersize',20); hold on; 
    xlabel('in silico PC')
    ylabel('in vitro PC')
end
%legend(mytitles)
ylim([0 11])
xlim([0 30])
[fitvals,s]=polyfit(isilmean,ivtmean,1)
x=0:0.1:30;
y=fitvals(1)*x+fitvals(2);
plot(x,y,'-k','linewidth',4); hold on; 
[yfit, dy]=polyconf(fitvals,x,s,'predopt','curve');
%plot(x,yfit - dy,'-k');
%plot(x,yfit + dy,'-k');
%patch([x, fliplr(x)], [yfit - dy fliplr(yfit + dy)], [211 211 211]/255, 'EdgeColor','none', 'FaceAlpha',0.5)
[f gof]=fit(isilmean',ivtmean','poly1')


subplot(1,4,3)
for c=1:3 errorbar(ivtmean(c),ivvmean(c),ivvste(c),ivvste(c),ivtste(c),ivtste(c),'o','color','k','markeredgecolor','k','markerfacecolor',mycolor(c,:),'markersize',20); hold on; 
    xlabel('in vitro PC')
    ylabel('in vivo PC')
end
%legend(mytitles)
ylim([0 5])
xlim([0 11])
[fitvals,s]=polyfit(ivtmean,ivvmean,1)
x=0:0.1:12;
y=fitvals(1)*x+fitvals(2);
plot(x,y,'-k','linewidth',4); hold on; 
[yfit, dy]=polyconf(fitvals,x,s,'predopt','curve');
%plot(x,yfit - dy,'-k');
%plot(x,yfit + dy,'-k');
%patch([x, fliplr(x)], [yfit - dy fliplr(yfit + dy)], [211 211 211]/255, 'EdgeColor','none', 'FaceAlpha',0.5)
[f gof]=fit(ivtmean',ivvmean','poly1')



type=[kblock; isilmean; ivtmean; ivvmean];
for t1=1:length(type)
    for t2=t1:length(type)
        [f gof]=fit(type(t1,:)',type(t2,:)','poly1')
        r2mat(t1,t2)=gof.rsquare;
    end
end

subplot(1,4,4)

xvalues = {'K Block','in silico','in vitro','in vivo'};
yvalues = {'K Block','in silico','in vitro','in vivo'};
h=heatmap(xvalues,yvalues,r2mat);
h.CellLabelFormat = '%.2f';
h.Title = 'R^2 Correlation';
caxis([0.95 1])

% Create colormap: 
%map = brewermap(8,'GnBu'); 
map = brewermap(8,'Greens'); 
map(1,:) = [1 1 1]; % optionally force first color to white 
colormap(map)

%c = gray;
%c = flipud(c);
%colormap(c);



