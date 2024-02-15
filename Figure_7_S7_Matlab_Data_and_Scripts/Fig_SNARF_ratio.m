


fileName = {'NPM_1_1.mat','NPM_1_2.mat','NPM_1_3.mat','NPM_2_1.mat','NPM_2_2.mat'};


%% caculate the intensity at each eavelengh

for ii = 1:length(fileName)

load(fileName{ii});
map_DFC = double(map_DFC); map_DFC(map_DFC==0) = nan;
map_GC = double(map_GC); map_GC(map_GC==0) = nan;
map_outside = double(map_outside); map_outside(map_outside==0) = nan;
map_FC = double(map_FC); map_FC(map_FC==0) = nan;
map_GC_outer = double(map_GC_outer); map_GC_outer(map_GC_outer==0) = nan;
map_GC_inner = double(map_GC_inner); map_GC_inner(map_GC_inner==0) = nan;

int_SNARF_GC(ii,:) = squeeze(nanmean(map_GC.*SNARF_GFP_image,[1,2])).';
int_SNARF_DFC(ii,:) = squeeze(nanmean(map_DFC.*SNARF_GFP_image,[1,2])).';
int_SNARF_outside(ii,:) = squeeze(nanmean(map_outside.*SNARF_GFP_image,[1,2])).';
int_SNARF_GC_outer(ii,:) = squeeze(nanmean(map_GC_outer.*SNARF_GFP_image,[1,2])).';
int_SNARF_GC_inner(ii,:) = squeeze(nanmean(map_GC_inner.*SNARF_GFP_image,[1,2])).';

int_SNARF_DFC(ii,:) = int_SNARF_DFC(ii,:)./int_SNARF_outside(ii,9);
int_SNARF_GC(ii,:) = int_SNARF_GC(ii,:)./int_SNARF_outside(ii,9);
int_SNARF_GC_outer(ii,:) = int_SNARF_GC_outer(ii,:)./int_SNARF_outside(ii,9);
int_SNARF_GC_inner(ii,:) = int_SNARF_GC_inner(ii,:)./int_SNARF_outside(ii,9);
int_SNARF_outside(ii,:) = int_SNARF_outside(ii,:)./int_SNARF_outside(ii,9);


wavelength = linspace(570,693,14);
Fig1 = figure('Units','inches','InnerPosition',[1,1,10,3]);
subplot(1,2,1); 
plot(wavelength,int_SNARF_GC(ii,:)); hold on;
plot(wavelength,int_SNARF_DFC(ii,:));
plot(wavelength,int_SNARF_outside(ii,:));
plot(wavelength,int_SNARF_GC_outer(ii,:),'k-');
plot(wavelength,int_SNARF_GC_inner(ii,:),'k-');
%plot([579,579],[0,30000],'k--');
%plot([646,646],[0,30000],'k--');
xlim([wavelength(1),wavelength(end)]); title('intensity');
legend('GC','DFC','outside','EdgeColor','none','color','none');
xlabel('wavelength(nm)'); ylabel('intensity');

end

%%
wavelength = linspace(570,693,14);
Fig1 = figure('Units','inches','InnerPosition',[1,1,1.9,1.9]); 

plot(wavelength,median(int_SNARF_DFC,1),'Color',"#01889F",'LineWidth',2); hold on;
plot(wavelength,median(int_SNARF_GC,1),'Color',"#7E2F8E",'LineWidth',2); 
plot(wavelength,median(int_SNARF_outside,1),'Color',"#808080",'LineWidth',2);

x = [wavelength,wavelength(end:-1:1)];
P1 = prctile(int_SNARF_GC,95,1); P2 = prctile(int_SNARF_GC(:,end:-1:1),5,1);
y_GC = [P1,P2];
pgon_GC = polyshape(x,y_GC);
P1 = prctile(int_SNARF_DFC,95,1); P2 = prctile(int_SNARF_DFC(:,end:-1:1),5,1);
y_DFC = [P1,P2];
pgon_DFC = polyshape(x,y_DFC);
P1 = prctile(int_SNARF_outside,95,1); P2 = prctile(int_SNARF_outside(:,end:-1:1),5,1);
y_outside = [P1,P2];
pgon_outside = polyshape(x,y_outside);
plot(pgon_DFC,'FaceColor',"#01889F",'FaceAlpha',0.5,'EdgeColor','none');
plot(pgon_GC,'FaceColor',"#7E2F8E",'FaceAlpha',0.5,'EdgeColor','none');
plot(pgon_outside,'FaceColor',"#808080",'FaceAlpha',0.5,'EdgeColor','none');

plot([579,579],[0,1.1],'k--','LineWidth',2);
plot([646,646],[0,1.1],'k--','LineWidth',2);
xlim([wavelength(1),wavelength(end)]);
legend('FC/DFC','GC','NP','EdgeColor','none','color','none'); ylim([0,1.1]);
xlabel('wavelength(nm)'); ylabel('Normalized intensity');

%%
Fig1 = figure('Units','inches','InnerPosition',[1,1,1.8,1.9]); 
X = categorical({'NP','GC(in)','GC(out)','FC/DFC'});
X = reordercats(X,{'FC/DFC','GC(in)','GC(out)','NP',});
ratio_GC = int_SNARF_GC(:,2)./int_SNARF_GC(:,9);
ratio_GC_inner = int_SNARF_GC_inner(:,2)./int_SNARF_GC_inner(:,9);
ratio_GC_outer = int_SNARF_GC_outer(:,2)./int_SNARF_GC_outer(:,9);
ratio_DFC = int_SNARF_DFC(:,2)./int_SNARF_DFC(:,9);
ratio_outside = int_SNARF_outside(:,2)./int_SNARF_outside(:,9);


pHCaculate = @(y) (y/82840).^(-1/6.808);


pH_GC = median(pHCaculate(ratio_GC));
pH_GC_inner = median(pHCaculate(ratio_GC_inner));
pH_GC_outer = median(pHCaculate(ratio_GC_outer));
pH_outside = median(pHCaculate(ratio_outside));
pH_DFC  = median(pHCaculate(ratio_DFC));

std_GC = std(pHCaculate(ratio_GC));
std_GC_inner = std(pHCaculate(ratio_GC_inner));
std_GC_outer = std(pHCaculate(ratio_GC_outer));
std_outside = std(pHCaculate(ratio_outside));
std_DFC = std(pHCaculate(ratio_DFC));
norm_factor = 7.21-pH_outside;

%b = bar(X,[ratio_outside,ratio_GC, ratio_DFC]); hold on;
bar(X(1),pH_outside+norm_factor,'FaceColor',"#808080");   hold on;
%bar(X(2),pH_GC+norm_factor,'FaceColor',"#7E2F8E");
bar(X(2),pH_GC_inner+norm_factor,'FaceColor',"#7E2F8E");
bar(X(3),pH_GC_outer+norm_factor,'FaceColor',"#7E2F8E");
bar(X(4),pH_DFC+norm_factor,'FaceColor',"#01889F"); ylim([6.2,7.4]);

errorbar(X,[pH_outside,pH_GC_inner,pH_GC_outer,pH_DFC]+norm_factor,[std_outside,std_GC_inner,std_GC_outer,std_DFC],"LineStyle","none",'CapSize',10,'LineWidth',2,'Color','k');
ylabel('pH');  
