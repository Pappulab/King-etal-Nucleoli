
%%
load('NPM_1_2.mat');
figure('Units','centimeters','InnerPosition', [12 2 20 8]);
subplot(1,4,[1,2]);
%imagesc((SNARF_GFP_image(:,:,ii)+SNARF_GFP_image(:,:,ii+1)+SNARF_GFP_image(:,:,ii+2))/3); colormap('turbo');axis image; colorbar; caxis([100,13000]); 
imagesc(mean(GFP_only_image(:,:),3)); colormap('turbo');axis image; colorbar; %caxis([100,1000]); 
viscircles(condensates_center_GC,condensates_radius_GC,'Color','g','LineWidth',0.1,'LineStyle','--'); 
viscircles(condensates_center_outside,condensates_radius_outside,'Color','w','LineWidth',0.1,'LineStyle','--'); 
viscircles(condensates_center_DFC,condensates_radius_DFC,'Color','m','LineWidth',0.1,'LineStyle','--'); 
title('GFP only image'); axis off;  axis image;
xlim([condensates_center_GC(1)-condensates_radius_GC-40,condensates_center_GC(1)+condensates_radius_GC+40]); 
ylim([condensates_center_GC(2)-condensates_radius_GC-40,condensates_center_GC(2)+condensates_radius_GC+40]); 

%
subplot(1,4,[1,2]+2);
imagesc(sum(SNARF_GFP_image(:,:,8),3)); colormap('hot');axis image; colorbar; caxis([10000,50000]); 
viscircles(condensates_center_GFP,condensates_radius_GFP,'Color','g','LineWidth',0.1,'LineStyle','--'); 
viscircles(condensates_center_outside,condensates_radius_outside,'Color','w','LineWidth',0.1,'LineStyle','--'); 
viscircles(condensates_center_DFC,condensates_radius_DFC,'Color','m','LineWidth',0.1,'LineStyle','--'); 
title('SNARF image (636 nm)');  axis off; axis image;
xlim([condensates_center_GC(1)-condensates_radius_GC-40,condensates_center_GC(1)+condensates_radius_GC+40]); 
ylim([condensates_center_GC(2)-condensates_radius_GC-40,condensates_center_GC(2)+condensates_radius_GC+40]); 

%% pixel-wise ratio image
load("pythonMap.mat");

map_all = map_DFC | map_GC; map_all = double(map_all);

Fig1 = figure('Units','inches','InnerPosition',[1,1,5*0.5,4*0.5]); 

ax1 = axes('Position',[0.11,0.7,0.27,0.27]); 
imagesc(GFP_only_image(:,:)./max(GFP_only_image(:,:),[],'all')); colormap(ax1,cmap_GFP(1:200,:)); axis image; axis off; 
colorbar('east','Position',[0.36,0.72,0.017,0.2]); %caxis([100,1000]);  hold on;
%viscircles(condensates_center_GC,condensates_radius_GC,'Color','m','LineWidth',0.1,'LineStyle','--'); 
%viscircles(condensates_center_outside,condensates_radius_outside,'Color','m','LineWidth',0.1,'LineStyle','--'); 
%viscircles(condensates_center_DFC,condensates_radius_DFC,'Color','m','LineWidth',0.1,'LineStyle','--'); 

line([condensates_center_GC(1)-8+condensates_radius_GC,condensates_center_GC(1)-8+condensates_radius_GC+1000/pixel_size],...
    [condensates_center_GC(2)+condensates_radius_GC+4,condensates_center_GC(2)+condensates_radius_GC+4],'Color',[1,1,0.97],'LineWidth',1);
%title('GFP only image'); axis off;  axis image; 
xlim([condensates_center_GC(1)-condensates_radius_GC-10,condensates_center_GC(1)+condensates_radius_GC+10]); 
ylim([condensates_center_GC(2)-condensates_radius_GC-10,condensates_center_GC(2)+condensates_radius_GC+10]); 
%exportgraphics(Fig1,'BarChart.pdf','ContentType','vector')


ratio = mean(SNARF_GFP_image(:,:,2),3)./mean(SNARF_GFP_image(:,:,9),3);
%b = 1.4256962; a = 3664.4114885; y_min = 0.0407970; y_max = 0.7062725;
%pHCaculate = @(y) -(1/b).*log(y./a);
pHCaculate = @(y) (y/82840).^(-1/6.808);
pH_ca = pHCaculate(ratio);
ax2 = axes('Position',[0.11,0.40,0.27,0.27]); 
imagesc(pH_ca); caxis([6.6,7.5]); %caxis([0,0.25]);
% viscircles(condensates_center_GC,condensates_radius_GC,'Color','g','LineWidth',0.2,'LineStyle','--'); 
% viscircles(condensates_center_outside,condensates_radius_outside,'Color','w','LineWidth',0.2,'LineStyle','--'); 
% viscircles(condensates_center_DFC,condensates_radius_DFC,'Color','m','LineWidth',0.2,'LineStyle','--'); 
colormap_cur = turbo(300);axis image; colormap(ax2,colormap_cur(end:-1:1,:));  colorbar('east','Position',[0.36,0.42,0.017,0.2]);
xlim([condensates_center_GC(1)-condensates_radius_GC-10,condensates_center_GC(1)+condensates_radius_GC+10]); 
ylim([condensates_center_GC(2)-condensates_radius_GC-10,condensates_center_GC(2)+condensates_radius_GC+10]); 
axis off;

ax3 = axes('Position',[0.11,0.1,0.27,0.27]);
map_condensates = map_all;
map_condensates(map_condensates==0) = nan;
ratio_bur = imgaussfilt(mean(SNARF_GFP_image(:,:,2),3).*map_all,1.5)./imgaussfilt(mean(SNARF_GFP_image(:,:,9),3).*map_all,1.5);
pH_ca_bur = pHCaculate(ratio_bur);
imagesc(pH_ca_bur); caxis([6.6,7.5]);%caxis([0,0.25]);
% viscircles(condensates_center_GC,condensates_radius_GC,'Color','g','LineWidth',0.2,'LineStyle','--'); 
% viscircles(condensates_center_outside,condensates_radius_outside,'Color','w','LineWidth',0.2,'LineStyle','--'); 
% viscircles(condensates_center_DFC,condensates_radius_DFC,'Color','m','LineWidth',0.2,'LineStyle','--'); 
axis image; colormap(ax3,[[0,0,0];colormap_cur(end:-1:1,:)]); colorbar('east','Position',[0.36,0.12,0.017,0.2]);
xlim([condensates_center_GC(1)-condensates_radius_GC-10,condensates_center_GC(1)+condensates_radius_GC+10]); 
ylim([condensates_center_GC(2)-condensates_radius_GC-10,condensates_center_GC(2)+condensates_radius_GC+10]); 
axis off;

%%

fileName = {'NPM_1_1.mat','NPM_1_2.mat','NPM_1_3.mat','NPM_2_1.mat','NPM_2_2.mat'};


for ii = 1:length(fileName)

load(fileName{ii});

ratio = mean(SNARF_GFP_image(:,:,2),3)./mean(SNARF_GFP_image(:,:,9),3);

[X,Y] = meshgrid(1:size(GFP_only_image,1),1:size(GFP_only_image,1));
distance = sqrt((X-condensates_center_GC(1)).^2+(Y-condensates_center_GC(2)).^2);
Fig1 = figure('Units','inches','InnerPosition',[1,1,5,5]);
scatter(distance(distance<condensates_radius_GC*1.3)*pixel_size,ratio(distance<condensates_radius_GC*1.3),'filled','MarkerEdgeAlpha',0.1,'MarkerFaceAlpha',0.1); hold on;
%scatter(distance(distance<condensates_radius_GC*1.2)*pixel_size,GFP_only_image(distance<condensates_radius_GC*1.2)/max(GFP_only_image,[],'all'),'filled','MarkerEdgeAlpha',0.2,'MarkerFaceAlpha',0.2);
hold on;
plot([condensates_radius_GC,condensates_radius_GC]*pixel_size,[0.0,0.2],'k--','LineWidth',1);
plot([condensates_radius_DFC,condensates_radius_DFC]*pixel_size,[0.0,0.2],'k--','LineWidth',1);
plot([condensates_peak_GC,condensates_peak_GC]*pixel_size,[0.0,0.2],'k--','LineWidth',1);

distance_set = linspace(0,condensates_radius_GC*1.3,30);
for jj = 1:length(distance_set)-1
    median_ratio(jj) = median(ratio(distance>=distance_set(jj) & distance<distance_set(jj+1)));
    mean_GFP(jj) = mean(GFP_only_image(distance>=distance_set(jj) & distance<distance_set(jj+1)))/max(GFP_only_image(distance<(condensates_radius_GC*1.2)),[],'all');
    
end
plot((distance_set(1:end-1)+distance_set(2:end))/2*pixel_size,median_ratio,'r*-','LineWidth',1);

yyaxis left
ylim([0.1,0.185]);
yyaxis right
plot((distance_set(1:end-1)+distance_set(2:end))/2*pixel_size,mean_GFP,'k','LineWidth',1);

plot([0,condensates_radius_GC*1.4*pixel_size],[(max(mean_GFP)+mean_GFP(end))/2,(max(mean_GFP)+mean_GFP(end))/2],'k');
plot([0,condensates_radius_GC*1.4*pixel_size],[(max(mean_GFP)+min(mean_GFP(1:20)))/2,(max(mean_GFP)+min(mean_GFP(1:20)))/2],'k');
pots_ratio_save{ii} = ratio(distance<condensates_radius_GC*1.3);
pots_x_save{ii}  = distance(distance<condensates_radius_GC*1.3)/condensates_radius_GC;
pots_x_save_abs{ii}  = distance(distance<condensates_radius_GC*1.3);
mean_GFP_save(ii,:) = mean_GFP;
median_ratio_save(ii,:) = median_ratio;
distance_set_save = linspace(0,1.3,30);

int_SNARF_outside(ii,:) = squeeze(nanmean(map_outside.*SNARF_GFP_image,[1,2])).';
int_SNARF_outside(ii,:) = int_SNARF_outside(ii,:)./int_SNARF_outside(ii,9);
end


% b = 1.4256962; a = 3664.4114885; y_min = 0.0407970; y_max = 0.7062725;
% pHCaculate = @(y) -(1/b)*log(y/a);
pHCaculate = @(y) (y/82840).^(-1/6.808);

ratio_outside = int_SNARF_outside(:,2)./int_SNARF_outside(:,9);
pH_outside = median(pHCaculate(ratio_outside));
norm_factor = 7.21-pH_outside;


%% plot individual one
load('NPM_1_2.mat'); idx_Nu = 2;


Fig1 = figure('Units','inches','InnerPosition',[1,1,1.9,1.9]);  hold on; box on;
pH = pHCaculate(median_ratio_save)+norm_factor;
scatter(pots_x_save_abs{idx_Nu}*pixel_size/1000,pHCaculate(pots_ratio_save{2})+norm_factor,1,'filled','MarkerEdgeAlpha',0.2,'MarkerFaceAlpha',0.2,'MarkerEdgeColor',[224, 177, 180]/255,'MarkerFaceColor',[224, 177, 180]/255); hold on;
[radius_set_output, pH_set_output,image_scatter] = histogram_scatter(pots_x_save_abs{idx_Nu}*pixel_size/1000,pHCaculate(pots_ratio_save{2})+norm_factor);
load('mapPurple.mat'); load('mapHorizon.mat');
%imagesc(radius_set_output, pH_set_output,image_scatter);  colormap(mapPurple([1,10:120],:));  colorbar('south','Position',[0.25,0.8,0.1,0.03]); 
%plot((distance_set_save(1:end-1)+distance_set_save(2:end))/2,pH,'Color',[224, 177, 180]/255,'LineWidth',2); hold on;
plot((distance_set_save(1:end-1)+distance_set_save(2:end))/2*condensates_radius_GC*pixel_size/1000,pH(idx_Nu,:),'r','LineWidth',2); hold on;
plot([1,1]*condensates_radius_GC*pixel_size/1000,[6.4,7.4],'k--','LineWidth',1);
plot([1,1]*condensates_radius_DFC*pixel_size/1000,[6.4,7.4],'k--','LineWidth',1);
plot([1,1]*condensates_peak_GC*pixel_size/1000,[6.4,7.4],'k--','LineWidth',1);

xlabel('Radial distance (\mum)');

yyaxis left;  
ylim([6.4,7.4]); ylabel('pH'); xlim([0,3.6]);
yyaxis right;   
%plot((distance_set_save(1:end-1)+distance_set_save(2:end))/2,mean_GFP_save,'Color',[166, 191, 178]/255,'LineWidth',2,'Marker','none','LineStyle','-'); hold on;

plot((distance_set_save(1:end-1)+distance_set_save(2:end))/2*condensates_radius_GC*pixel_size/1000,mean_GFP_save(idx_Nu,:)./max(mean_GFP_save(idx_Nu,:)),'Color',[12, 112, 57]/255,'LineWidth',2,'Marker','none');
GFP_temp = mean(mean_GFP_save,1);
%plot([0,1.3],[(max(GFP_temp)+GFP_temp(end))/2,(max(GFP_temp)+GFP_temp(end))/2],'k');
%plot([0,1.3],[(max(GFP_temp)+min(GFP_temp(1:20)))/2,(max(GFP_temp)+min(GFP_temp(1:20)))/2],'k');
ylabel('Normalized intensity');  set(gca,'YColor',[53, 92, 64]/255);
%% plot individual one
load('NPM_1_2.mat'); idx_Nu = 2;


Fig1 = figure('Units','inches','InnerPosition',[1,1,1.9,1.9]);  hold on; box on;
pH = pHCaculate(median_ratio_save)+norm_factor;
%scatter(pots_x_save_abs{idx_Nu}*pixel_size/1000,pHCcaculate(pots_ratio_save{2})+norm_factor,2,'filled','MarkerEdgeAlpha',0.3,'MarkerFaceAlpha',0.3,'MarkerEdgeColor',[224, 177, 180]/255,'MarkerFaceColor',[224, 177, 180]/255); hold on;
[radius_set_output, pH_set_output,image_scatter] = histogram_scatter(pots_x_save_abs{idx_Nu}*pixel_size/1000,pHCaculate(pots_ratio_save{2})+norm_factor);
load('mapPurple.mat'); load('mapHorizon.mat');
imagesc(radius_set_output, pH_set_output,image_scatter);  colormap(mapPurple([1,10:120],:));  colorbar('south','Position',[0.25,0.8,0.1,0.03]); 
%plot((distance_set_save(1:end-1)+distance_set_save(2:end))/2,pH,'Color',[224, 177, 180]/255,'LineWidth',2); hold on;
plot((distance_set_save(1:end-1)+distance_set_save(2:end))/2*condensates_radius_GC*pixel_size/1000,pH(idx_Nu,:),'r','LineWidth',2); hold on;
plot([1,1]*condensates_radius_GC*pixel_size/1000,[6.2,7.4],'k--','LineWidth',1);
plot([1,1]*condensates_radius_DFC*pixel_size/1000,[6.2,7.4],'k--','LineWidth',1);
plot([1,1]*condensates_peak_GC*pixel_size/1000,[6.2,7.4],'k--','LineWidth',1);

xlabel('Radial distance (\mum)');

yyaxis left;  
ylim([6.4,7.4]); ylabel('pH'); xlim([0,3.6]);
yyaxis right;   
%plot((distance_set_save(1:end-1)+distance_set_save(2:end))/2,mean_GFP_save,'Color',[166, 191, 178]/255,'LineWidth',2,'Marker','none','LineStyle','-'); hold on;

plot((distance_set_save(1:end-1)+distance_set_save(2:end))/2*condensates_radius_GC*pixel_size/1000,mean_GFP_save(idx_Nu,:)./max(mean_GFP_save(idx_Nu,:)),'Color',[12, 112, 57]/255,'LineWidth',2,'Marker','none');
GFP_temp = mean(mean_GFP_save,1);
%plot([0,1.3],[(max(GFP_temp)+GFP_temp(end))/2,(max(GFP_temp)+GFP_temp(end))/2],'k');
%plot([0,1.3],[(max(GFP_temp)+min(GFP_temp(1:20)))/2,(max(GFP_temp)+min(GFP_temp(1:20)))/2],'k');
ylabel('Normalized intensity');  set(gca,'YColor',[53, 92, 64]/255); 

%% plot individual one
Fig1 = figure('Units','inches','InnerPosition',[1,1,1.9,1.9]); 
pH = pHCaculate(median_ratio_save)+norm_factor;
%plot((distance_set_save(1:end-1)+distance_set_save(2:end))/2,pH,'Color',[224, 177, 180]/255,'LineWidth',2); hold on;
plot((distance_set_save(1:end-1)+distance_set_save(2:end))/2,median(pH,1),'r','LineWidth',2); hold on;
P1 = prctile(pH,95,1);
P2 = prctile(pH(:,end:-1:1),5,1);
x = [(distance_set_save(1:end-1)+distance_set_save(2:end))/2]; x=[x,x(end:-1:1)];
y_pH = [P1,P2];
pgon_pH = polyshape(x,y_pH);
plot(pgon_pH,'FaceColor',[237, 61, 50]/255,'FaceAlpha',0.5,'EdgeColor','none');
plot([1,1],[6.4,7.4],'k--','LineWidth',1);
plot([0.6117,0.6117],[6.4,7.4],'k--','LineWidth',1);
plot([0.82931,0.82931],[6.4,7.4],'k--','LineWidth',1);



yyaxis left
ylim([6.4,7.4]); ylabel('pH');  xlim([0,1.3]); xlabel('Radial distance (normalized)');
yyaxis right
%plot((distance_set_save(1:end-1)+distance_set_save(2:end))/2,mean_GFP_save,'Color',[166, 191, 178]/255,'LineWidth',2,'Marker','none','LineStyle','-'); hold on;

P1 = prctile(mean_GFP_save,95,1);
P2 = prctile(mean_GFP_save(:,end:-1:1),5,1);
x = [(distance_set_save(1:end-1)+distance_set_save(2:end))/2]; x=[x,x(end:-1:1)];
y_GFP = [P1,P2];
pgon_GFP= polyshape(x,y_GFP./max(y_GFP));
plot((distance_set_save(1:end-1)+distance_set_save(2:end))/2,median(mean_GFP_save,1)./max(y_GFP),'Color',[12, 112, 57]/255,'LineWidth',2,'Marker','none');
plot(pgon_GFP,'FaceColor',[12, 112, 57]/255,'FaceAlpha',0.5,'EdgeColor','none');
ylim([0,1]);

%plot([0,1.3],[(max(GFP_temp)+GFP_temp(end))/2,(max(GFP_temp)+GFP_temp(end))/2],'k');
%plot([0,1.3],[(max(GFP_temp)+min(GFP_temp(1:20)))/2,(max(GFP_temp)+min(GFP_temp(1:20)))/2],'k');
ylabel('Normalized intensity');  set(gca,'YColor',[53, 92, 64]/255);  
%% functions
function im_ch = read_tiff(filename, imageN)

imageR = Tiff(filename);
for ii = 1:imageN
    imageR.setDirectory(ii);
    im_ch(:,:,ii) = imageR.read();
end

end



function [radius_set_output, pH_set_output,image_scatter] = histogram_scatter(radius,pH)

pH_set = 6:0.02:8;
radius_set = 0:0.1:5;

radius_set_output = (radius_set(1:end-1)+radius_set(2:end))/2;
pH_set_output = (pH_set(1:end-1)+pH_set(2:end))/2;

for ii = 2:length(pH_set)
   for jj = 2:length(radius_set)
        idx = radius<radius_set(jj) & radius>radius_set(jj-1) & pH<pH_set(ii) & pH>pH_set(ii-1);
        image_scatter(ii,jj) = sum(idx);

   end
end


end