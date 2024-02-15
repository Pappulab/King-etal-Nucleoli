
%%
load('NPM_1_2.mat');
figure('Units','centimeters','InnerPosition', [12 2 20 8]);
subplot(1,4,[1,2]); 
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

%% example images for figure
load("pythonMap.mat");

Fig1 = figure('Units','inches','InnerPosition',[1,1,5*0.45,8*0.45]); 

ax1 = axes('Position',[0.01,0.7,0.27,0.27]);
imagesc(mean(GFP_only_image(:,:),3)); colormap(ax1,cmap_GFP(1:200,:)); axis image; axis off; %colorbar; %caxis([100,1000]);  hold on;


line([condensates_center_GC(1)-8+condensates_radius_GC,condensates_center_GC(1)-8+condensates_radius_GC+1000/pixel_size],...
    [condensates_center_GC(2)+condensates_radius_GC+4,condensates_center_GC(2)+condensates_radius_GC+4],'Color',[1,1,0.97],'LineWidth',1);
xlim([condensates_center_GC(1)-condensates_radius_GC-10,condensates_center_GC(1)+condensates_radius_GC+10]); 
ylim([condensates_center_GC(2)-condensates_radius_GC-10,condensates_center_GC(2)+condensates_radius_GC+10]); 


ax2 = axes('Position',[0.31,0.7,0.27,0.27]);
imagesc(SNARF_GFP_image(:,:,2)); colormap(ax2,cmap_red(:,:)); axis image;  axis off; %colorbar; %caxis([100,1000]);  hold on;
xlim([condensates_center_GC(1)-condensates_radius_GC-10,condensates_center_GC(1)+condensates_radius_GC+10]); 
ylim([condensates_center_GC(2)-condensates_radius_GC-10,condensates_center_GC(2)+condensates_radius_GC+10]); 

ax3 = axes('Position',[0.61,0.7,0.27,0.27]);
imagesc(SNARF_GFP_image(:,:,9)); colormap(ax3,cmap_red(:,:)); axis image; axis off;  %colorbar; %caxis([100,1000]);  hold on;

xlim([condensates_center_GC(1)-condensates_radius_GC-10,condensates_center_GC(1)+condensates_radius_GC+10]); 
ylim([condensates_center_GC(2)-condensates_radius_GC-10,condensates_center_GC(2)+condensates_radius_GC+10]); 

%% functions
function im_ch = read_tiff(filename, imageN)

imageR = Tiff(filename);
for ii = 1:imageN
    imageR.setDirectory(ii);
    im_ch(:,:,ii) = imageR.read();
end

end
