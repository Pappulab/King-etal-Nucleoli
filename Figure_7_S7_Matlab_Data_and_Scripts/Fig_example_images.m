% 
% %% input
% % fileFoder and fileName
% % define circular regions of the nucleoli by giving 'condensates_center' &
% % 'condensates_radius'
% %% ****input folder name********
% fileFolder = 'F:\OneDrive - Washington University in St. Louis\github\SNARF\230310\tiff file\';
% scanning_SNARF_name = '230223_npm1_pH6_lambda_n1.tif';
% %*** if didn't capture GFP only image, leave the name to be 'scanning_GFP_name = Nan;'*** 
% scanning_GFP_name = '230223_npm1_pH6_gfp_n1.tif';
% 
% condensates_Nub = 2;
% %% read image
% 
% SNARF_GFP_image = double(read_tiff([fileFolder, scanning_SNARF_name],14));
% if isnan(scanning_GFP_name)==0
% GFP_only_image = double(read_tiff([fileFolder,scanning_GFP_name],1));
% else
%     GFP_only_image = SNARF_GFP_image*0;
% end
% % SNARF_GFP_image_temp = zeros(932,932,14);
% % %SNARF_GFP_image_temp((1:810)+31,(1:810)+31,:) = SNARF_GFP_image;
% % %SNARF_GFP_image_temp((1:810)+36,(1:810)+47,:) = SNARF_GFP_image;
% % SNARF_GFP_image_temp((1:810)+105,(1:810)+87,:) = SNARF_GFP_image;
% % SNARF_GFP_image = SNARF_GFP_image_temp;
% SNARF_GFP_image = imresize(SNARF_GFP_image,size(GFP_only_image,1)/size(SNARF_GFP_image,1));
% 
% pixel_size = 1000/11.7635;% nm
% diffraction_list = 300; %nm;
% 
% %% define condensates region
% %* manually input 'condensates_center' & 'condensates_radius'
% % *this part will generate images at each wavelength, and save the image to
% % folder 'analyse_image_save'
% if condensates_Nub==1
% 
% condensates_center_GC = [232,250];  % in unit of pixel
% condensates_radius_GC =32; % in unit of pixel
% condensates_center_outside = [200,200]; % in unit of pixel
% condensates_radius_outside = 20;
% condensates_center_DFC = [232,250];  % in unit of pixel
% condensates_radius_DFC = 21; % in unit of pixel
% 
% condensates_center_GFP = [232,250]; % in unit of pixel
% condensates_radius_GFP = 35; % in unit of pixel
% 
% map_FC = GFP_only_image(:,:)>1000000000;
% 
% elseif condensates_Nub==2
% 
% condensates_center_GC = [378,275];  % in unit of pixel
% condensates_radius_GC = 33; % in unit of pixel
% condensates_center_outside = [300,200]; % in unit of pixel
% condensates_radius_outside = 20;
% condensates_center_DFC = [378,275];  % in unit of pixel
% condensates_radius_DFC = 23; % in unit of pixel
% 
% condensates_center_GFP = [378,275]; % in unit of pixel
% condensates_radius_GFP = 34; % in unit of pixel
% 
% map_FC = GFP_only_image(:,:)>1000000000;
% 
% elseif condensates_Nub==3
% 
% condensates_center_GC = [662,827];  % in unit of pixel
% condensates_radius_GC = 34; % in unit of pixel
% condensates_center_outside = [600,800]; % in unit of pixel
% condensates_radius_outside = 20;
% condensates_center_DFC = [662,827];   % in unit of pixel
% condensates_radius_DFC = 23; % in unit of pixel
% 
% condensates_center_GFP = [662,827];  % in unit of pixel
% condensates_radius_GFP = 34; % in unit of pixel
% 
% map_FC = GFP_only_image(:,:)>1000000000;
% 
% end
% 
% x = 1:size(GFP_only_image,2);
% y = 1:size(GFP_only_image,1);
% [X,Y] = meshgrid(x,y);
% distance = sqrt((X-condensates_center_GC(1)).^2+(Y-condensates_center_GC(2)).^2);
% map_all = distance*0;
% map_all(distance<=condensates_radius_GC)=1; 
% 
% distance = sqrt((X-condensates_center_outside(1)).^2+(Y-condensates_center_outside(2)).^2);
% map_outside_w_all = distance*0;
% map_outside_w_all(distance<=condensates_radius_outside)=1; 
% 
% distance = sqrt((X-condensates_center_DFC(1)).^2+(Y-condensates_center_DFC(2)).^2);
% map_DFC_w_FC = distance*0;
% map_DFC_w_FC(distance<=condensates_radius_DFC)=1; 
% 
% map_DFC = map_FC==0 & map_DFC_w_FC==1;
% map_GC = map_DFC_w_FC==0 & map_all==1;
% map_outside = map_all==0 & map_outside_w_all==1;
% 
% 
% %create folder for save images
% fileSave = [fileFolder,'analyse_image_save_condensate2'];
% if ~exist(fileSave, 'dir')
%    mkdir(fileSave)
% end

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

%% example images for figure
load("pythonMap.mat");

Fig1 = figure('Units','inches','InnerPosition',[1,1,5*0.45,8*0.45]); 

ax1 = axes('Position',[0.01,0.7,0.27,0.27]);
imagesc(mean(GFP_only_image(:,:),3)); colormap(ax1,cmap_GFP(1:200,:)); axis image; axis off; %colorbar; %caxis([100,1000]);  hold on;
%viscircles(condensates_center_GC,condensates_radius_GC,'Color','m','LineWidth',0.1,'LineStyle','--'); 
%viscircles(condensates_center_outside,condensates_radius_outside,'Color','m','LineWidth',0.1,'LineStyle','--'); 
%viscircles(condensates_center_DFC,condensates_radius_DFC,'Color','m','LineWidth',0.1,'LineStyle','--'); 

line([condensates_center_GC(1)-8+condensates_radius_GC,condensates_center_GC(1)-8+condensates_radius_GC+1000/pixel_size],...
    [condensates_center_GC(2)+condensates_radius_GC+4,condensates_center_GC(2)+condensates_radius_GC+4],'Color',[1,1,0.97],'LineWidth',1);
%title('GFP only image'); axis off;  axis image; 
xlim([condensates_center_GC(1)-condensates_radius_GC-10,condensates_center_GC(1)+condensates_radius_GC+10]); 
ylim([condensates_center_GC(2)-condensates_radius_GC-10,condensates_center_GC(2)+condensates_radius_GC+10]); 
%exportgraphics(Fig1,'BarChart.pdf','ContentType','vector')


ax2 = axes('Position',[0.31,0.7,0.27,0.27]);
imagesc(SNARF_GFP_image(:,:,2)); colormap(ax2,cmap_red(:,:)); axis image;  axis off; %colorbar; %caxis([100,1000]);  hold on;
%viscircles(condensates_center_GC,condensates_radius_GC,'Color','m','LineWidth',0.1,'LineStyle','--'); 
%viscircles(condensates_center_outside,condensates_radius_outside,'Color','m','LineWidth',0.1,'LineStyle','--'); 
%viscircles(condensates_center_DFC,condensates_radius_DFC,'Color','m','LineWidth',0.1,'LineStyle','--'); 

%line([condensates_center_GC(1)-8+condensates_radius_GC,condensates_center_GC(1)-8+condensates_radius_GC+1000/pixel_size],...
%    [condensates_center_GC(2)+condensates_radius_GC+8,condensates_center_GC(2)+condensates_radius_GC+8],'Color',[1,1,0.97],'LineWidth',2);
%title('GFP only image'); axis off;  axis image; 
xlim([condensates_center_GC(1)-condensates_radius_GC-10,condensates_center_GC(1)+condensates_radius_GC+10]); 
ylim([condensates_center_GC(2)-condensates_radius_GC-10,condensates_center_GC(2)+condensates_radius_GC+10]); 

ax3 = axes('Position',[0.61,0.7,0.27,0.27]);
imagesc(SNARF_GFP_image(:,:,9)); colormap(ax3,cmap_red(:,:)); axis image; axis off;  %colorbar; %caxis([100,1000]);  hold on;
%viscircles(condensates_center_GC,condensates_radius_GC,'Color','m','LineWidth',0.1,'LineStyle','--'); 
%viscircles(condensates_center_outside,condensates_radius_outside,'Color','m','LineWidth',0.1,'LineStyle','--'); 
%viscircles(condensates_center_DFC,condensates_radius_DFC,'Color','m','LineWidth',0.1,'LineStyle','--'); 

%line([condensates_center_GC(1)-8+condensates_radius_GC,condensates_center_GC(1)-8+condensates_radius_GC+1000/pixel_size],...
%    [condensates_center_GC(2)+condensates_radius_GC+8,condensates_center_GC(2)+condensates_radius_GC+8],'Color',[1,1,0.97],'LineWidth',2);
%title('GFP only image'); axis off;  axis image; 
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
