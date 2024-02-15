scriptPath = mfilename('fullpath');
[fileFolder, ~, ~] = fileparts(scriptPath);
scanning_SNARF_name = 'NCL_2uM_1ngPre_Lamda_Frame1.tif';

%% read image

SNARF_GFP_image = double(read_tiff([fileFolder, scanning_SNARF_name],14));
ratio = SNARF_GFP_image(:,:,2)./SNARF_GFP_image(:,:,9);

b = 1.4256962; a = 3664.4114885; y_min = 0.0407970; y_max = 0.7062725;
pHCaculate = @(y) -(1/b).*log(y./a);
pH_ca = pHCaculate(ratio);

figure(); imagesc(pH_ca); caxis([6.5,7.4]); axis image;




%% functions
function im_ch = read_tiff(filename, imageN)

imageR = Tiff(filename);
for ii = 1:imageN
    imageR.setDirectory(ii);
    im_ch(:,:,ii) = imageR.read();
end

end