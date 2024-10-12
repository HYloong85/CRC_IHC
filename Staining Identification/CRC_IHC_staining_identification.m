clc;
clear;
openslide_load_library();
load('net.mat');
svsList = dir('./data/TMA*.svs');
maskList = dir('./mask_new/mask_out_*.mat');
marginList = dir('./output_margin/margin*.xlsx');
n = length(svsList);
for i = 1:n 
    %% If there are multiple tissue specimens annotated with bounding boxes in an image, record the coordinates of each bounding box.
    % xmlDoc = xml_read(['./dataTMA',svsname,'.xml']);
    % [m, ~] = size(xmlDoc.Annotation.Regions.Region);
    % pts = zeros(m, 4);
    % for j = 1:m
    %     pts(j, :) = [xmlDoc.Annotation.Regions.Region(j).Vertices.Vertex(1).ATTRIBUTE.X,...
    %         xmlDoc.Annotation.Regions.Region(j).Vertices.Vertex(1).ATTRIBUTE.Y,...
    %         xmlDoc.Annotation.Regions.Region(j).Vertices.Vertex(3).ATTRIBUTE.X,...
    %         xmlDoc.Annotation.Regions.Region(j).Vertices.Vertex(3).ATTRIBUTE.Y];
    % end
    % %Obtain annotated images based on the coordinates of two vertices；
    % x = round(pts(j, 1));
    % y = round(pts(j, 2));
    % w = round(pts(j, 3) - pts(j, 1));
    % h = round(pts(j, 4) - pts(j, 2));
    
    img = openslide_open(['./data/', svsList(i).name]);
    [w, h] = openslide_get_level_dimensions(img,0);
    image = openslide_read_region(img, 0, 0, w,h, 0);
    %Convert four channels to three channels；
    image = image(:, :, 2:4);
    imR = image(:, :, 1);
    imG = image(:, :, 2);
    imB = image(:, :, 3);
    rgb = [imR(:)'; imG(:)'; imB(:)'];
    
    
    % Use model prediction；
    predict = net(rgb);
    [~, predictClass] = max(predict);
   


    %% Calculate the proportion of stained tissue
     margin = xlsread(['./output_margin/', marginList(i).name]);
     mask = [];
    for z = 1:n
        if (mean(maskList(z).name(10:16) == svsList(i).name(1:7)) == 1) 
            mask = load(['./mask_new/', maskList(z).name]);
            mask = mask.mask_new;
            break;
        end
    end
    [height, width] = size(mask);
    image_region = zeros(height*72, width*72);
    for p = 1:height
        for q = 1:width
            image_region(72*(p-1)+1 : 72*(p-1)+72, 72*(q-1)+1 : 72*(q-1)+72) = mask(p, q);
        end
    end
    right = ones(height*72, 152).*mask(height, width);
    down = ones(152, width*72 + 152).*mask(height, width);
    image_region = [image_region, right];
    image_region = [image_region; down];
    [height, width] = size(image_region);
    margin_up = ones(margin(1, 1), width);
    margin_down = ones(margin(1, 2), width);
    margin_left = ones(height + margin(1, 1) + margin(1, 2), margin(1, 1));
    margin_right = ones(height + margin(1, 1) + margin(1, 2), margin(1, 3));
    image_region = [margin_up; image_region; margin_down];
    image_region = [margin_left, image_region, margin_right];

    image_region = image_region(:)';
    
    predictClass_GLA = predictClass(image_region(:)==2);  %GLA；
    predictClass_NOR_STR = predictClass(image_region(:)==3);  %NORM_STR；
    predictClass_TUM = predictClass(image_region(:)==4);  %TUM；
    predictClass_TUM_STR = predictClass(image_region(:)==5);  %TUM_STR；
    
    n0 = sum(predictClass~=1)/length(predictClass);  %Proportion of tissues；
    n1 = sum(predictClass_GLA==3)/sum(predictClass_GLA~=1);  %Proportion of staining in the gland；
    n2 = sum(predictClass_GLA~=1)/sum(predictClass~=1);  %Proportion of glands in all tissues；
    n3 = sum(predictClass_NOR_STR==3)/sum(predictClass_NOR_STR~=1);  %Proportion of staining in the normal stroma；
    n4 = sum(predictClass_NOR_STR~=1)/sum(predictClass~=1);  %Proportion of normal stroma in all tissues；
    n5 = sum(predictClass_TUM==3)/sum(predictClass_TUM~=1);  %Proportion of staining in the tumor；
    n6 = sum(predictClass_TUM~=1)/sum(predictClass~=1);  %Proportion of tumors in all tissues；
    n7 = sum(predictClass_TUM_STR==3)/sum(predictClass_TUM_STR~=1);  %Proportion of staining in the tumor stroma；
    n8 = sum(predictClass_TUM_STR~=1)/sum(predictClass~=1);  %The proportion of tumor stroma in all tissues；
    tissuePercent(1, :) = [n0, n1, n2, n3, n4, n5, n6, n7, n8];
    save(['./tissuePercent_', svsList(i).name(1:end-4), '.mat'], 'tissuePercent');


    %% Visualization of stained tissue classification results
    predictClass = cat(3, predictClass, predictClass, predictClass);
    predictClass(1,predictClass(1, :, 1) == 1,1) = 255;
    predictClass(1,predictClass(1, :, 2) == 1,2) = 255;
    predictClass(1,predictClass(1, :, 3) == 1,3) = 255;
    predictClass(1,predictClass(1, :, 1) == 2,1) = 139;
    predictClass(1,predictClass(1, :, 2) == 2,2) = 137;
    predictClass(1,predictClass(1, :, 3) == 2,3) = 137;
    predictClass(1,predictClass(1, :, 1) == 3,1) = 128;
    predictClass(1,predictClass(1, :, 2) == 3,2) = 64;
    predictClass(1,predictClass(1, :, 3) == 3,3) = 0;
    
    img = uint8(zeros(h,w,3));
    for z = 1:3
        for p = 1:w
            img(:, p, z) = predictClass(1, (p-1)*h+1:p*h, z);
        end
    end
    imshow(img);
    imwrite(img, ['./',svsList(i).name(1:end-4),'_stain.tiff']);

end

