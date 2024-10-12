clc;
clear;
datapath = ('./output_masks/');
outputpath = ('./mask_new/');
maskList = dir([datapath, 'mask_out*.mat']);
N = length(maskList);
addpath(genpath('.\function\'));

%% Divide the stroma into normal stroma and tumor stroma
border = 20;
for i = 1:N
    % Four types of tissues: 1:Disturb  2:GLA  3:STR  4:TUM
    mask = load([datapath, maskList(i).name]);
    mask = mask.mask_out;
    [x, y, ~] = size(mask);
    mask_new = zeros(x, y);
    for m = 1:x
        for n = 1:y
            mask_new(m, n) = max(find(mask(m, n, :) == max(mask(m, n, :)))); 
        end
    end  
    % If it is paracancerous tissue, remove interference from surrounding glands;
    prop = zeros(1,4);
    for m = 1:4
        prop(m) = length(mask_new(mask_new(:) == m)) / length(mask_new(:));
    end
    if mod(str2double(maskList(i).name(25:end-4)), 8) == 7 || mod(str2double(maskList(i).name(25:end-4)), 8) == 0
        if prop(3) > 0.1
            mask_new = str2nor_tum_imdilate(mask_new);
        end
        prop = zeros(1,4);
        for m = 1:4
            prop(m) = length(mask_new(mask_new(:) == m)) / length(mask_new(:));
        end
    end  
    
    if mod(str2double(maskList(i).name(25:end-4)), 8) == 1 || mod(str2double(maskList(i).name(25:end-4)), 8) == 2
        mask_new(mask_new(:) == 3) = 3;  %If the tissue specimen is normal tissue, then all stroma is normal stroma；
    else
        if mod(str2double(maskList(i).name(25:end-4)), 8) == 5 || mod(str2double(maskList(i).name(25:end-4)), 8) == 6
            mask_new(mask_new(:) == 3) = 5;  %If the tissue specimen is the tumor center tissue, then all stroma is normal stroma；
        else
            if  prop(2) > prop(4)*3
                mask_new(mask_new(:) == 3) = 3;  %If the proportion of glands in the tissue specimen is much larger than the proportion of tumors, then it is all normal stroma；
            else
                if prop(2)*3 < prop(4)
                    mask_new(mask_new(:) == 3) = 5;  %If the proportion of tumors in the tissue specimen is much larger than the proportion of glands, then it is all tumor stroma；
                else
                    % If it aren't for the above situation , check the adjacent tissue categories of the interstitium
                    for m = 1:x
                        for n = 1:y
                            side = str2nor_tum_side(m, n, x, y, border);%side: left; right; up; down;
                            
                            result = str2nor_tum(mask_new(m, n),mask_new, side, prop);
                            mask_new(m, n) = result;
                        end
                    end
                end
            end
        end
    end
    save ([outputpath, maskList(i).name], 'mask_new');
    
end

%% Save the visualization classification results of the new five types of tissues；
%% 1:Disturb  2:GLA  3:NORM_STR  4:TUM  5:TUM_STR
prop_new = zeros(1,5);
for m = 1:5
    prop_new(m) = length(mask_new(mask_new(:) == m)) / length(mask_new(:));
end  

img = uint8(zeros(x, y, 3));
for m = 1:x
    for n = 1:y
        if mask_new(m, n) == 1
            img(m, n, :) = 255; % Disturb: White
        end
        if mask_new(m, n) == 2
            img(m, n, 1) = 255; % GLA: Red；
        end
        if mask_new(m, n) == 3
            img(m, n, 2) = 255; % NORM_STR: Green；
        end
        if mask_new(m, n) == 4
            img(m, n, 3) = 255; % TUM: Blue；
        end
        if mask_new(m, n) == 5
            img(m, n, 1) = 255;
            img(m, n, 3) = 255; % TUM_STR: Pink；
        end
    end
end
imwrite(img,['./img_new/', maskList(i).name(1:end-4), '.tiff']);
subplot(1,2,2),imshow(img);
