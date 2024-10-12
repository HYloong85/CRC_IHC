clc;
clear;
close all;
openslide_load_library();

load('vgg19_model.mat');

addpath(genpath('./subroutines/'));
mkdir('./output_images');
mkdir('./output_masks');
mkdir('./output_rgbout');
mkdir('./output_allstats');
mkdir('./output_margin');

datapath = ('./data/');
% xmlList = dir([datapath, 'TMA*.xml']);
svsList = dir([datapath, 'TMA*.svs']);
N = length(svsList);

colors = [255, 255, 255; 255, 0 , 0; 255, 255, 0; 0, 0, 255]/255; % white, red, yellow, blue
tissuenames = {'Disturb', 'GLA', 'STR', 'TUM'};
marginnames = {'Left', 'Down', 'Right'};
verbose = true;
border = [76 76]; 
                   
layerNum = numel(myNet.Layers);
bsize = myNet.Layers(1).InputSize;
bsize = bsize(1:2)-2*border;
rmov = ceil(border(1)/bsize(1)); 

for i = 1:N
    t0 = tic;
    t = 1;
    svs = openslide_open([datapath, svsList(i).name]);

%% If there are multiple tissue specimens annotated with bounding boxes in an image, record the coordinates of each bounding box.
%     xmlDoc = xml_read([datapath, xmlList(i).name]);
%     %Number of annotation boxes in an SVS and record their coordinates;£»
%     anno = length(xmlDoc.Annotation.Regions.Region);
%     pts = zeros(anno, 4);
%     for j = 1:anno
%         pts(j, :) = [xmlDoc.Annotation.Regions.Region(j).Vertices.Vertex(1).ATTRIBUTE.X,...
%             xmlDoc.Annotation.Regions.Region(j).Vertices.Vertex(1).ATTRIBUTE.Y,...
%             xmlDoc.Annotation.Regions.Region(j).Vertices.Vertex(3).ATTRIBUTE.X,...
%             xmlDoc.Annotation.Regions.Region(j).Vertices.Vertex(3).ATTRIBUTE.Y];
%     end
%     
%     missed_tissue_in = missed_tissue;
%     if str2double(svsList(i).name(4:5)) == 1
%         missed_tissue_in(13, :) = -1;
%     end
%     if str2double(svsList(i).name(4:5)) == 8
%         missed_tissue_in(2, :) = missed_tissue_in(3, :);
%     end
%     t = 0;    
%     for j = 1:anno
%         t1 = tic;
%         
%         if ( mod(i, 15)==0) && (j==1||j==10||j==11||j==12) 
%             continue;
%         end
%         t = t+1;
%         
%         if missed_tissue_in(str2double(svsList(i).name(6:7)), j) < 0
%             allstats(t, :) = -1;
%             margin(t, :) = -1;
%             fprintf('Skipped %s, tissue %d;\n', svsList(i).name(1:end-4), t);
%             continue;
%         end
%         
%         x = round(pts(j, 1));
%         y = round(pts(j, 2));
%         width = round(pts(j, 3) - pts(j, 1));
%         height = round(pts(j, 4) - pts(j, 2));
    [width, height] = openslide_get_level_dimensions(svs,0);
    image = openslide_read_region(svs, 0, 0, width, height, 0);
    image = image(:, :, 2:4);
    
    ind = floor(2000/bsize(1));
    
    max_col = floor(double(width) / ((ind - 3) * bsize(1)));
    max_row = floor(double(height) / ((ind - 3) * bsize(1)));
    nPat = (max_col-1) * (max_row-1);
    
    mask_total = cell(1, nPat);
    count1 = 0;
    
    %% Exclude whether it exceeds the boundary
    subwidth_error = width - ((ind - 3) * bsize(1)) * (max_col - 2);
    right2_error = floor(double(subwidth_error)/bsize(1))-rmov;
    subheight_error = height - ((ind - 3) * bsize(1)) * (max_row - 2);
    right1_error = floor(double(subheight_error)/bsize(1))-rmov;
    
    if ((right2_error+1)*bsize(1)+border(1) > subwidth_error) || ((right1_error+1)*bsize(1)+border(1) > subheight_error)
    error(i,j) = 1;
    allstats(t, :) = 0;
    margin(t, :) = 0;
    fprintf('Errored %s, tissue %d;\n', svsList(i).name(1:end-4), t);
    continue;
    end

    %% Use a sliding window to select a specified size of patch for prediction in each tissue microarray
    for row = 1:max_row-1
        for col = 1:max_col-1
            
            count1 = count1 + 1;
            
            subheight = 2000;
            subwidth = 2000;
            left = rmov;
            right1 = ind-rmov;
            right2 = ind-rmov;
            
            if col == max_col - 1 && row ~= max_row - 1
                subwidth = width - ((ind - 3) * bsize(1)) * (max_col - 2);
                right2 = floor(double(subwidth)/bsize(1))-rmov;
            end
            if row == max_row - 1 && col ~= max_col - 1
                subheight = height - ((ind - 3) * bsize(1)) * (max_row - 2);
                right1 = floor(double(subheight)/bsize(1))-rmov;
            end
            if row == max_row - 1 && col == max_col - 1
                subwidth = width - ((ind - 3) * bsize(1)) * (max_col - 2);
                subheight = height - ((ind - 3) * bsize(1)) * (max_row - 2);
                right1 = floor(double(subheight)/bsize(1))-rmov;
                right2 = floor(double(subwidth)/bsize(1))-rmov;
            end
                result = image((row-1) * ((ind - 3) * bsize(1))+1:(row-1) * ((ind - 3) * bsize(1))+subheight,...
                    (col-1) * ((ind - 3) * bsize(1))+1:(col-1) * ((ind - 3) * bsize(1))+subwidth,:);
            
            count2 = 0;
            imgset = [];
            for p = left:right1
                for q = left:right2
                    count2 = count2 + 1;
                    imgset(:, :, :, count2) = result(p*bsize(1)-border(1)+1:(p+1)*bsize(1)+border(1),...
                        q*bsize(1)-border(1)+1:(q+1)*bsize(1)+border(1), :);
                end
            end
            
            [labels, score] = classify(myNet, imgset);
            
            mid1 = zeros((right1-1)*(right2-1), 1, 4);
            mid3 = zeros(right1-1, right2-1, 4);
            for g = 1:4
                mid1(:, :, g) = score(:, g);
            end
            mid2 = reshape(mid1, right2-1, right1-1, 4);
            for k = 1:4
                mid3(:, :, k) = mid2(:, :, k)';
            end
            mask_total{count1} = mid3;
        end
    end
        
    %% Combine patches in sequence to form the original image
    mask = cell(max_row-1, max_col-1);
    count3 = 0;
    for m = 1:max_row-1
        for n = 1:max_col-1
            count3 = count3 + 1;
            mask{m, n} = mask_total{count3};
        end
    end
    mask_out = cell2mat(mask);

    [rgbout, currstats] = mask9toRGB(mask_out, colors);
    allstats(t, :) = currstats(:);

    % Record the edges of the image that are not covered by the sliding window
    margin_left = left * bsize(1) - border(1);
    margin_down = subheight - ((right1 + 1) * bsize(1) + border(1));
    margin_right = subwidth - ((right2 + 1) * bsize(1) + border(1));
    margin(t, :) = [margin_left, margin_down, margin_right];
    
    save (['./output_masks/mask_out_', svsList(i).name(1:end-4),...
        '_tissue_', num2str(t), '.mat'], 'mask_out');
    save (['./output_rgbout/rgbout_', svsList(i).name(1:end-4),...
        '_tissue_', num2str(t), '.mat'], 'rgbout');
    imwrite(rgbout, ['./output_images/rgbout_',svsList(i).name(1:end-4),...
        '_tissue_', num2str(t), '.tiff']);
    
    fprintf('Successfully finished %s, tissue %d;\n',...
        svsList(i).name(1:end-4), t);
    
%     end

    statstable1 = array2table(allstats, 'VariableNames', tissuenames);
    writetable(statstable1, ['./output_allstats/allstats_',...
        svsList(i).name(1:end-4),'.xlsx']);
    
    statstable2 = array2table(margin, 'VariableNames', marginnames);
    writetable(statstable2, ['./output_margin/margin_',...
        svsList(i).name(1:end-4),'.xlsx']);
    
    fprintf('Successfully finished %s, time %f;\n',...
        svsList(i).name(1:end-4), toc(t0));
    
end