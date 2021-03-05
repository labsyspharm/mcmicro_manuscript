function [] = PlotMask(mask_path,ROI,ROI_inputX,ROI_inputY,table_CellID,...
    table_CD8,table_CD20,table_Keratin,name_output)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%% Input
mask_path = '/Users/denisschapiro/Desktop/TNP_Data/CyCIF/TNP_pilot_cycif/segmentation/cellMask.tif';
ROI = 1;
ROI_inputX = 17100:20600;
ROI_inputY = 8900:11700;
table_CellID = CyCIF.CellID;
table_CD8 = CyCIF.CD8_pos;
table_CD20 = CyCIF.CD20_pos;
table_Keratin = CyCIF.Keratin_pos;
name_output = 'cycif_mask.tif';

%% 
mask_full = ...
    double(imread(mask_path));

% Select intersect
if ROI == 1
    mask_ROI = mask_full(ROI_inputX,ROI_inputY);
    figure;imshow(mask_ROI);
else
    mask_ROI = mask_full
end

% Create RGB
image_empty =ones(size(mask_ROI,1),size(mask_ROI,2),3);

% Get cell ids for CD20
cell_ids_CD8_ROI = intersect(unique(mask_ROI),table_CellID(table_CD8 ==1));
cell_ids_CD20_ROI = intersect(unique(mask_ROI),table_CellID(table_CD20 == 1));
cell_ids_keratin_ROI = intersect(unique(mask_ROI),table_CellID(table_Keratin == 1));

for i=1:size(cell_ids_CD8_ROI,1)
    row = []; col = [];
    [row,col] = find(mask_ROI==cell_ids_CD8_ROI(i));
    for k=1:size(row,1)
        image(row(k),col(k),1) = 65535;
        image(row(k),col(k),2) = 0;
        image(row(k),col(k),3) = 0;
    end
end

for i=1:size(cell_ids_keratin_ROI,1)
    row = []; col = [];
    [row,col] = find(mask_ROI==cell_ids_keratin_ROI(i));
    for k=1:size(row,1)
        image(row(k),col(k),1) = 0;
        image(row(k),col(k),2) = 65535;
        image(row(k),col(k),3) = 0;
    end
end

for i=1:size(cell_ids_CD20_ROI,1)
    row = []; col = [];
    [row,col] = find(mask_ROI==cell_ids_CD20_ROI(i));
    for k=1:size(row,1)
        image(row(k),col(k),1) = 0;
        image(row(k),col(k),2) = 0;
        image(row(k),col(k),3) = 65535;
    end
end

figure, imshow(image)
imwrite(image,name_output)
end

