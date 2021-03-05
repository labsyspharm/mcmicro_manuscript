%% Plot CyCIF mask on CODEX ROI
% Load mask for CyCIF

%% CyCIF
mask_cycif = ...
    double(imread('/Users/denisschapiro/Desktop/TNP_Data/CyCIF/New_seg/cellMask.tif'));
mask_cycif_old = ...
    double(imread('/Users/denisschapiro/Desktop/TNP_Data/CyCIF/TNP_pilot_cycif/segmentation/cellMask.tif'));

mask_raw_ROI_DAPI = double(imread('/Users/denisschapiro/Desktop/TNP_Data/CyCIF/New_seg/cellMask.tif',1,...
    'PixelRegion',{[17100,20600],[8900,11700]}));

border= double(boundarymask(mask_raw_ROI));
imwrite(border,'ROI_mask_CYCIF_border.tif')


image_raw_ROI_DAPI = imread('/Users/denisschapiro/Desktop/TNP_Data/CyCIF/registration/TNP_pilot_cycif.ome.tif',1,...
    'PixelRegion',{[17100,20600],[8900,11700]});
image_raw_ROI_Keratin = imread('/Users/denisschapiro/Desktop/TNP_Data/CyCIF/registration/TNP_pilot_cycif.ome.tif',12,...
    'PixelRegion',{[17100,20600],[8900,11700]});
image_raw_ROI_CD8 = imread('/Users/denisschapiro/Desktop/TNP_Data/CyCIF/registration/TNP_pilot_cycif.ome.tif',15,...
    'PixelRegion',{[17100,20600],[8900,11700]});
image_raw_ROI_CD20 = imread('/Users/denisschapiro/Desktop/TNP_Data/CyCIF/registration/TNP_pilot_cycif.ome.tif',19,...
    'PixelRegion',{[17100,20600],[8900,11700]});

image_mask_vis = image_raw_ROI_DAPI;
image_mask_vis(:,:,2) = image_raw_ROI_Keratin;
image_mask_vis(:,:,3) = image_raw_ROI_CD8;
image_mask_vis(:,:,4) = image_raw_ROI_CD20;


imwrite(image_raw_ROI,'image_ROI_DAPI.tif')

for i=1:4
imwrite(image_mask_vis(:,:,i), 'image_ROI_all4.tif', 'writemode', 'append');
end

figure;imshow(image_ROI);
imwrite(image_ROI,'napari_vis_cycif_DAPI.tif');

% % Select intersect:
% mask_cycif_ROI = mask_cycif(17100:20600,8900:11700);
% figure;imshow(mask_cycif_ROI);

% Select intersect:
mask_cycif_old_ROI = mask_cycif_old(17100:20600,8900:11700);
figure;imshow(mask_cycif_old_ROI);

% Create RGB
%image =zeros(size(mask_cycif_ROI,1),size(mask_cycif_ROI,2),3);
image =ones(size(mask_cycif_ROI,1),size(mask_cycif_ROI,2),3);

% Get cell ids for CD20
cell_ids_CD8_ROI = intersect(unique(mask_cycif_old_ROI),CyCIF.CellID(CyCIF.CD8_pos ==1));
cell_ids_CD20_ROI = intersect(unique(mask_cycif_old_ROI),CyCIF.CellID(CyCIF.CD20_pos == 1));
cell_ids_keratin_ROI = intersect(unique(mask_cycif_old_ROI),CyCIF.CellID(CyCIF.Keratin_pos == 1));

for i=1:size(cell_ids_CD8_ROI,1)
    row = []; col = [];
    [row,col] = find(mask_cycif_old_ROI==cell_ids_CD8_ROI(i));
    for k=1:size(row,1)
        image(row(k),col(k),1) = 65535;
        image(row(k),col(k),2) = 0;
        image(row(k),col(k),3) = 0;
    end
end

for i=1:size(cell_ids_keratin_ROI,1)
    row = []; col = [];
    [row,col] = find(mask_cycif_old_ROI==cell_ids_keratin_ROI(i));
    for k=1:size(row,1)
        image(row(k),col(k),1) = 0;
        image(row(k),col(k),2) = 65535;
        image(row(k),col(k),3) = 0;
    end
end

for i=1:size(cell_ids_CD20_ROI,1)
    row = []; col = [];
    [row,col] = find(mask_cycif_old_ROI==cell_ids_CD20_ROI(i));
    for k=1:size(row,1)
        image(row(k),col(k),1) = 0;
        image(row(k),col(k),2) = 0;
        image(row(k),col(k),3) = 65535;
    end
end

figure, imshow(image)
%imwrite(image,'napari_vis_cycif.tif')
imwrite(image,'cycif_mask.tif')


%% Input CODEX
mask_path = '/Users/denisschapiro/Desktop/TNP_Data/CyCIF/TNP_pilot_cycif/segmentation/cellMask.tif';
ROI = 1;
ROI_inputX = 17100:20600;
ROI_inputY = 8900:11700;
table_CellID = CyCIF.CellID;
table_CD8 = CyCIF.CD8_pos;
table_CD20 = CyCIF.CD20_pos;
table_Keratin = CyCIF.Keratin_pos;
name_output = 'cycif_mask.tif';

% Run CODEX plot
PlotMask(mask_path,ROI,ROI_inputX,ROI_inputY,table_CellID,...
    table_CD8,table_CD20,table_Keratin,name_output)


%% CODEX
mask_codex =...
    double(imread('/Users/denisschapiro/Desktop/TNP_Data/TNP_pilot_codex/segmentation/PilotTonsil_5_z08/cellMask.tif'));

% Create RGB
image =ones(size(mask_codex,1),size(mask_codex,2),3);

% Get cell ids for CD20
cell_ids_CD8_ROI = intersect(unique(mask_codex),CyCIF.CellID(CODEX.CD8_pos ==1));
cell_ids_CD20_ROI = intersect(unique(mask_codex),CyCIF.CellID(CODEX.CD20_pos == 1));
cell_ids_keratin_ROI = intersect(unique(mask_codex),CyCIF.CellID(CODEX.Keratin_pos == 1));

for i=1:size(cell_ids_CD8_ROI,1)
    row = []; col = [];
    [row,col] = find(mask_cycif_old_ROI==cell_ids_CD8_ROI(i));
    for k=1:size(row,1)
        image(row(k),col(k),1) = 65535;
        image(row(k),col(k),2) = 0;
        image(row(k),col(k),3) = 0;
    end
end

for i=1:size(cell_ids_keratin_ROI,1)
    row = []; col = [];
    [row,col] = find(mask_cycif_old_ROI==cell_ids_keratin_ROI(i));
    for k=1:size(row,1)
        image(row(k),col(k),1) = 0;
        image(row(k),col(k),2) = 65535;
        image(row(k),col(k),3) = 0;
    end
end

for i=1:size(cell_ids_CD20_ROI,1)
    row = []; col = [];
    [row,col] = find(mask_cycif_old_ROI==cell_ids_CD20_ROI(i));
    for k=1:size(row,1)
        image(row(k),col(k),1) = 0;
        image(row(k),col(k),2) = 0;
        image(row(k),col(k),3) = 65535;
    end
end

figure, imshow(image)
%imwrite(image,'napari_vis_cycif.tif')
imwrite(image,'cycif_mask.tif')



%% mcmicro_compare_segmentation --archive
addpath(pwd)
ilastik = readtable('ilastik-PilotTonsil_5_z08.csv');
unmicst = readtable('unmicst-PilotTonsil_5_z08.csv');

scatter(ilastik.X_centroid,ilastik.Y_centroid);
hold on
scatter(unmicst.X_centroid,unmicst.Y_centroid);
hold off

ilastik

% Merge tables
ilastik_matrix = [table2array(ilastik),ones(size(ilastik_matrix,1),1)];
unmicst_matrix = [table2array(unmicst),2*ones(size(unmicst,1),1)];

merge = [ilastik_matrix;unmicst_matrix];
merge_table = array2table(merge,'VariableNames',ilastik.Properties.VariableNames);