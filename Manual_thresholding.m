function [CODEX,CyCIF,mIHC] = Manual_thresholding(CODEX,CyCIF,mIHC,intersect_new_cycif,intersect_new_mIHC)
%UNTITLED5 Summary of this function goes here

% Scatter size
size_scatter = 0.7;

%% CD20, PanCK, CD8 and CD45
% CD20
%histogram(log(tSNE_codex_data(:,1)));
CODEX.CD20_pos = (log(CODEX.CD20)>7.9);
CODEX_CD20_pos_index = find(CODEX.CD20_pos);

% Image = zeros(size(mask_codex,1),size(mask_codex,2),3);
%
% Image(mask_codex == 10);
%
% cell_index = (find(mask_codex == 10));
% [row,col] = ind2sub(size(mask_codex),cell_index);

figure()
scatter(CODEX.Y_position(CODEX.CD20_pos == 1),...
    -CODEX.X_position(CODEX.CD20_pos == 1),...
    size_scatter,'b');
hold on

%PanCK
%histogram(log(tSNE_codex_data(:,2)));
CODEX.Keratin_pos = (log(CODEX.PanCK)>7.5);

scatter(CODEX.Y_position(CODEX.Keratin_pos == 1),...
    -CODEX.X_position(CODEX.Keratin_pos == 1),...
    size_scatter,'g');
hold on

%CD8
%histogram(log(tSNE_codex_data(:,3)));
CODEX.CD8_pos = (log(CODEX.CD8)>6.2);

scatter(CODEX.Y_position(CODEX.CD8_pos == 1),...
    -CODEX.X_position(CODEX.CD8_pos == 1),...
    size_scatter,'r');
hold on

%CD45
%histogram(log(tSNE_codex_data(:,4)));
% CODEX.CD45_pos = (log(CODEX.CD45)>6.1);
%
% scatter(CODEX.Y_position(CODEX.CD45_pos == 1),...
%     -CODEX.X_position(CODEX.CD45_pos == 1),...
%     [],'c');
hold off

%% CD20, PanCK, CD8 and CD45
% CD20
%histogram(log(tSNE_cycif_data(:,1)));
CyCIF.CD20_pos = (log(CyCIF.CD20_488)>10.0);
CD20_index_cycif = intersect(intersect_new_cycif,find(CyCIF.CD20_pos == 1));

figure()
scatter(-CyCIF.X_position(CD20_index_cycif),...
    -CyCIF.Y_position(CD20_index_cycif),...
    size_scatter,'b');
hold on

%PanCK
%histogram(log(tSNE_cycif_data(:,2)));
CyCIF.Keratin_pos = (log(CyCIF.Keratin_570)>7.5);
PanCK_index_cycif = intersect(intersect_new_cycif,find(CyCIF.Keratin_pos == 1));

scatter(-CyCIF.X_position(PanCK_index_cycif),...
    -CyCIF.Y_position(PanCK_index_cycif),...
    size_scatter,'g');
hold on

%CD8
%histogram(log(tSNE_cycif_data(:,3)));
CyCIF.CD8_pos = (log(CyCIF.CD8a_488)>9);

CD8_index_cycif = intersect(intersect_new_cycif,find(CyCIF.CD8_pos == 1));

scatter(-CyCIF.X_position(CD8_index_cycif),...
    -CyCIF.Y_position(CD8_index_cycif),...
    size_scatter,'r');
hold on

%CD45
%histogram(log(tSNE_cycif_data(:,4)));
% CyCIF.CD45_pos = (log(CyCIF.CD45_647)>9.5);
%
% CD45_index_cycif = intersect(intersect_new_cycif,find(CyCIF.CD45_pos == 1));
%
% scatter(-CyCIF.X_position(CD45_index_cycif),...
%     -CyCIF.Y_position(CD45_index_cycif),...
%     [],'c');
hold off

%% CD20, PanCK, CD8 and CD45
% CD20
%histogram(log(tSNE_mIHC_data(:,1)));
mIHC.CD20_pos = (log(mIHC.CD20_cellMask)>4.2);
CD20_index_mIHC = intersect(intersect_new_mIHC,find(mIHC.CD20_pos == 1));

figure()
scatter(mIHC.Y_centroid(CD20_index_mIHC),...
    mIHC.X_centroid(CD20_index_mIHC),...
    size_scatter,'b');
hold on

%PanCK
%histogram(log(tSNE_mIHC_data(:,2)));
mIHC.Keratin_pos = (log(mIHC.PANCK_cellMask)>3.5);
PanCK_index_mIHC = intersect(intersect_new_mIHC,find(mIHC.Keratin_pos == 1));

scatter(mIHC.Y_centroid(PanCK_index_mIHC),...
    mIHC.X_centroid(PanCK_index_mIHC),...
    size_scatter,'g');
hold on

%CD8
%histogram(log(tSNE_mIHC_data(:,3)));
mIHC.CD8_pos= (log(mIHC.CD8_cellMask)>4.4);
CD8_index_mIHC = intersect(intersect_new_mIHC,find(mIHC.CD8_pos == 1));

scatter(mIHC.Y_centroid(CD8_index_mIHC),...
    mIHC.X_centroid(CD8_index_mIHC),...
    size_scatter,'r');
hold on

%CD45
%histogram(log(tSNE_mIHC_data(:,4)));
% mIHC.CD45_pos= (log(mIHC.CD45_cellMask)>4.7);
% CD45_index_mIHC = intersect(intersect_new_mIHC,find(mIHC.CD45_pos == 1));
%
% scatter(mIHC.Y_centroid(CD45_index_mIHC),...
%     mIHC.X_centroid(CD45_index_mIHC),...
%     [],'c');
hold off

end

