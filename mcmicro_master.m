%% Load CODEX, CyCIF and mIHC
CODEX = readtable('PilotTonsil_5_z08.csv');
CyCIF = readtable('TNP_pilot_cycif.csv');
mIHC = readtable('unmicst-mIHC_ROI.csv');

%% Extract relevant ROIs
% CODEX (WARNING!!! - old X and Y results)
figure()
scatter(CODEX.Y_position,-CODEX.X_position,[],log(CODEX.CD8));

% CyCIF (WARNING!!! - old X and Y results)
% scatter(CyCIF.Y_position,-CyCIF.X_position,[],log(CyCIF.Keratin_570));
% Manual selection
intersect_new_cycif = find(CyCIF.Y_position > 8900 & CyCIF.Y_position < 11700 & ...
    CyCIF.X_position > 17100 & CyCIF.X_position < 20600);
figure()
scatter(-CyCIF.X_position(intersect_new_cycif),...
    -CyCIF.Y_position(intersect_new_cycif),[],...
log(CyCIF.CD8a_488(intersect_new_cycif)));

% mIHC
%figure()
%scatter(mIHC.Y_centroid,mIHC.X_centroid,[],log(mIHC.CD8_cellMask));
intersect_new_mIHC = find(mIHC.Y_centroid> 1500 & mIHC.Y_centroid < 5500 & ...
    mIHC.X_centroid > 2500 & mIHC.X_centroid < 6000);
figure()
scatter(mIHC.Y_centroid(intersect_new_mIHC),mIHC.X_centroid(intersect_new_mIHC)...
    ,[],log(mIHC.CD8_cellMask(intersect_new_mIHC)));

%% Create tSNE
% Select markers of interest
tSNE_codex_data =   [CODEX.CD20,CODEX.PanCK,CODEX.CD8];
tSNE_cycif_data =   [CyCIF.CD20_488,CyCIF.Keratin_570,CyCIF.CD8a_488];
tSNE_mIHC_data =    [mIHC.CD20_cellMask,mIHC.PANCK_cellMask,mIHC.CD8_cellMask];

% Select ROIs 
tSNE_cycif_data_ROI = tSNE_cycif_data(intersect_new_cycif,:);
tSNE_mIHC_data_ROI = tSNE_mIHC_data(intersect_new_mIHC,:);

% Combine to run tSNE function and for visualization 
tSNE_combine = [tSNE_codex_data;tSNE_cycif_data_ROI;tSNE_mIHC_data_ROI];
tSNE_code = [zeros(size(tSNE_codex_data,1),1);ones(size(tSNE_cycif_data_ROI,1),1);2*ones(size(tSNE_mIHC_data_ROI,1),1)];

%% Run tSNE
[tsne_xy,y] = tsne(tSNE_combine);

% Run tSNE normalized
tSNE_codex_data_norm = normalize(tSNE_codex_data);
tSNE_cycif_data_ROI_norm = normalize(tSNE_cycif_data_ROI);
tSNE_mIHC_data_ROI_norm = normalize(tSNE_mIHC_data_ROI);

tSNE_combine_norm = [tSNE_codex_data_norm;tSNE_cycif_data_ROI_norm;tSNE_mIHC_data_ROI_norm];

subsample_rand = randperm(size(tSNE_combine_norm,1));

[tsne_xy_norm,y_norm] = tsne(tSNE_combine_norm(subsample_rand(1:10000),1:3));

%% Visualize tSNE
figure()
scatter(tsne_xy(find(tSNE_code==0),1),tsne_xy(find(tSNE_code==0),2),[],[0, 0.4470, 0.7410]);
hold on
scatter(tsne_xy(find(tSNE_code==1),1),tsne_xy(find(tSNE_code==1),2),[],[0.4940, 0.1840, 0.5560]);
hold on
scatter(tsne_xy(find(tSNE_code==2),1),tsne_xy(find(tSNE_code==2),2),[],[0.25, 0.25, 0.25]);
hold off

% Normalized
tSNE_code_subsample = tSNE_code(subsample_rand(1:10000));
figure()
scatter(tsne_xy_norm(find(tSNE_code_subsample==0),1),tsne_xy_norm(find(tSNE_code_subsample==0),2),[],[0, 0.4470, 0.7410]);
hold on
scatter(tsne_xy_norm(find(tSNE_code_subsample==1),1),tsne_xy_norm(find(tSNE_code_subsample==1),2),[],[0.4940, 0.1840, 0.5560]);
hold on
scatter(tsne_xy_norm(find(tSNE_code_subsample==2),1),tsne_xy_norm(find(tSNE_code_subsample==2),2),[],[0, 0, 0]);
hold off

%% Manual cell type
% Function call
[CODEX,CyCIF,mIHC] = Manual_thresholding(CODEX,CyCIF,mIHC,intersect_new_cycif,intersect_new_mIHC);

tSNE_code_CD20 = [CODEX.CD20_pos == 1;CyCIF.CD20_pos(intersect_new_cycif)==1;...
    mIHC.CD20_pos(intersect_new_mIHC)==1];

tSNE_code_CD20_subsample = tSNE_code_CD20(subsample_rand(1:10000));

tSNE_code_Keratin = [CODEX.Keratin_pos == 1;CyCIF.Keratin_pos(intersect_new_cycif) ==1;...
    mIHC.Keratin_pos(intersect_new_mIHC)==1];

tSNE_code_Keratin_subsample = tSNE_code_Keratin(subsample_rand(1:10000));

tSNE_code_CD8 = [CODEX.CD8_pos == 1;CyCIF.CD8_pos(intersect_new_cycif) ==1;...
    mIHC.CD8_pos(intersect_new_mIHC)==1];

tSNE_code_CD8_subsample = tSNE_code_CD8(subsample_rand(1:10000));

% Visualize cell types on tSNE
figure()
scatter(tsne_xy(:,1),tsne_xy(:,2),[],[0.8 0.8 0.8]);
hold on
scatter(tsne_xy(find(tSNE_code_CD20),1),tsne_xy(find(tSNE_code_CD20),2),[],[0.9290 0.6940 0.1250]);
hold on
scatter(tsne_xy(find(tSNE_code_Keratin),1),tsne_xy(find(tSNE_code_Keratin),2),[],[0.4660 0.6740 0.1880]);
hold on
scatter(tsne_xy(find(tSNE_code_CD8),1),tsne_xy(find(tSNE_code_CD8),2),[],[0.6350 0.0780 0.1840]);
hold off

figure()
scatter(tsne_xy_norm(:,1),tsne_xy_norm(:,2),[],[0.8 0.8 0.8]);
hold on
scatter(tsne_xy_norm(find(tSNE_code_CD20_subsample),1),tsne_xy_norm(find(tSNE_code_CD20_subsample),2),[],[0.9290 0.6940 0.1250]);
hold on
scatter(tsne_xy_norm(find(tSNE_code_Keratin_subsample),1),tsne_xy_norm(find(tSNE_code_Keratin_subsample),2),[],[0.4660 0.6740 0.1880]);
hold on
scatter(tsne_xy_norm(find(tSNE_code_CD8_subsample),1),tsne_xy_norm(find(tSNE_code_CD8_subsample),2),[],[0.6350 0.0780 0.1840]);
hold off

% Plot cell types
% figure()
% % CD20
% CD20_index_cycif = intersect(intersect_new_cycif,find(CyCIF.CD20_pos == 1));
% scatter(-CyCIF.X_position(CD20_index_cycif),...
%     -CyCIF.Y_position(CD20_index_cycif),...
%     [],'r');
% hold on
% 
% % Keratin
% Kerain_index_cycif = intersect(intersect_new_cycif,find(CyCIF.Keratin_pos == 1));
% scatter(-CyCIF.X_position(Kerain_index_cycif),...
%     -CyCIF.Y_position(Kerain_index_cycif),...
%     [],'g');
% hold on
% 
% % CD8
% CD8_index_cycif = intersect(intersect_new_cycif,find(CyCIF.CD8_pos == 1));
% scatter(-CyCIF.X_position(CD8_index_cycif),...
%     -CyCIF.Y_position(CD8_index_cycif),...
%     [],'b');
% hold off
% 
% 
% % CODEX
% % CD20 CODEX
% % --> Let's use 8 as cut-off
% CODEX.CD20_pos = (log(CODEX.CD20)>8);
% % Keratin CODEX
% % --> Let's use 7.5 as cut-off
% CODEX.Keratin_pos = (log(CODEX.PanCK)>7.5);
% % CD8 CODEX
% % --> Let's use 6.5 as cut-off
% CODEX.CD8_pos = (log(CODEX.CD8)>6.5);
% 
% % Plot cell types
% figure()
% % CD31
% scatter(CODEX.Y_position(CODEX.CD31_pos == 1),...
%     -CODEX.X_position(CODEX.CD31_pos == 1),...
%     [],'r');
% hold on
% 
% % Keratin
% scatter(CODEX.Y_position(CODEX.Keratin_pos == 1),...
%     -CODEX.X_position(CODEX.Keratin_pos == 1),...
%     [],'g');
% hold on
% 
% % CD8
% scatter(CODEX.Y_position(CODEX.CD8_pos == 1),...
%     -CODEX.X_position(CODEX.CD8_pos == 1),...
%     [],'b');
% hold off
% 
% % mIHC
% % CD31 mIHC
% mIHC.CD8_pos = (log(mIHC.CD8_cellMask)>4);
% % CD8
% scatter(mIHC.Y_centroid(mIHC.CD8_pos == 1),...
%     -mIHC.X_centroid(mIHC.CD8_pos == 1),...
%     [],'b');


%% Calculate distribution
%CODEX
CODEX_CD20_cells = size(find(CODEX.CD20_pos == 1),1);
CODEX_Keratin_cells = size(find(CODEX.Keratin_pos== 1),1);
CODEX_CD8_cells = size(find(CODEX.CD8_pos == 1),1);
%CODEX_CD45_cells = size(find(CODEX.CD45_pos == 1),1);

figure()
X = categorical({'CD20','Keratin','CD8'});
X = reordercats(X,{'CD20','Keratin','CD8'});
Y = [CODEX_CD20_cells,CODEX_Keratin_cells,CODEX_CD8_cells];
bar(X,Y);title('CODEX')

%CyCIF
CyCIF_CD20_cells = size(find(CyCIF.CD20_pos == 1),1);
CyCIF_Keratin_cells = size(find(CyCIF.Keratin_pos== 1),1);
CyCIF_CD8_cells = size(find(CyCIF.CD8_pos == 1),1);
%CyCIF_CD45_cells = size(find(CyCIF.CD45_pos == 1),1);

figure()
X = categorical({'CD20','Keratin','CD8'});
X = reordercats(X,{'CD20','Keratin','CD8'});
Y = [CyCIF_CD20_cells,CyCIF_Keratin_cells,CyCIF_CD8_cells];
bar(X,Y);title('CYCIF')

%mIHC
mIHC_CD20_cells = size(find(mIHC.CD20_pos == 1),1);
mIHC_Keratin_cells = size(find(mIHC.Keratin_pos== 1),1);
mIHC_CD8_cells = size(find(mIHC.CD8_pos == 1),1);
%mIHC_CD45_cells = size(find(mIHC.CD45_pos == 1),1);

figure()
X = categorical({'CD20','Keratin','CD8'});
X = reordercats(X,{'CD20','Keratin','CD8'});
Y = [mIHC_CD20_cells,mIHC_Keratin_cells,mIHC_CD8_cells];
bar(X,Y);title('mIHC')

figure()
X = categorical({'CD20','Keratin','CD8'});
X = reordercats(X,{'CD20','Keratin','CD8'});

% Calculate ratio for visualization
all_cells_CyCIF = CyCIF_CD20_cells+CyCIF_Keratin_cells+CyCIF_CD8_cells;
all_cells_CODEX = CODEX_CD20_cells+CODEX_Keratin_cells+CODEX_CD8_cells;
all_cells_mIHC = mIHC_CD20_cells+mIHC_Keratin_cells+mIHC_CD8_cells;



Y = [CyCIF_CD20_cells/all_cells_CyCIF,CODEX_CD20_cells/all_cells_CODEX,...
    mIHC_CD20_cells/all_cells_mIHC;...
    CyCIF_Keratin_cells/all_cells_CyCIF,CODEX_Keratin_cells/all_cells_CODEX,...
    mIHC_Keratin_cells/all_cells_mIHC;...
    CyCIF_CD8_cells/all_cells_CyCIF,CODEX_CD8_cells/all_cells_CODEX,...
    mIHC_CD8_cells/all_cells_mIHC];

bar(X,Y);






%% Backlog
CODEX_CD8_norm = norm_vector(CODEX.CD8);
CODEX_CD31_norm = norm_vector(CODEX.CD31);

CyCIF_CD8_norm = norm_vector(CyCIF.CD8a_488);
CyCIF_CD31_norm = norm_vector(CyCIF.CD31_647);

scatter(CODEX_CD8_norm,CODEX_CD31_norm);
scatter(CyCIF_CD8_norm,CyCIF_CD31_norm);

figure()
histogram(CODEX_CD8_norm,'FaceColor','blue','FaceAlpha',0.1)
set(gca,'xscale','log')
hold on
histogram(CyCIF_CD8_norm,'FaceColor','red','FaceAlpha',0.1)
set(gca,'xscale','log')
hold off

figure()
histogram(CODEX_CD31_norm,'FaceColor','blue','FaceAlpha',0.1)
set(gca,'xscale','log')
hold on
histogram(CyCIF_CD31_norm,'FaceColor','red','FaceAlpha',0.1)
set(gca,'xscale','log')
hold off

