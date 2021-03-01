%mcmicro_compare_segmentation
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