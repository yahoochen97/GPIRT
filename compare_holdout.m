addpath("results");
file_name = 'SupremeCourt_holdout.csv';
opts = detectImportOptions(file_name);
data = readtable(file_name);
measures = unique(data.measure);

fig = figure(1);
tiledlayout(1,3,'Padding', 'none', 'TileSpacing', 'compact');
for k=1:numel(measures)
   nexttile;
   
   
end


