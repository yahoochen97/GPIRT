addpath("results");
file_name = 'SupremeCourt_holdout_results.csv';
opts = detectImportOptions(file_name);
data = readtable(file_name);
data = data(strcmp(data.type, "0")==0,:);
data = sortrows(data, ["model"]);
measures = unique(data.measure);
MODELS = unique(data.model);
colors = ["b","r","k"];

fig = figure(1);
tiledlayout(1,3,'Padding', 'none', 'TileSpacing', 'compact');
for k=1:numel(measures)
   nexttile;
   measure = measures{k};
   for p=1:numel(MODELS)
      MODEL = MODELS{p};
      tmp = data(strcmp(data.measure, measure) & strcmp(data.model, MODEL),:);
      means = str2double(tmp.mean);
      sds = str2double(tmp.sd);
      horizon = str2double(tmp.type);
      errorbar(horizon+0.1*(p-2), means, sds,'square', ...
          'MarkerSize',8,...
          'MarkerEdgeColor',colors(p),...
          'MarkerFaceColor',colors(p),...
          'LineWidth', 3); hold on;
   end
   xlim([0.0,max(horizon)+1]);
   xticks(1:max(horizon));
   xlabel('Forecast horizon (years)','FontSize', 12);
   ylabel("Averaged " + measure,'FontSize', 12);
   legend(upper(MODELS), 'Location','north','FontSize',12);
   legend boxoff;
end

set(fig, 'PaperSize', [15 5]); 

filename = "./results/SupremeCourt_holdout_compare.pdf";
print(fig, filename, '-dpdf','-r300', '-fillpage');
close;

