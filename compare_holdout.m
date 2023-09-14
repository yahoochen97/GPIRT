addpath("results");
file_name = 'SupremeCourt_holdout_results.csv';
opts = detectImportOptions(file_name);
data = readtable(file_name);
data = data(strcmp(data.type, "0")==0,:);
data = sortrows(data, ["model"]);
measures = unique(data.measure);
MODELS = unique(data.model);
MODELS = {'doirt','SEirt','gpirt'};
colors = {[160 160 160]/255,[204 102 0]/255,[255 51 255]/255};
LABELS = ["Pred accuracy", "Log likelihood", "Posterior std"];
measures = {'acc';'sd'};
LABELS = ["Pred accuracy of votes", "Post std of ideology"];

fig = figure(1);
tiledlayout(1,numel(measures),'Padding', 'none', 'TileSpacing', 'compact');
for k=1:numel(measures)
   nexttile;
   measure = measures{k};
   for p=1:numel(MODELS)
      MODEL = MODELS{p};
      tmp = data(strcmp(data.measure, measure) & strcmp(data.model, MODEL),:);
      means = str2double(tmp.mean);
      sds = str2double(tmp.sd);
      horizon = str2double(tmp.type);
      errorbar(horizon+0.2*(p-2), means, sds,'square', ...
          'Color', colors{p}, ...
          'MarkerSize',8,...
          'MarkerEdgeColor',colors{p},...
          'MarkerFaceColor',colors{p},...
          'LineWidth', 3); hold on;
%       h = bar(horizon+0.2*(p-2), means, 'BarWidth', 0.1);
%       set(h,'FaceColor',);
   end
   xlim([0.0,max(horizon)+1]);
   xticks(1:max(horizon));
   xlabel('Forecast horizon (years)','FontSize', 16);
   ylabel(LABELS(k),'FontSize', 16);
   if k==1
    legend({'D-IRT','GD-GPIRT(RBF)','GD-GPIRT(Matérn)'}, 'Location','north','FontSize',12);
    legend boxoff;
   end
end

set(fig, 'PaperSize', [10 5]); 

filename = "./results/SupremeCourt_holdout_compare.pdf";
print(fig, filename, '-dpdf','-r300', '-fillpage');
close;

