addpath("./results");
file_name = 'TAPS_holdout_summary.csv';
% file_name = 'TAPS_holdout_summary_rebuttal.csv';
opts = detectImportOptions(file_name);
results = readtable(file_name);
METRICS = ["test_acc","test_lls"];
MODELS = unique(results.model);
DRS = unique(results.dropratio);
HORIZON = max(results.horizon);

% results = results(results.dropratio==10,:);

fig = figure(1);
tiledlayout(1,numel(METRICS),'Padding', 'none', 'TileSpacing', 'compact');   

for i=1:numel(MODELS)
    for k=1:numel(METRICS)
%     nexttile;
    disp(METRICS(k));
        disp(MODELS(i));
        y = zeros(numel(DRS),HORIZON);
        yerr = zeros(numel(DRS),HORIZON);
        tmp = results(strcmp(results.metric,METRICS(k)) & strcmp(results.model,MODELS(i)),:);
        for j=1:numel(DRS)
            for h=1:HORIZON
                tmp1 = tmp.v(tmp.horizon==h & tmp.dropratio==DRS(j),:);
                y(j,h) = mean(tmp1(tmp1~=0)); 
                yerr(j,h) = std(tmp1(tmp1~=0));
            end
        end
%         errorbar(1:HORIZON,y,yerr); hold on;
%         errorbar([1,5,10], y(4,[1,5,10]),yerr(4,[1,5,10])); hold on;
%         disp(y(4,[1,5,10]));
%         disp(yerr(4,[1,5,10])/5);
        disp(y(1,[1,5,10]));
        disp(yerr(1,[1,5,10])/5);
    end
    
%     xlabel('Congress session','FontSize', 18);
%     ylabel(char(score_types(k)),'FontSize', 12);
%     ylabel(YLABELS(k), 'FontSize', 18);
%     ylim([-1.0*(3-k),1.0*(3-k)]);
    
%     legend({"DO-IRT","GD-GPIRT"}, 'Location','southeast','FontSize',14, 'NumColumns',1);
%     legend boxoff;
end
    
% set(fig, 'PaperPosition', [0 0 20 4]); 
% set(fig, 'PaperSize', [20 4]); 
% 
% filename = "./results/abortion_dynamic.pdf";
% print(fig, filename, '-dpdf','-r300','-fillpage');
% close;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

METRICS = ["train_acc","train_lls"];
MODELS = unique(results.model);
DRS = unique(results.dropratio);
HORIZON = max(results.horizon);

for i=1:numel(MODELS)
    for k=1:numel(METRICS)
        disp(METRICS(k));
        disp(MODELS(i));
        y = zeros(numel(DRS),HORIZON);
        yerr = zeros(numel(DRS),HORIZON);
        tmp = results(strcmp(results.metric,METRICS(k)) & strcmp(results.model,MODELS(i)),:);
        for j=1:numel(DRS)
            h = 0;
            tmp1 = tmp.v(tmp.horizon==h & tmp.dropratio==DRS(j),:);
        end
        disp( mean(tmp1(tmp1~=0)));
%         disp( std(tmp1(tmp1~=0)/5));
    end
    
end