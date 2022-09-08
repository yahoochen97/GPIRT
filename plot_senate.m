% read results
addpath("./results");
file_name = 'all_nominate_data_90.csv';
opts = detectImportOptions(file_name);
all_nominate = readtable(file_name);
all_nominate = sortrows(all_nominate, ["party"]);

% plot three ideology scores
% 2001-2003, 2007-2009 , 2019-2021
session_ids = 107:116;
session_ids = 92:101;

fig = figure(1);
tiledlayout(2,5);
for i=1:numel(session_ids)
    nexttile;
    session_id = session_ids(i);
    x = all_nominate.nominate(all_nominate.session==session_id);
    y = all_nominate.gpirt(all_nominate.session==session_id);
    party = all_nominate.party(all_nominate.session==session_id);
   
    if numel(unique(party))==3
        colors = 'bkr';
    else
        colors = 'br';
    end
    h = gscatter(x,y,party, colors,'s^o', 8);
    for n = 1:length(h)
      set(h(n), 'MarkerFaceColor', colors(n));
    end
    xlim([-1.0,1.0]);
    ylim([-1.2,1.2]);
    xlabel('NOMINATE Dimension 1 Ideology','FontSize', 16);
    ylabel('GPIRT Ideology','FontSize', 16);
    xticks([-0.5,0,0.5]);
    title( int2str(session_id) + "th U.S. Congress", 'FontSize', 12);
    legend('Location','northwest','FontSize',12);
    legend boxoff;
    
end
set(fig, 'PaperPosition', [0 0 30 10]); 
set(fig, 'PaperSize', [30 10]); 

filename = "./results/senate_nominate_90.pdf";
print(fig, filename, '-dpdf','-r300');
close;

file_name = 'gpirt_abortion_results.csv';
opts = detectImportOptions(file_name);
all_nominate = readtable(file_name);
all_nominate = sortrows(all_nominate, ["party"]);

% plot three ideology scores
% 92 to 110 sessions
session_ids = 92:2:110;

fig = figure(1);
tiledlayout(2,5);
for i=1:numel(session_ids)
    nexttile;
    session_id = session_ids(i);
    x = all_nominate.nominate(all_nominate.session==session_id);
    y = all_nominate.gpirt(all_nominate.session==session_id);
    party = all_nominate.party(all_nominate.session==session_id);
   
    if numel(unique(party))==3
        colors = 'bkr';
    else
        colors = 'br';
    end
    h = gscatter(x,y,party, colors,'s^o', 8);
    for n = 1:length(h)
      set(h(n), 'MarkerFaceColor', colors(n));
    end
    xlim([-1.0,1.0]);
    ylim([-0.8,0.8]);
    xlabel('NOMINATE Dimension 1 Ideology','FontSize', 16);
    ylabel('GPIRT Ideology','FontSize', 16);
    xticks([-0.5,0,0.5]);
    title( int2str(session_id) + "th U.S. Congress", 'FontSize', 12);
    legend('Location','northwest','FontSize',12);
    legend boxoff;
    
end
set(fig, 'PaperPosition', [0 0 30 10]); 
set(fig, 'PaperSize', [30 10]); 

filename = "./results/abortion_nominate.pdf";
print(fig, filename, '-dpdf','-r300');
close;

file_name = 'gpirt_CIRI_results.csv';
opts = detectImportOptions(file_name);
all_nominate = readtable(file_name);

session_ids = 1981:3:2010;

fig = figure(1);
tiledlayout(2,5);
for i=1:numel(session_ids)
    nexttile;
    session_id = session_ids(i);
    x = all_nominate.CIRI_theta(all_nominate.session==session_id);
    y = all_nominate.gpirt(all_nominate.session==session_id);
    y = sign(corr(x,y))*y;
    continents = all_nominate.continent(all_nominate.session==session_id);
    colors = 'gbrkc';
    h = gscatter(x,y,continents, colors,'s^odv', 6);
    for n = 1:length(h)
      set(h(n), 'MarkerFaceColor', colors(n));
    end
    xlim([-3.0,3.0]);
    ylim([-1.5,1.5]);
    xlabel('CIRI score','FontSize', 16);
    ylabel('GPIRT score','FontSize', 16);
    xticks((-3.0):0.5:3.0);
    title("Y " + int2str(session_id), 'FontSize', 12);
    legend('Location','southeast','FontSize',12);
    legend boxoff;
    
end
set(fig, 'PaperPosition', [0 0 30 10]); 
set(fig, 'PaperSize', [30 10]); 

filename = "./results/gpirt_CIRI_results.pdf";
print(fig, filename, '-dpdf','-r300');
close;


% file_name = 'gpirt_icc_result.csv';
% opts = detectImportOptions(file_name);
% all_icc = readtable(file_name);
% all_icc.session = str2double(all_icc.session);
% all_icc.nominate = str2double(all_icc.nominate);
% all_icc.lik = str2double(all_icc.lik);
% 
% hs = [108, 109, 110, 111, 112, 113];
% hs = [91, 92, 93, 94, 95, 96];
% 
% fig = figure(2);
% tiledlayout(2, 3);
% for i=1:numel(hs)
%     nexttile;
%     session_id = hs(i);
%     x = all_icc.nominate(all_icc.session==session_id);
%     y = all_icc.lik(all_icc.session==session_id);
% 
%     plot(x,y, 'LineWidth', 4);
%     xlim([-1.,1.0]);
%     ylim([0,1]);
%     xticks([-1.0,-0.5,0,0.5,1.0]);
%     xlabel('GPIRT Ideology','FontSize', 16);
%     ylabel('Pr(y = 1)','FontSize', 16);
%     title( int2str(session_id) + "th U.S. Congress", 'FontSize', 12);
% end
% set(fig, 'PaperPosition', [0 0 18 10]); 
% set(fig, 'PaperSize', [18 10]); 
% 
% filename = "./results/senate_irf.pdf";
% print(fig, filename, '-dpdf','-r300');
% close;
