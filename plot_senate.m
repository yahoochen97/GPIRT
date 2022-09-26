% read results
addpath("./results");
file_name = 'all_nominate_data_90.csv';
opts = detectImportOptions(file_name);
all_nominate = readtable(file_name);
all_nominate = sortrows(all_nominate, ["party"]);

% plot three ideology scores
% 2001-2003, 2007-2009 , 2019-2021
session_ids = 90:99;

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
        shapes = 's^o';
    else
        colors = 'br';
        shapes = 'so';
    end
    h = gscatter(x,y,party, colors,shapes, 8);
    for n = 1:length(h)
      set(h(n), 'MarkerFaceColor', colors(n));
    end
    xlim([-1.0,1.0]);
    ylim([-3.0,3.0]);
    xlabel('D-NOMINATE score','FontSize', 16);
    ylabel('GD-GPIRT score','FontSize', 16);
    xticks([-1.0:0.5:1.0]);
    yticks([-3.0:0.5:3.0]);
    title( int2str(session_id) + "th U.S. Congress", 'FontSize', 12);
    legend('Location','northwest','FontSize',12);
    legend boxoff;
end

set(fig, 'PaperPosition', [0 0 30 10]); 
set(fig, 'PaperSize', [30 10]); 

filename = "./results/senate_nominate_90.pdf";
print(fig, filename, '-dpdf','-r300');
close;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

file_name = 'gpirt_abortion_dynamic.csv';
opts = detectImportOptions(file_name);
all_nominate = readtable(file_name);
all_nominate = sortrows(all_nominate, ["type"]);

% plot three ideology scores
% 92 to 110 sessions
session_ids = 92:110;

senator_ids = unique(all_nominate.icpsr);
senator_names = [];
for icpsr = senator_ids'
    senator_names = [senator_names,string(unique(all_nominate.bioname(all_nominate.icpsr==icpsr,1)))];
end

colors = ["b","k","r","g","magenta"];
score_types = unique(all_nominate.type);
shapes = ["-s", "-o"];

fig = figure(1);
tiledlayout(1,2);   
for k=1:2
    clear h;
    nexttile;
    for i=1:numel(senator_ids)
        senator_id = senator_ids(i);
        x = all_nominate.session(all_nominate.icpsr==senator_id & strcmp(all_nominate.type, score_types(k)));
        y = all_nominate.score(all_nominate.icpsr==senator_id & strcmp(all_nominate.type, score_types(k)));
        h{i} = plot(x,y,strcat(shapes(k),colors(i)),'MarkerSize',8); hold on;
    end
    
    xlabel('Congress session','FontSize', 12);
    ylim([-1.0*(3-k),1.0*(3-k)]);
    ylabel(char(score_types(k)),'FontSize', 12);
    legend(senator_names, 'Location','southeast','FontSize',6);
    legend boxoff;
end
    
set(fig, 'PaperPosition', [0 0 20 5]); 
set(fig, 'PaperSize', [20 5]); 

filename = "./results/abortion_dynamic.pdf";
print(fig, filename, '-dpdf','-r300');
close;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% file_name = 'gpirt_abortion_results.csv';
% opts = detectImportOptions(file_name);
% all_nominate = readtable(file_name);
% all_nominate = sortrows(all_nominate, ["party"]);
% icpsr_ids = unique(all_nominate.icpsr);
% % plot three ideology scores
% % 92 to 110 sessions
% session_ids = 92:1:110;
% score_types = ["nominate","gpirt"];
% shapes = ["b-s", "b-o"];
% 
% all_nominate.gpirt = all_nominate.gpirt ./ std(all_nominate.gpirt);
% all_nominate.nominate = all_nominate.nominate ./ std(all_nominate.nominate);
% 
% select_icpsr = [];
% 
% fig = figure(1);
% tiledlayout(2,1);   
% for k=1:1
%     clear h;
%     nexttile;
%     for i=1:numel(icpsr_ids)
%         icpsr = icpsr_ids(i);
%         x = all_nominate.session(all_nominate.icpsr==icpsr);
%         y = all_nominate.(score_types(k))(all_nominate.icpsr==icpsr);
%         bioname = unique(all_nominate.bioname(all_nominate.icpsr==icpsr));
%         if(std(y)>=0.2), select_icpsr=[select_icpsr,icpsr];end
%         plot(x,y,shapes(k),'MarkerSize',1); hold on;
%     end
%     
%     if k==2, xlabel('Congress session','FontSize', 16); end
%     ylim([-3.0,3.0]);
%     ylabel(char(score_types(k)),'FontSize', 16);
%     hleg = legend('Location','southwest','FontSize',8);
%     set(hleg,'visible','off')
% end

% fig = figure(1);
% tiledlayout(6,6);
% for icpsr=select_icpsr
%     nexttile;
%     x = all_nominate.session(all_nominate.icpsr==icpsr);
%     y1 = all_nominate.(score_types(1))(all_nominate.icpsr==icpsr);
%     y2 = all_nominate.(score_types(2))(all_nominate.icpsr==icpsr);
%     bioname = unique(all_nominate.bioname(all_nominate.icpsr==icpsr));
%     plot(x,y1,'b');title(bioname);
%     legend(score_types(1),'Location','best','FontSize',8);
% end
%     
% set(fig, 'PaperPosition', [0 0 15 15]); 
% set(fig, 'PaperSize', [15 15]); 
% 
% filename = "./results/abortion_dynamic.pdf";
% print(fig, filename, '-dpdf','-r300');
% close;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

file_name = 'gpirt_abortion_results.csv';
opts = detectImportOptions(file_name);
all_nominate = readtable(file_name);
all_nominate = sortrows(all_nominate, ["party"]);
% plot three ideology scores
% 92 to 110 sessions
session_ids = 92:4:110;

fig = figure(1);
tiledlayout(1,5);
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
    
    text(0.6,-2.8,"{\rho} = " + round(corr(x,y),3), 'FontSize',14);
    for n = 1:length(h)
      set(h(n), 'MarkerFaceColor', colors(n));
    end
    xlim([-1.0,1.0]);
    ylim([-3.2,3.2]);
    xlabel('D-NOMINATE Ideology','FontSize', 16);
    ylabel('GD-GPIRT Ideology','FontSize', 16);
    xticks([-1.0, -0.5,0,0.5, 1.0]);
    title( int2str(session_id) + "th U.S. Congress", 'FontSize', 12);
    legend('Location','northwest','FontSize',12);
    legend boxoff;
    
end
set(fig, 'PaperPosition', [0 0 30 5]); 
set(fig, 'PaperSize', [30 5]); 

filename = "./results/abortion_nominate.pdf";
print(fig, filename, '-dpdf','-r300');
close;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

file_name = 'gpirt_CIRI_results.csv';
opts = detectImportOptions(file_name);
all_nominate = readtable(file_name);

session_ids = 1981:6:2010;

fig = figure(1);
tiledlayout(2,5);
for i=1:numel(session_ids)
    nexttile;
    session_id = session_ids(i);
    x = all_nominate.CIRI_theta(all_nominate.session==session_id);
    y = all_nominate.gpirt(all_nominate.session==session_id);
    % y = sign(corr(x,y))*y;
    % y = y./ std(y);
    continents = all_nominate.continent(all_nominate.session==session_id);
    colors = 'gbrkc';
    h = gscatter(x,y,continents, colors,'s^odv', 6);
    for n = 1:length(h)
      set(h(n), 'MarkerFaceColor', colors(n));
    end
    xlim([-3.5,3.5]);
    ylim([-3.5,3.5]);
    xlabel('DO-IRT score','FontSize', 16);
    ylabel('GD-GPIRT score','FontSize', 16);
    xticks((-3.5):0.5:3.5);
    title("Y " + int2str(session_id) + " (all)", 'FontSize', 12);
    legend('Location','southeast','FontSize',12);
    legend boxoff;
    
end

for i=1:numel(session_ids)
    nexttile;
    session_id = session_ids(i);
    x = all_nominate.CIRI_theta(all_nominate.session==session_id);
    y = all_nominate.gpirt(all_nominate.session==session_id);
    % y = sign(corr(x,y))*y;
    % y = y./ std(y);
    continents = all_nominate.continent(all_nominate.session==session_id);
    colors = 'gbrkc';
    h = gscatter(x,y,continents, colors,'s^odv', 6);
    for n = 1:length(h)
      set(h(n), 'MarkerFaceColor', colors(n));
    end
    xlim([-1.0,1.0]);
    ylim([-1.0,1.0]);
    xlabel('DO-IRT score','FontSize', 16);
    ylabel('GD-GPIRT score','FontSize', 16);
    xticks((-1.0):0.5:1.0);
    title("Y " + int2str(session_id)  + " (zoomed in)", 'FontSize', 12);
    legend('Location','southeast','FontSize',12);
    legend boxoff;
    
end

set(fig, 'PaperPosition', [0 0 30 10]); 
set(fig, 'PaperSize', [30 10]); 

filename = "./results/gpirt_CIRI_results.pdf";
print(fig, filename, '-dpdf','-r300');
close;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

file_name = 'gpirt_CIRI_dynamic.csv';
opts = detectImportOptions(file_name);
all_nominate = readtable(file_name);
all_nominate = sortrows(all_nominate, ["type"]);

% plot three ideology scores
% 92 to 110 sessions
session_ids = 1981:2010;

country_ids = unique(all_nominate.id);
country_names = [];
for id = country_ids'
    country_names = [country_names,
        string(unique(all_nominate.country(all_nominate.id==id,1)))];
end

colors = ["k", "b","r","magenta"];
score_types = flip(unique(all_nominate.type));
shapes = ["-s", "-o"];

fig = figure(1);
tiledlayout(1,2);   
for k=1:2
    clear h;
    nexttile;
    for i=1:numel(country_ids)
        id = country_ids(i);
        x = all_nominate.session(all_nominate.id==id & strcmp(all_nominate.type, score_types(k)));
        y = all_nominate.score(all_nominate.id==id & strcmp(all_nominate.type, score_types(k)));
        h{i} = plot(x,y,strcat(shapes(k),colors(i)),'MarkerSize',8); hold on;
    end
    
    xlabel('year','FontSize', 12);
    % ylim([-1.0*(3-k),1.0*(3-k)]);
    ylim([-2,1.0]);
    ylabel(char(score_types(k)),'FontSize', 12);
    legend(country_names, 'Location','northwest','FontSize',8);
    legend boxoff;
end
    
set(fig, 'PaperPosition', [0 0 20 5]); 
set(fig, 'PaperSize', [20 5]); 

filename = "./results/CIRI_dynamic.pdf";
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
