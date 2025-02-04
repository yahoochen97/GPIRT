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
tiledlayout(2,5,'Padding', 'none', 'TileSpacing', 'compact');
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
print(fig, filename, '-dpdf','-r300', '-fillpage');
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
YLABELS = ["GD-GPIRT", "D-NOMINATE"];
shapes = ["-s", "-o"];
% senator_names(1) = "R. Byrd (D)";
% senator_names(2) = "E. Kennedy (D)";
% senator_names(3) = "L. Weicker (R)";
% senator_names(4) = "J. Helms (R)";
% senator_names(5) = "C. Grassley (R)";
% senator_names(1) = "Q. Burdick (D)";
% senator_names(2) = "R. Byrd (D)";
% senator_names(3) = "E. Kennedy (D)";
% senator_names(4) = "C. Pell (D)";
% senator_names(5) = "L. Weicker (R)";
% senator_names(6) = "J. Biden (D)";
% senator_names(7) = "J. Helms (R)";
% senator_names(8) = "C. Grassley (R)";
% senator_names(9) = "M. Wallop (R)";
% senator_names(10) = "G. Mitchell (D)";
% senator_names(11) = "T. Gorton (R)";
senator_names(1) = "R. Byrd (D)";
senator_names(2) = "E. Kennedy (D)";
senator_names(3) = "J. Biden (D)";
senator_names(4) = "O. Hatch (R)";
senator_names(5) = "M. Wallop (R)";
senator_names(6) = "G. Mitchell (D)";


% fig = figure(1);
% tiledlayout(1,2,'Padding', 'none', 'TileSpacing', 'compact');   
% for k=1:2
%     clear h;
%     nexttile;
%     for i=1:numel(senator_ids)
%         senator_id = senator_ids(i);
%         x = all_nominate.session(all_nominate.icpsr==senator_id & strcmp(all_nominate.type, score_types(k)));
%         y = all_nominate.score(all_nominate.icpsr==senator_id & strcmp(all_nominate.type, score_types(k)));
%         h{i} = plot(x,y,strcat(shapes(k),colors(i)),'MarkerSize',8, 'LineWidth', 4); hold on;
%     end
%     
%     xlabel('Congress session','FontSize', 18);
%     ylabel(char(score_types(k)),'FontSize', 12);
%     ylabel(YLABELS(k), 'FontSize', 18);
%     ylim([-1.0*(3-k),1.0*(3-k)]);
%     
%     legend(senator_names, 'Location','southeast','FontSize',14, 'NumColumns',5);
%     legend boxoff;
% end
%     
% set(fig, 'PaperPosition', [0 0 20 4]); 
% set(fig, 'PaperSize', [20 4]); 
% 
% filename = "./results/abortion_dynamic.pdf";
% print(fig, filename, '-dpdf','-r300','-fillpage');
% close;


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
tiledlayout(1,5,'Padding', 'none', 'TileSpacing', 'compact');
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
print(fig, filename, '-dpdf','-r300', '-fillpage');
close;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

session_ids = 92:2:110;
fig = figure(1);
% subplot(4,1,1);
h1 = subplot(2,1,1);
set(h1, 'OuterPosition', [0,0.75, 1, .25]);
set(h1, 'Position', [0.05,0.80, 0.95, .20]);
% tiledlayout(4,1,'Padding', 'none', 'TileSpacing', 'compact');
polarization = [];
% nexttile;
for i=1:numel(session_ids)
    session_id = session_ids(i);
    x = all_nominate.gpirt(all_nominate.session==session_id);
    party = all_nominate.party(all_nominate.session==session_id);
    polarization = [polarization, sum(strcmp(party,'Republicans')==(x>=0))/numel(x)];
end
p1 = plot((1:1:numel(session_ids))-1,polarization,'-',...
    'Color',[0.5,0.5,0.5], 'LineWidth',2);

ylim([0.4,1.0]);
YTICKS = (1:1:numel(session_ids))-1;
xticks(YTICKS);
xlim([-0.6,(numel(session_ids))+0.2]);
xtickformat('%dth');
% xticklabels(session_ids);
a = get(gca,'XTickLabel');
set(gca,'box','off','ytick',[0.5:0.25:1.0],...
    'XTickLabel',(session_ids-min(session_ids))*2+1971,'fontsize',8);
lgd = legend([p1],{'Polarization ratio'},...
 'Location','northwest','FontSize',12, 'NumColumns' ,1);
legend boxoff;
lgd.Position(2) = 0.960;

% nexttile([3 1]);
% subplot(4,1,[2,3,4]);
h2 = subplot(2,1,2);
set(h2, 'OuterPosition', [0, 0, 1, .70]);
set(h2, 'Position', [0.05,0.05, 0.95, .65]);
for i=1:numel(session_ids)
    session_id = session_ids(i);
    x = all_nominate.gpirt(all_nominate.session==session_id);
    party = all_nominate.party(all_nominate.session==session_id);
    polarization = [polarization, sum(strcmp(party,'Republicans')==(x>=0))/numel(x)];
   
    [f,Xi] = ksdensity(x(strcmp(party,'Democrats')));
    f = f(abs(Xi)<3);
    Xi = Xi(abs(Xi)<3);
    MARGIN = 0.025;
%     plot(f+i-1+MARGIN,Xi, 'Color', 'b');
    p1 = fill(i-1-f-MARGIN,Xi,'b','FaceAlpha',0.5,'LineStyle','none');
    hold on;
    [f,Xi] = ksdensity(x(strcmp(party,'Republicans')));
    f = f(abs(Xi)<3);
    Xi = Xi(abs(Xi)<3);
%     plot(i-1-f-MARGIN,Xi,'Color', 'r');
    p2 = fill(i-1+f+MARGIN,Xi,'r','FaceAlpha',0.5,'LineStyle','none');
%     text(i-1-0.1,(2*polarization(end)-1)*3+0.2,sprintf('%.2f',round(polarization(end),2)), 'FontSize',20);     
end

% yyaxis left;
ylim([-3.0,3.8]);
ylabel('GD-GPIRT score','FontSize', 12);
ylabel('Pro-life','FontSize', 16);
% yticks((-3.0):1.0:3.0);
YTICKS = (1:1:numel(session_ids))-1;
xticks(YTICKS);
xlim([-0.6,(numel(session_ids))+0.2]);
xtickformat('%dth');
% xticklabels(session_ids);
a = get(gca,'XTickLabel');
set(gca,'box','off','ytick',[],...
    'XTickLabel',(session_ids-min(session_ids))*2+1971,'fontsize',8);
% annotation('textarrow', [0 0],[0.8 1.0], 'String','pro-life');
% [-1.0,2.0],...
%     'YTickLabel',["pro-choice","pro-life"],'fontsize',8,...
%     'YTickLabelRotation',90,...

% the arrows
annotation('arrow', [0.05 0.05],[0.6 0.7], 'LineStyle', 'none');

% unique(all_nominate.icpsr)'

select_icpsr = [1366, 10808,12032,14105,14226];
select_icpsr = [1366, 10808, 14101, 14503, 14511, 14713];

i=1;
xadjust = [-0.8, 0.0, 0.0,0, 0.8,1.0 ];
yadjust = [-0.25, -0.2, 0.25,0.2,0.6,-0.6 ];
for icpsr=select_icpsr
    x = all_nominate.gpirt(all_nominate.icpsr==icpsr & any(all_nominate.session==session_ids,2));
    y = (all_nominate.session(all_nominate.icpsr==icpsr & any(all_nominate.session==session_ids,2)));
    y = (y-min(session_ids))/2;
    x_std = std(x);
    if(x_std>0.0)
        p4 = line(y,x,'LineStyle','--','Color',[0.3,0.3,0.3], 'LineWidth',2);
        hold on;
        idx = int8(numel(x)/2);
        text(y(idx)+xadjust(i),x(idx)+yadjust(i),senator_names(i), 'FontSize',10);
        i = i + 1;
    end
end

% title( int2str(session_id) + "th U.S. Congress", 'FontSize', 12);
lgd = legend([p1,p2],{'Democrats','Republicans'},...
    'Location','southwest','FontSize',12, 'NumColumns' ,2);
legend boxoff;
lgd.Position(1) = 0.03;
lgd.Position(2) = 0.650;

set(fig, 'PaperPosition', [0 0 10 4]); 
set(fig, 'PaperSize', [10 4]); 

filename = "./results/abortion_gpirt.pdf";
print(fig, filename, '-dpdf','-r300', '-fillpage');
close;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
session_ids = 92:2:110;
fig = figure(1);
% tiledlayout(4,1,'Padding', 'none', 'TileSpacing', 'compact');
% nexttile;
h1 = subplot(2,1,1);
set(h1, 'OuterPosition', [0,0.75, 1, .25]);
set(h1, 'Position', [0.05,0.80, 0.95, .20]);
polarization = [];
for i=1:numel(session_ids)
    session_id = session_ids(i);
    x = all_nominate.nominate(all_nominate.session==session_id);
    party = all_nominate.party(all_nominate.session==session_id);
    polarization = [polarization, sum(strcmp(party,'Republicans')==(x>=0))/numel(x)];
end
p1 = plot((1:1:numel(session_ids))-1,polarization,'-','Color',...
    [0.5,0.5,0.5], 'LineWidth',2);

ylim([0.8,1.0]);
YTICKS = (1:1:numel(session_ids))-1;
xticks(YTICKS);
xlim([-0.6,(numel(session_ids))+0.2]);
xtickformat('%dth');
% xticklabels(session_ids);
a = get(gca,'XTickLabel');
set(gca,'box','off','ytick',[0.8:0.1:1.0],...
    'XTickLabel',(session_ids-min(session_ids))*2+1971,'fontsize',8);
lgd = legend({'Polarization ratio'},'Location','northwest','FontSize',12, 'NumColumns' ,1);
legend boxoff;
lgd.Position(2) = 0.960;

% nexttile([3 1]);
h2 = subplot(2,1,2);
set(h2, 'OuterPosition', [0, 0, 1, .70]);
set(h2, 'Position', [0.05,0.05, 0.95, .65]);
polarization = [];
for i=1:numel(session_ids)
    session_id = session_ids(i);
    x = all_nominate.nominate(all_nominate.session==session_id);
    party = all_nominate.party(all_nominate.session==session_id);
    polarization = [polarization, sum(strcmp(party,'Republicans')==(x>=0))/numel(x)];
   
    [f,Xi] = ksdensity(x(strcmp(party,'Democrats')));
    f = f(abs(Xi)<3)/4;
    Xi = Xi(abs(Xi)<3);
    MARGIN = 0.025;
%     plot(f+i-1+MARGIN,Xi, 'Color', 'b');
    p1 = fill(i-1-f-MARGIN,Xi,'b','FaceAlpha',0.5,'LineStyle','none');
    hold on;
    [f,Xi] = ksdensity(x(strcmp(party,'Republicans')));
    f = f(abs(Xi)<1)/4;
    Xi = Xi(abs(Xi)<1);
%     plot(i-1-f-MARGIN,Xi,'Color', 'r');
    p2 = fill(i-1+f+MARGIN,Xi,'r','FaceAlpha',0.5,'LineStyle','none');
%     text(i-1-0.1,(2*polarization(end)-1)*1+0.1,sprintf('%.2f',round(polarization(end),2)), 'FontSize',20);     
end

% yyaxis left;
ylim([-1.0,1.2]);
ylabel('Pro-life','FontSize', 16);
% yticks((-3.0):1.0:3.0);
YTICKS = (1:1:numel(session_ids))-1;
xticks(YTICKS);
xlim([-0.8,(numel(session_ids))+0.2]);
xtickformat('%dth');
% xticklabels(session_ids);
a = get(gca,'XTickLabel');
set(gca,'box','off','ytick',[],...
    'XTickLabel',(session_ids-min(session_ids))*2+1971,'fontsize',8);

% the arrows
annotation('arrow', [0.05 0.05],[0.6 0.7], 'LineStyle', 'none');

% yyaxis right;
% yticks(YTICKS);
% ylim([-0.2,(numel(session_ids)-0.8)]);
% yticklabels(flip(round(polarization,2)));
% ylabel('Polarization ratio','FontSize', 16);
% set(get(gca,'ylabel'),'Rotation',270, 'Position', [4.3,4.5],'Color', 'k');
% set(get(gca,'YTickLabel'),'Color', 'k');

% unique(all_nominate.icpsr)'

i=1;
xadjust = [-0.5, 0.0, -0.1, 0.1,0,0.1];
yadjust = [0.1, -0.15, 0.15,0.1,0.15,-0.15];
for icpsr=select_icpsr
    x = all_nominate.nominate(all_nominate.icpsr==icpsr & any(all_nominate.session==session_ids,2));
    y = (all_nominate.session(all_nominate.icpsr==icpsr & any(all_nominate.session==session_ids,2)));
    y = (y-min(session_ids))/2;
    x_std = std(x);
    if(x_std>0.0)
        p4 = line(y,x,'LineStyle','--','Color',[0.3,0.3,0.3], 'LineWidth',2);
        hold on;
        idx = int8(numel(x)/2);
        if icpsr<=10808, idx = 1;end
        text(y(idx)+xadjust(i),x(idx)+yadjust(i),senator_names(i), 'FontSize',10);
        i = i + 1;
    end
end

% title( int2str(session_id) + "th U.S. Congress", 'FontSize', 12);
lgd = legend([p1,p2],{'Democrats','Republicans'},...
    'Location','northwest','FontSize',12, 'NumColumns' ,2);
legend boxoff;
lgd.Position(2) = 0.650;

set(fig, 'PaperPosition', [0 0 10 4]); 
set(fig, 'PaperSize', [10 4]); 

filename = "./results/abortion_nominate.pdf";
print(fig, filename, '-dpdf','-r300', '-fillpage');
close;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% file_name = 'gpirt_CIRI_results.csv';
% opts = detectImportOptions(file_name);
% all_nominate = readtable(file_name);
% 
% session_ids = 1981:2010;
% CIRI_cor = [];
% for i=1:numel(session_ids)
%     session_id = session_ids(i);
%     x = all_nominate.CIRI_theta(all_nominate.session==session_id);
%     y = all_nominate.gpirt(all_nominate.session==session_id);
%     CIRI_cor = [CIRI_cor,corr(x,y)]; 
% end
% 
% fig = figure(1);
% tiledlayout(2,5);
% 
% session_ids = 1981:6:2010;
% for i=1:numel(session_ids)
%     nexttile;
%     session_id = session_ids(i);
%     x = all_nominate.CIRI_theta(all_nominate.session==session_id);
%     y = all_nominate.gpirt(all_nominate.session==session_id);
%     % y = sign(corr(x,y))*y;
%     % y = y./ std(y);
%     continents = all_nominate.continent(all_nominate.session==session_id);
%     colors = 'gbrkc';
%     h = gscatter(x,y,continents, colors,'s^odv', 6);
%     for n = 1:length(h)
%       set(h(n), 'MarkerFaceColor', colors(n));
%     end
%     xlim([-3.5,3.5]);
%     ylim([-3.5,3.5]);
%     xlabel('DO-IRT score','FontSize', 16);
%     ylabel('GD-GPIRT score','FontSize', 16);
%     xticks((-3.5):0.5:3.5);
%     title("Y " + int2str(session_id) + " (all)", 'FontSize', 12);
%     legend('Location','southeast','FontSize',12);
%     legend boxoff;
%     
% end
% 
% for i=1:numel(session_ids)
%     nexttile;
%     session_id = session_ids(i);
%     x = all_nominate.CIRI_theta(all_nominate.session==session_id);
%     y = all_nominate.gpirt(all_nominate.session==session_id);
%     % y = sign(corr(x,y))*y;
%     % y = y./ std(y);
%     continents = all_nominate.continent(all_nominate.session==session_id);
%     colors = 'gbrkc';
%     h = gscatter(x,y,continents, colors,'s^odv', 6);
%     for n = 1:length(h)
%       set(h(n), 'MarkerFaceColor', colors(n));
%     end
%     xlim([-1.0,1.0]);
%     ylim([-1.0,1.0]);
%     xlabel('DO-IRT score','FontSize', 16);
%     ylabel('GD-GPIRT score','FontSize', 16);
%     xticks((-1.0):0.5:1.0);
%     title("Y " + int2str(session_id)  + " (zoomed in)", 'FontSize', 12);
%     legend('Location','southeast','FontSize',12);
%     legend boxoff;
%     
% end
% 
% set(fig, 'PaperPosition', [0 0 30 10]); 
% set(fig, 'PaperSize', [30 10]); 
% 
% filename = "./results/gpirt_CIRI_results.pdf";
% print(fig, filename, '-dpdf','-r300', '-fillpage');
% close;

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

% colors = ["k", "b","r","magenta"];
colors = {[166,206,227]/255,...
    [31,120,180]/255,...
    [178,223,138]/255,...
    [51,160,44]/255};


score_types = flip(unique(all_nominate.type));
YLABELS = ["GD-GPIRT", "DO-IRT"];
shapes = ["-", "-"];

fig = figure(1);
tiledlayout(1,2,'Padding', 'none', 'TileSpacing', 'compact');   
for k=1:2
    clear h;
    ax = nexttile;
    for i=1:numel(country_ids)
        id = country_ids(i);
        x = all_nominate.session(all_nominate.id==id & strcmp(all_nominate.type, score_types(k)));
        y = all_nominate.score(all_nominate.id==id & strcmp(all_nominate.type, score_types(k)));
        y_sd = all_nominate.sd(all_nominate.id==id & strcmp(all_nominate.type, score_types(k)));
        
%         h{i} = errorbar(x,y,y_sd, shapes(k),'Color',colors{i},...
%             'LineWidth', 2); hold on;
        f = [y+2*y_sd; flipdim(y-2*y_sd,1)];
        h{i} = fill([x; flipdim(x,1)], f, colors{i},...
            'facealpha', 0.3, 'edgecolor', 'none'); hold on;
    end
    
%     for i=1:4
%         h{i}.Color = [h{i}.Color 0.7];  % alpha=0.7
%     end
    
%     xlabel('year','FontSize', 10);
    % ylim([-1.0*(3-k),1.0*(3-k)]);
    ylim([-2.0,2.0]);
    ylabel("Respect for human right",'FontSize', 10);
    title(YLABELS(k),'FontSize', 10);
    legend([h{1},h{2},h{3},h{4}],country_names, 'Location','northwest','FontSize',8, 'NumColumns', 4);
    legend boxoff;
    box(ax,'off');
    set(gca,'ytick',[]);
end

annotation('arrow', [0.025 0.025],[0.85 0.95], 'LineStyle', 'none');
annotation('arrow', [0.54 0.54],[0.85 0.95], 'LineStyle', 'none');
    
set(fig, 'PaperPosition', [0 0 10 3]); 
set(fig, 'PaperSize', [10 3]); 

filename = "./results/CIRI_dynamic.pdf";
print(fig, filename, '-dpdf','-r300', '-fillpage');
close;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

file_name = 'gpirt_Supreme_Court_dynamic.csv';
opts = detectImportOptions(file_name);
all_nominate = readtable(file_name);
% all_nominate = sortrows(all_nominate, ["score"]);

close all;
fig = figure(1);
x = all_nominate.score(strcmp(all_nominate.type,"GPIRT"));
y = all_nominate.score(strcmp(all_nominate.type,"Martin-Quinn"));
JUSTICES = all_nominate.name(strcmp(all_nominate.type,"GPIRT"));
[x, idx] = sort(x);
y = y(idx);
JUSTICES = JUSTICES(idx);

n_justice = numel(unique(JUSTICES));
% darkblue 0, 0, 139 dodger blue 30, 144, 255
% orange 255,140,0 darkred 	139, 0, 0
colors = [linspace(0,30,6)',...
    linspace(0,144,6)', ...
    linspace(139,255,6)']/255;

colors = [colors; [linspace(255,139,10)',...
    linspace(140,0,10)', ...
    linspace(0,0,10)']/255];

h = gscatter(y,x,JUSTICES,'k','o',8);
for i=1:n_justice
    h(i).Color=colors(i,:);
    h(i).MarkerFaceColor=colors(i,:);
end
xlim([-4.1,4.2]);
ylim([-2.1,1.8]);
xlabel('Martin-Quinn','FontSize', 30);
ylabel('GD-GPIRT','FontSize', 30);
% xticks([-1.0:0.5:1.0]);
% yticks([-3.0:0.5:3.0]);
% title( int2str(session_id) + "th U.S. Congress", 'FontSize', 12);
legend('Location','northwest','FontSize',20);
legend boxoff;

set(fig, 'PaperPosition', [0 0 10 10]); 
set(fig, 'PaperSize', [10 10]); 

filename = "./results/gpirt_supreme_court_compare.pdf";
print(fig, filename, '-dpdf','-r300', '-fillpage');
close;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compare martin quinn and gd-gpirt in magnitudes of estimates scores

file_name = 'gpirt_Supreme_Court_dynamic.csv';
opts = detectImportOptions(file_name);
all_nominate = readtable(file_name);
% all_nominate = sortrows(all_nominate, ["score"]);

close all;
years = sort(unique(all_nominate.year));
JUSTICES = unique(all_nominate.name(strcmp(all_nominate.type,"GPIRT")));
n_justice = numel(unique(JUSTICES));

slopes_1 = zeros(n_justice,1);
slopes_2 = zeros(n_justice,1);
for i=1:n_justice
    x = all_nominate.score(strcmp(all_nominate.name, JUSTICES(i)) & strcmp(all_nominate.type,"GPIRT"));
    y = all_nominate.score(strcmp(all_nominate.name, JUSTICES(i)) & strcmp(all_nominate.type,"Martin-Quinn"));
    slopes_1(i) = fitlm(1:numel(x),x).Coefficients{2,1};
    slopes_2(i) = fitlm(1:numel(y),y).Coefficients{2,1};
end

[h,p] = ttest(abs(slopes_1),abs(slopes_2));

diff_slopes = mean(abs(slopes_1))-mean(abs(slopes_2));

diff_slopes_std = std(abs(slopes_1)-abs(slopes_2)) / sqrt(numel(slopes_1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TAPS 2014 plot trajectories of 50 confidence and 50 unconfidence ppl
% Add a vertical line for the 2016 election

file_name = 'gpirt_TAPS2014_dynamic.csv';
opts = detectImportOptions(file_name);
all_nominate = readtable(file_name);
waves = sort(unique(all_nominate.wave));

sort_scores_tbl = sortrows(grpstats(all_nominate,...
    "WUSTLID","mean","DataVars","score"),"mean_score");
wustlids = sort_scores_tbl.WUSTLID;

close all;
fig = figure(1);

% 50 least confident

for i=1:50
    tmp = all_nominate(all_nominate.WUSTLID==wustlids(i),:);
%     scatter(waves, tmp.score,8,[205,183,181]/255, 'filled'); hold on;
    plot(waves, tmp.score, 'color', [.5 .5 .5 .3], 'linewidth', 1); hold on;
end
tmp = all_nominate(any(all_nominate.WUSTLID==wustlids(1:50)',2),:);
tmp = groupsummary(tmp, "wave", "mean", "score");
plot(waves, tmp.mean_score, '-','color',[232, 27, 35, 204]/255, 'LineWidth',3); hold on;

% 50 most confident
for i=(numel(wustlids)-49):numel(wustlids)
    tmp = all_nominate(all_nominate.WUSTLID==wustlids(i),:);
%     scatter(waves, tmp.score, 8, [178,58,238]/255, 'filled'); hold on;
    plot(waves, tmp.score, 'color', [.5 .5 .5 .3], 'linewidth', 1); hold on;
end
tmp = all_nominate(any(all_nominate.WUSTLID==wustlids((end-49):end)',2),:);
tmp = groupsummary(tmp, "wave", "mean", "score");
plot(waves, tmp.mean_score, '-','color',[178,58,238,204]/255, 'LineWidth',3); hold on;

% wave=1 is Jan 2014, wave=35 Nov 2016 Election
xline(35.5, '--', 'color',[0.5,0.5,0.5]);
h = text([36], [min(all_nominate.score)], {'Presidential election'});
set(h,'Rotation',90, 'Color', [0.5,0.5,0.5]);
% wave=18 Jun 2015 Election
xline(18.5, '--', 'color',[0.5,0.5,0.5]);
h = text([19], [min(all_nominate.score)], {'Announcement of Trump''s candidacy'});
set(h,'Rotation',90, 'Color', [0.5,0.5,0.5]);
xlim([0.5,41]);
XTICK = [1,7,13,19,25,31,37];
XTICKLABELS = ["Jan 2014", "Jul 2014", "Jan 2015", "Jul 2015",...
    "Jan 2016", "Jul 2016", "Jan 2017"];
ylim([min(all_nominate.score)-0.1,max(all_nominate.score)+0.1])
YTICK = [(min(all_nominate.score)+max(all_nominate.score))/2];
YTICKLABELS = ["Economic Confidene"];

set(gca, 'xtick', XTICK, ...
         'xticklabels', XTICKLABELS,...
         'XTickLabelRotation',0,...
         'ytick', YTICK, ...
         'yticklabels', YTICKLABELS,...
         'YTickLabelRotation',90,...
         'box', 'off', ...
         'tickdir', 'out', ...
    'FontSize',12);
h=gca; h.XAxis.TickLength = [0 0]; h.YAxis.TickLength = [0 0]; 

axp = get(gca,'Position');

% determine startpoint and endpoint for the arrows 
xs=axp(1);
xe=axp(1)+axp(3)+0.01;
ys=axp(2);
ye=axp(2)+axp(4)+0.01;

% make the arrows
annotation('arrow', [xs xe],[ys ys]);
annotation('arrow', [xs xs],[ys ye]);

filename = "./figures/TAPS2014/gpirt_TAPS2014" + ".pdf";
set(fig, 'PaperPosition', [-1 -0.05 10 4]); 
set(fig, 'PaperSize', [8.2 3.8]);
print(fig, filename, '-dpdf','-r300');
close;