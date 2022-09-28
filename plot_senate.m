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
tiledlayout(4,1,'Padding', 'none', 'TileSpacing', 'compact');
polarization = [];
nexttile;
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
set(gca,'box','off','ytick',[0.4:0.1:1.0],...
    'XTickLabel',(session_ids-min(session_ids))*2+1971,'fontsize',8);
legend([p1],{'Polarization ratio'},...
    'Location','northwest','FontSize',12, 'NumColumns' ,1);
legend boxoff;

nexttile([3 1]);
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
ylim([-3.0,4.0]);
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
annotation('arrow', [0.0205 0.0205],[0.65 0.75], 'LineWidth', 0.01);

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
legend([p1,p2],{'Democrats','Republicans'},...
    'Location','northeast','FontSize',12, 'NumColumns' ,2);
legend boxoff;

set(fig, 'PaperPosition', [0 0 10 4]); 
set(fig, 'PaperSize', [10 4]); 

filename = "./results/abortion_gpirt.pdf";
print(fig, filename, '-dpdf','-r300', '-fillpage');
close;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
session_ids = 92:2:110;
fig = figure(1);
tiledlayout(4,1,'Padding', 'none', 'TileSpacing', 'compact');
nexttile;
polarization = [];
for i=1:numel(session_ids)
    session_id = session_ids(i);
    x = all_nominate.gpirt(all_nominate.session==session_id);
    party = all_nominate.party(all_nominate.session==session_id);
    polarization = [polarization, sum(strcmp(party,'Republicans')==(x>=0))/numel(x)];
end
p1 = plot((1:1:numel(session_ids))-1,polarization,'-','Color',...
    [0.5,0.5,0.5], 'LineWidth',2);

ylim([0.4,1.0]);
YTICKS = (1:1:numel(session_ids))-1;
xticks(YTICKS);
xlim([-0.6,(numel(session_ids))+0.2]);
xtickformat('%dth');
% xticklabels(session_ids);
a = get(gca,'XTickLabel');
set(gca,'box','off','ytick',[0.4:0.1:1.0],...
    'XTickLabel',(session_ids-min(session_ids))*2+1971,'fontsize',8);
legend([p1],{'Polarization ratio'},...
    'Location','northwest','FontSize',12, 'NumColumns' ,1);
legend boxoff;

nexttile([3 1]);
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
annotation('arrow', [0.0205 0.0205],[0.65 0.75], 'LineWidth', 0.01);

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
legend([p1,p2],{'Democrats','Republicans'},...
    'Location','northeast','FontSize',12, 'NumColumns' ,2);
legend boxoff;

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

colors = ["k", "b","r","magenta"];
score_types = flip(unique(all_nominate.type));
YLABELS = ["GD-GPIRT", "DO-IRT"];
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
        h{i} = plot(x,y,strcat(shapes(k),colors(i)),'MarkerSize',8, 'LineWidth', 4); hold on;
    end
    
    xlabel('year','FontSize', 18);
    % ylim([-1.0*(3-k),1.0*(3-k)]);
    ylim([-2.5,1.5]);
    ylabel(char(score_types(k)),'FontSize', 18);
    ylabel(YLABELS(k),'FontSize', 18);
    legend(country_names, 'Location','northwest','FontSize',12, 'NumColumns', 4);
    legend boxoff;
end
    
set(fig, 'PaperPosition', [0 0 10 4]); 
set(fig, 'PaperSize', [10 4]); 

filename = "./results/CIRI_dynamic.pdf";
print(fig, filename, '-dpdf','-r300', '-fillpage');
close;
