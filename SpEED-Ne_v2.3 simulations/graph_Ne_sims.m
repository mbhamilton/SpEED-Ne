

% code to analyze Ne simulation runs

% load saved .mat file into workspace
load('27Sept_medusa_test.mat');

% Ne(:,1) - r_squared_c AFW correction_factor_har_mean_S
% Ne(:,2) - r_squared_c AFT correction_factor_har_mean_S
% Ne(:,3) - r_squared_delta AFW r_delta_correction_factor_har_mean_S
% Ne(:,4) - r_squared_delta AFT r_delta_correction_factor_har_mean_S
% Ne(:,5) - r_squared_c AFT Waples regression
% Ne(:,6) - r_squared_c AFW median_r_sqr_c_AFW_perm
% Ne(:,7) - r_squared_c AFT median_r_sqr_c_AFT_perm


print -deps2 '/Users/matthewhamilton/Dropbox/R project/manuscript/time forward Ne simulation results/28_Sept_medusa_test_50.eps'



%boxplot(Ne(:,1,1),Ne(:,2,1),Ne(:,3,1),Ne(:,4,1),Ne(:,5,1),Ne(:,6,1),Ne(:,7,1)); 

figure(1);
boxplot(Ne,'Colors','k','Symbol','k+','DataLim',[0,250],'ExtremeMode','compress','Labels',{'r^2_C AFT','r^2_C AFW','r^2_\Delta AFW','r^2_\Delta AFT','Waples','r^2_C perm AFT','r^2_C perm AFW'})
hold on;
ax = gca;
ax.YGrid = 'on';
ax.YTick = -100:50:400;
hold off;

h = findobj(gca, 'type', 'text');
set(h, 'Interpreter', 'tex');
set(h,'FontSize',14);
hold off;


print -deps2 myplot.eps

% -deps .eps black and white 
% -depsc .eps color
% -deps2 .eps level 2 black and white
% -depsc2 .eps level 2 color 
% -dpdf .pdf color file format


'Interpreter','tex',

'DataLim',[-250,250]

'ExtremeMode', 'compress' % 'clip' (default)

boxplot(Ne)
set(axes,'TickLabelInterpreter','none','XTick',[1 2 3 4 5 6 7],...
    'XTickLabel',...
    {'r^2_{C} AFT','r^2_{C} AFW','r^2_{\Delta} AFW','r^2_{\Delta} AFT','Waples','r^2_{C} perm AFT','r^2_{C} perm AFW'},...
    'YGrid','on','YTick',...
    [-500 -450 -400 -350 -300 -250 -200 -150 -100 -50 0 50 100 150 200 250 300 350 400 450 500]);

boxplot(Ne)
set(axes,'TickLabelInterpreter','tex','XTick',[1 2 3 4 5 6 7],...
    'XTickLabel',...
    {'r^2_{comp} AFT','r^2_{comp} AFW','r^2_{delta} AFW','r^2_{delta} AFT','Waples','r^2_{comp} perm AFT','r^2_{comp} perm AFW'},...
    'YGrid','on','YTick',...
    [-500 -450 -400 -350 -300 -250 -200 -150 -100 -50 0 50 100 150 200 250 300 350 400 450 500]);

boxplot(Ne)
set(axes,'TickLabelInterpreter','tex','XTick',[1 2 3 4 5 6 7],...
    'XTickLabel',...
    {'r^2_{comp} AFT','r^2_{comp} AFW','r^2_{delta} AFW','r^2_{delta} AFT','Waples','r^2_{comp} perm AFT','r^2_{comp} perm AFW'});


