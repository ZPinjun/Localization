 clc;clear;close all


%%
figure(1)
centerX = 300;
centerY = 300;
width =590;
height = 420;
set(gcf,'position',[centerX, centerY,width, height])
fs = 8;
txt1 = 'BS 1';
txt2 = 'BS 2';
FS = 9;
d = -0.32;

load('Data_CCRB_vs_pos.mat')
subplot(2,2,1)
imagesc(Pu_range_x, Pu_range_y, data_planar_peb); hold on
[C,h] = contour(Pu_range_x, Pu_range_y, data_planar_peb, [0.05,0.1,0.2,0.3], 'ShowText','on', 'LineColor', 'w');
clabel(C,h,'FontSize',fs,'Color','k','LabelSpacing',100);
plot([-10,-10], [10.1,8], 'r', 'Linewidth', 2);
plot([10-sqrt(2),10], [10,10-sqrt(2)], 'r', 'Linewidth', 2);
caxis([min(min(data_planar_peb)) max(max(data_planar_peb))]);
colorbar
xlabel('y [m]');
ylabel('x [m]');
title('(a) PEB of planar array [m]','FontWeight','Normal','FontName','Times New Roman','FontSize',12,'Units', 'normalized', 'Position', [0.5, d, 0]);
t = text(9,8,txt1,'Color','r','FontSize',FS,'HorizontalAlignment','left');
set(t,'Rotation',90);
text(-9,9,txt2,'Color','r','FontSize',FS,'HorizontalAlignment','left');

subplot(2,2,2)
imagesc(Pu_range_x, Pu_range_y, data_cuboidal_peb); hold on
[C,h] = contour(Pu_range_x, Pu_range_y, data_cuboidal_peb, [0.05,0.1,0.2,0.3], 'ShowText','on', 'LineColor', 'w');
clabel(C,h,'FontSize',fs,'Color','k');
plot([-10,-10], [10.1,8], 'r', 'Linewidth', 2);
plot([10-sqrt(2),10], [10,10-sqrt(2)], 'r', 'Linewidth', 2);
caxis([min(min(data_planar_peb)) max(max(data_planar_peb))]);
colorbar
xlabel('y [m]');
ylabel('x [m]');
title('(b) PEB of cuboidal array [m]','FontWeight','Normal','FontName','Times New Roman','FontSize',12,'Units', 'normalized', 'Position', [0.5, d, 0]);
t = text(9,8,txt1,'Color','r','FontSize',FS,'HorizontalAlignment','left');
set(t,'Rotation',90);
text(-9,9,txt2,'Color','r','FontSize',FS,'HorizontalAlignment','left');

subplot(2,2,3)
imagesc(Pu_range_x, Pu_range_y, data_planar_oeb); hold on
[C,h] = contour(Pu_range_x, Pu_range_y, data_planar_oeb, [1.4,2,3], 'ShowText','on', 'LineColor', 'w');
clabel(C,h,'FontSize',fs,'Color','k','LabelSpacing',80);
plot([-10,-10], [10.1,8], 'r', 'Linewidth', 2);
plot([10-sqrt(2),10], [10,10-sqrt(2)], 'r', 'Linewidth', 2);
caxis([min(min(data_planar_oeb)) 3.5]);
colorbar
xlabel('y [m]');
ylabel('x [m]');
title('(c) OEB of planar array [deg]','FontWeight','Normal','FontName','Times New Roman','FontSize',12,'Units', 'normalized', 'Position', [0.5, d, 0]);
t = text(9,8,txt1,'Color','r','FontSize',FS,'HorizontalAlignment','left');
set(t,'Rotation',90);
text(-9,9,txt2,'Color','r','FontSize',FS,'HorizontalAlignment','left');

subplot(2,2,4)
imagesc(Pu_range_x, Pu_range_y, data_cuboidal_oeb); hold on
[C,h] = contour(Pu_range_x, Pu_range_y, data_cuboidal_oeb, [1.4,2,3], 'ShowText','on', 'LineColor', 'w');
clabel(C,h,'FontSize',fs,'Color','k');
plot([-10,-10], [10.1,8], 'r', 'Linewidth', 2);
plot([10-sqrt(2),10], [10,10-sqrt(2)], 'r', 'Linewidth', 2);
caxis([min(min(data_planar_oeb)) 3.5]);
colorbar
xlabel('y [m]');
ylabel('x [m]');
title('(d) OEB of cuboidal array [deg]','FontWeight','Normal','FontName','Times New Roman','FontSize',12,'Units', 'normalized', 'Position', [0.5, d, 0]);
t = text(9,8,txt1,'Color','r','FontSize',FS,'HorizontalAlignment','left');
set(t,'Rotation',90);
text(-9,9,txt2,'Color','r','FontSize',FS,'HorizontalAlignment','left');





%%
figure(2)
centerX = 300;
centerY = 300;
width =590;
height = 420;
set(gcf,'position',[centerX, centerY,width, height])
fs = 8;
txt1 = 'BS 1';
txt2 = 'BS 2';
FS = 9;
d = -0.32;

load('Data_CCRB_vs_ori.mat')
subplot(2,2,1)
h = imagesc(Ru_range_y, Ru_range_z, data_planar_peb); hold on
set(h,'alphadata',~isnan(data_planar_peb));
% contour(Ru_range_y, Ru_range_z, data_planar_peb.', 'ShowText','on', 'LineColor', 'w');
colorbar
caxis([min(min(data_cuboidal_peb)) 0.2]);
xlabel('\gamma [deg]');
ylabel('\beta [deg]');
xticks([0 90 180 270 360]);
yticks([0 90 180 270 360]);
ax.TickLabelInterpreter = 'latex';
title('(a) PEB of planar array [m]','FontWeight','Normal','FontName','Times New Roman','FontSize',12,'Units', 'normalized', 'Position', [0.5, d, 0]);

subplot(2,2,2)
h = imagesc(Ru_range_y, Ru_range_z, data_cuboidal_peb); hold on
set(h,'alphadata',~isnan(data_cuboidal_peb));
% contour(Ru_range_y, Ru_range_z, data_cuboidal_peb.', 'ShowText','on', 'LineColor', 'w');
colorbar
caxis([min(min(data_cuboidal_peb)) 0.2]);
xlabel('\gamma [deg]');
ylabel('\beta [deg]');
xticks([0 90 180 270 360]);
yticks([0 90 180 270 360]);
ax.TickLabelInterpreter = 'latex';
title('(b) PEB of cuboidal array [m]','FontWeight','Normal','FontName','Times New Roman','FontSize',12,'Units', 'normalized', 'Position', [0.5, d, 0]);

subplot(2,2,3)
h = imagesc(Ru_range_y, Ru_range_z, data_planar_oeb); hold on
set(h,'alphadata',~isnan(data_planar_oeb));
% contour(Ru_range_y, Ru_range_z, data_planar_oeb.', 'ShowText','on', 'LineColor', 'w');
colorbar
caxis([min(min(data_cuboidal_oeb)) 4]);
xlabel('\gamma [deg]');
ylabel('\beta [deg]');
xticks([0 90 180 270 360]);
yticks([0 90 180 270 360]);
ax.TickLabelInterpreter = 'latex';
title('(c) OEB of planar array [deg]','FontWeight','Normal','FontName','Times New Roman','FontSize',12,'Units', 'normalized', 'Position', [0.5, d, 0]);

subplot(2,2,4)
h = imagesc(Ru_range_y, Ru_range_z, data_cuboidal_oeb); hold on
set(h,'alphadata',~isnan(data_cuboidal_oeb));
% contour(Ru_range_y, Ru_range_z, data_cuboidal_oeb.', 'ShowText','on', 'LineColor', 'w');
colorbar
caxis([min(min(data_cuboidal_oeb)) 4]);
xlabel('\gamma [deg]');
ylabel('\beta [deg]');
xticks([0 90 180 270 360]);
yticks([0 90 180 270 360]);
ax.TickLabelInterpreter = 'latex';
title('(d) OEB of cuboidal array [deg]','FontWeight','Normal','FontName','Times New Roman','FontSize',12,'Units', 'normalized', 'Position', [0.5, d, 0]);






%%
figure(3)
centerX = 300;
centerY = 300;
width =660;
height = 480;
set(gcf,'position',[centerX, centerY,width, height])

load('Data_CCRB_vs_Kr.mat')
subplot(2,1,1)
plot(Kr_test, data_planar_peb(1,:),'bo-'); hold on
plot(Kr_test, data_cuboidal_peb(1,:),'rs-'); hold off
set(gca(),'XScale','log');
set(gca(),'YScale','log');
xlabel('Kr');
ylabel('PEB');
legend('Planar','Cuboidal')
grid on

subplot(2,1,2)
plot(Kr_test, data_planar_oeb(1,:),'bo-'); hold on
plot(Kr_test, data_cuboidal_oeb(1,:),'rs-'); hold off
set(gca(),'XScale','log');
set(gca(),'YScale','log');
xlabel('Kr');
ylabel('OEB');
legend('Planar','Cuboidal')
grid on


%%
load('Data_Coverage_CCRB_vBS.mat')
figure(4)
centerX = 300;
centerY = 300;
width =660;
height = 480;
set(gcf,'position',[centerX, centerY,width, height])

subplot(2,1,1)
[f,x] =ecdf(data_planar_2BS(:,1),'Bounds','on'); 
x(isinf(x)) = []; f = f(1:length(x));
x = x(1:10:end); f = f(1:10:end);
plot(x,1-f,'r--','linewidth',1.5); hold on
[f,x] =ecdf(data_planar_3BS(:,1),'Bounds','on');
x(isinf(x)) = []; f = f(1:length(x));
x = x(1:10:end); f = f(1:10:end);
plot(x,1-f,'c--','linewidth',1.5); hold on
[f,x] =ecdf(data_planar_4BS(:,1),'Bounds','on'); 
x(isinf(x)) = []; f = f(1:length(x));
x = x(1:10:end); f = f(1:10:end);
plot(x,1-f,'b--','linewidth',1.5); hold on
[f,x] = ecdf(data_cuboidal_2BS(:,1),'Bounds','on'); 
x = x(1:100:end); f = f(1:100:end);
plot(x,1-f,'r-','linewidth',1.5); hold on
[f,x] = ecdf(data_cuboidal_3BS(:,1),'Bounds','on'); 
x = x(1:100:end); f = f(1:100:end);
plot(x,1-f,'c-','linewidth',1.5); hold on
[f,x] = ecdf(data_cuboidal_4BS(:,1),'Bounds','on'); 
x = x(1:100:end); f = f(1:100:end);
plot(x,1-f,'b-','linewidth',1.5); hold on
set(gca(),'XScale','log');
set(gca(),'YScale','log');
xlabel('PEB threshold [m]');
ylabel('CCDF');
title('Empirical CCDF of PEB');
xlim([1e-2, 3]);
ylim([1e-2, 1]);
grid on
legend('2 BSs, (planar)','3 BSs, (planar)','4 BSs, (planar)','2 BSs, (cuboidal)','3 BSs, (cuboidal)','4 BSs, (cuboidal)');
yticks([1e-4,1e-3,1e-2,1e-1,1]);

subplot(2,1,2)
[f, x] = ecdf(data_planar_2BS(:,2),'Bounds','on'); 
x = x(1:100:end); f = f(1:100:end);
plot(x,1-f,'r--','linewidth',1.5); hold on
[f, x] = ecdf(data_planar_3BS(:,2),'Bounds','on'); 
x = x(1:100:end); f = f(1:100:end);
plot(x,1-f,'c--','linewidth',1.5); hold on
[f, x] = ecdf(data_planar_4BS(:,2),'Bounds','on'); 
x = x(1:100:end); f = f(1:100:end);
plot(x,1-f,'b--','linewidth',1.5); hold on
[f, x] = ecdf(data_cuboidal_2BS(:,2),'Bounds','on'); 
x = x(1:100:end); f = f(1:100:end);
plot(x,1-f,'r-','linewidth',1.5); hold on
[f, x] = ecdf(data_cuboidal_3BS(:,2),'Bounds','on'); 
x = x(1:100:end); f = f(1:100:end);
plot(x,1-f,'c-','linewidth',1.5); hold on
[f, x] = ecdf(data_cuboidal_4BS(:,2),'Bounds','on'); 
x = x(1:100:end); f = f(1:100:end);
plot(x,1-f,'b-','linewidth',1.5); hold on
set(gca(),'XScale','log');
set(gca(),'YScale','log');
xlabel('OEB threshold [deg]');
ylabel('CCDF');
title('Empirical CCDF of OEB');
xlim([1e-1, 30]);
ylim([1e-2, 1]);
grid on
legend('2 BSs, (planar)','3 BSs, (planar)','4 BSs, (planar)','2 BSs, (cuboidal)','3 BSs, (cuboidal)','4 BSs, (cuboidal)');
xticks([1 5 10 30])
xticklabels({'1','5','10','30'})




%
load('Data_Coverage_CCRB_vDirec.mat')
figure(5)
centerX = 300;
centerY = 300;
width =660;
height = 480;
set(gcf,'position',[centerX, centerY,width, height])
subplot(2,1,1)
[f,x] =ecdf(data_planar_v1(:,1),'Bounds','on'); 
x(isinf(x)) = []; f = f(1:length(x));
x = x(1:50:end); f = f(1:50:end);
plot(x,1-f,'r--','linewidth',1.5); hold on
[f,x] =ecdf(data_planar_v2(:,1),'Bounds','on'); 
x(isinf(x)) = []; f = f(1:length(x));
x = x(1:50:end); f = f(1:50:end);
plot(x,1-f,'c--','linewidth',1.5); hold on
[f,x] =ecdf(data_planar_v3(:,1),'Bounds','on'); 
x(isinf(x)) = []; f = f(1:length(x));
x = x(1:50:end); f = f(1:50:end);
plot(x,1-f,'b--','linewidth',1.5); hold on
[f,x] = ecdf(data_cuboidal_v1(:,1),'Bounds','on'); 
x(isinf(x)) = []; f = f(1:length(x));
x = x(1:5:end); f = f(1:5:end);
plot(x,1-f,'r-','linewidth',1.5); hold on
[f,x] = ecdf(data_cuboidal_v2(:,1),'Bounds','on'); 
x(isinf(x)) = []; f = f(1:length(x));
x = x(1:5:end); f = f(1:5:end);
plot(x,1-f,'c-','linewidth',1.5); hold on
[f,x] = ecdf(data_cuboidal_v3(:,1),'Bounds','on'); 
x = x(1:100:end); f = f(1:100:end);
plot(x,1-f,'b-','linewidth',1.5); hold on
set(gca(),'XScale','log');
set(gca(),'YScale','log');
xlabel('PEB threshold [m]');
ylabel('CCDF');
title('Empirical CCDF of PEB');
ylim([1e-2, 1]);
grid on
legend('2D, V=120','2D, V=150','2D, V=180','3D, V=120','3D, V=150','3D, V=180');
yticks([1e-4,1e-3,1e-2,1e-1,1]);
subplot(2,1,2)
[f, x] = ecdf(data_planar_v1(:,2),'Bounds','on'); 
x(isinf(x)) = []; f = f(1:length(x));
x = x(1:50:end); f = f(1:50:end);
plot(x,1-f,'r--','linewidth',1.5); hold on
[f, x] = ecdf(data_planar_v2(:,2),'Bounds','on'); 
x(isinf(x)) = []; f = f(1:length(x));
x = x(1:50:end); f = f(1:50:end);
plot(x,1-f,'c--','linewidth',1.5); hold on
[f, x] = ecdf(data_planar_v3(:,2),'Bounds','on'); 
x(isinf(x)) = []; f = f(1:length(x));
x = x(1:50:end); f = f(1:50:end);
plot(x,1-f,'b--','linewidth',1.5); hold on
[f, x] = ecdf(data_cuboidal_v1(:,2),'Bounds','on'); 
x(isinf(x)) = []; f = f(1:length(x));
x = x(1:50:end); f = f(1:50:end);
plot(x,1-f,'r-','linewidth',1.5); hold on
[f, x] = ecdf(data_cuboidal_v2(:,2),'Bounds','on'); 
x(isinf(x)) = []; f = f(1:length(x));
x = x(1:50:end); f = f(1:50:end);
plot(x,1-f,'c-','linewidth',1.5); hold on
[f, x] = ecdf(data_cuboidal_v3(:,2),'Bounds','on'); 
x = x(1:100:end); f = f(1:100:end);
plot(x,1-f,'b-','linewidth',1.5); hold on
set(gca(),'XScale','log');
set(gca(),'YScale','log');
xlabel('OEB threshold [deg]');
ylabel('CCDF');
title('Empirical CCDF of OEB');
ylim([1e-2, 1]);
grid on
legend('2D, V=120','2D, V=150','2D, V=180','3D, V=120','3D, V=150','3D, V=180');



