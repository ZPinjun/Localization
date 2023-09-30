clc; clear; close all

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

load('Data_capacity_vs_pos.mat')
data_planar_capacity = data_planar_capacity/1e9;
data_cuboidal_capacity = data_cuboidal_capacity/1e9;

subplot(2,2,1)
imagesc(Pu_range_x, Pu_range_y, data_planar_capacity); hold on
% [C,h] = contour(Pu_range_x, Pu_range_y, data_planar_capacity, [16,17,18,19,20], 'ShowText','on', 'LineColor', 'w');
title('(a) Capacity (planar) over different position','FontWeight','Normal','FontName','Times New Roman','FontSize',12,'Units', 'normalized', 'Position', [0.5, d, 0]);
colorbar
caxis([min(min(data_cuboidal_capacity)) max(max(data_planar_capacity))]);
% clabel(C,h,'FontSize',fs,'Color','k');
scatter(0, 0,'rx','Linewidth', 2);
xlabel('y [m]');
ylabel('x [m]');
plot([-10,-10], [10.1,8], 'r', 'Linewidth', 2);
plot([10-sqrt(2),10], [10,10-sqrt(2)], 'r', 'Linewidth', 2);
t = text(9,8,txt1,'Color','r','FontSize',FS,'HorizontalAlignment','left');
set(t,'Rotation',90);
text(-9,9,txt2,'Color','r','FontSize',FS,'HorizontalAlignment','left');

subplot(2,2,2)
imagesc(Pu_range_x, Pu_range_y, data_cuboidal_capacity); hold on
% [C,h] = contour(Pu_range_x, Pu_range_y, data_cuboidal_capacity, [16,17,18,19,20], 'ShowText','on', 'LineColor', 'w');
title('(b) Capacity (cuboidal) over different position','FontWeight','Normal','FontName','Times New Roman','FontSize',12,'Units', 'normalized', 'Position', [0.5, d, 0]);
colorbar
caxis([min(min(data_cuboidal_capacity)) max(max(data_planar_capacity))]);
% clabel(C,h,'FontSize',fs,'Color','k');
scatter(0, 0,'rx','Linewidth', 2);
xlabel('y [m]');
ylabel('x [m]');
plot([-10,-10], [10.1,8], 'r', 'Linewidth', 2);
plot([10-sqrt(2),10], [10,10-sqrt(2)], 'r', 'Linewidth', 2);
t = text(9,8,txt1,'Color','r','FontSize',FS,'HorizontalAlignment','left');
set(t,'Rotation',90);
text(-9,9,txt2,'Color','r','FontSize',FS,'HorizontalAlignment','left');

load('Data_capacity_vs_ori.mat')
data_planar_capacity = data_planar_capacity/1e9;
data_cuboidal_capacity = data_cuboidal_capacity/1e9;
subplot(2,2,3)
imagesc(Ru_range_y, Ru_range_z, data_planar_capacity); hold on
title('(c) Capacity (planar) over different orientation','FontWeight','Normal','FontName','Times New Roman','FontSize',12,'Units', 'normalized', 'Position', [0.5, d, 0]);
colorbar
caxis([min(min(data_planar_capacity)) max(max(data_planar_capacity))]);
xlabel('\gamma [deg]');
ylabel('\beta [deg]');

subplot(2,2,4)
imagesc(Ru_range_y, Ru_range_z, data_cuboidal_capacity); hold on
title('(d) Capacity (cuboidal) over different orientation','FontWeight','Normal','FontName','Times New Roman','FontSize',12,'Units', 'normalized', 'Position', [0.5, d, 0]);
colorbar
caxis([min(min(data_planar_capacity)) max(max(data_planar_capacity))]);
xlabel('\gamma [deg]');
ylabel('\beta [deg]');






figure(2)
centerX = 300;
centerY = 300;
width =590;
height = 420;
set(gcf,'position',[centerX, centerY,width, height])

load('Data_capacity_vs_Kr.mat')
data_capacity_planar= data_capacity_planar/1e9;
data_capacity_cuboidal = data_capacity_cuboidal/1e9;
semilogx(RangeKr,data_capacity_planar(3,:),'b^--'); hold on
semilogx(RangeKr,data_capacity_planar(2,:),'bs--'); hold on
semilogx(RangeKr,data_capacity_planar(1,:),'bv--'); hold on
semilogx(RangeKr,data_capacity_cuboidal(3,:),'r^-'); hold on
semilogx(RangeKr,data_capacity_cuboidal(2,:),'rs-'); hold on
semilogx(RangeKr,data_capacity_cuboidal(1,:),'rv-'); hold off
legend('Planar,P=15mW','Planar,P=10mW','Planar,P=5mW','Cuboidal,P=15mW','Cuboidal,P=10mW','Cuboidal,P=5mW');
grid on
xlabel('Kr');
ylabel('Ergodic capacity [Gbps]');

figure(3)
centerX = 300;
centerY = 300;
width =590;
height = 420;
set(gcf,'position',[centerX, centerY,width, height])

load('Data_Coverage_capacity.mat')

[f,x] =ecdf(data_planar_v1(:,1),'Bounds','on'); x = x([1:50,51:100:end]); f = f([1:50,51:100:end]);
plot(x/1e9,f,'r--','linewidth',1.5); hold on
[f,x] =ecdf(data_planar_v2(:,1),'Bounds','on'); x = x([1:50,51:100:end]); f = f([1:50,51:100:end]);
plot(x/1e9,f,'c--','linewidth',1.5); hold on
[f,x] =ecdf(data_planar_v3(:,1),'Bounds','on'); x = x([1:50,51:100:end]); f = f([1:50,51:100:end]);
plot(x/1e9,f,'b--','linewidth',1.5); hold on
[f,x] = ecdf(data_cuboidal_v1(:,1),'Bounds','on'); x = x(1:100:end); f = f(1:100:end);
plot(x/1e9,f,'r-','linewidth',1.5); hold on
[f,x] = ecdf(data_cuboidal_v2(:,1),'Bounds','on'); x = x(1:100:end); f = f(1:100:end);
plot(x/1e9,f,'c-','linewidth',1.5); hold on
[f,x] = ecdf(data_cuboidal_v3(:,1),'Bounds','on'); x = x(1:100:end); f = f(1:100:end);
plot(x/1e9,f,'b-','linewidth',1.5); hold off
%set(gca(),'XScale','log');
set(gca(),'YScale','log');
xlabel('Ergodic capacity threshold [Gbps]');
ylabel('1-Ergodic capacity coverage');
title('Ergodic capacity coverage vs. threshold');
% xlim([6, 30]);
ylim([1e-2, 1]);
grid on
legend('Planar, Kr=1','Planar, Kr=4','Planar, Kr=16','Cuboidal, Kr=1','Cuboidal, Kr=4','Cuboidal, Kr=16')


