clc; clear; close all

figure(10)
centerX = 300;
centerY = 300;
width =590;
height = 420;
set(gcf,'position',[centerX, centerY,width, height])
txt1 = 'BS 1';
txt2 = 'BS 2';
FS = 9;
d = -0.32;

load('Data_Pout_vs_pos.mat')
f1 = subplot(2,2,1);
imagesc(Pu_range_x, Pu_range_y, data_planar_sumrate/1e9); hold on
title('sum-rate');   
colormap(f1,parula)
colorbar
xlabel('y [m]');
ylabel('x [m]');
title('(a) maximum sum-rate over BSs [Gbps]','Interpreter','latex','FontSize',12,'Units', 'normalized', 'Position', [0.5, d, 0]);
plot([-10,-10], [10.1,8], 'r', 'Linewidth', 2);
plot([10-sqrt(2),10], [10,10-sqrt(2)], 'r', 'Linewidth', 2);
t = text(9,8,txt1,'Color','r','FontSize',FS,'HorizontalAlignment','left');
set(t,'Rotation',90);
text(-9,9,txt2,'Color','r','FontSize',FS,'HorizontalAlignment','left');
xticks([-10 -5 0 5 10]);
yticks([-10 -5 0 5 10]);

f2 = subplot(2,2,2);
imagesc(Pu_range_x, Pu_range_y, data_planar_mtilde); hold on
title('m_tilde');
C = [0.6  0.2   1;
    0 0.6 0];
C(end+1,:) = 1;
colormap(f2, flipud(C));
colorbar('Ticks', [0 1 2])
xlabel('y [m]');
ylabel('x [m]');
caxis([0 2]);
title('(b) Selection of $\tilde{m}$','Interpreter','latex','FontSize',12,'Units', 'normalized', 'Position', [0.5, d, 0]);
plot([-10,-10], [10.1,8], 'r', 'Linewidth', 2);
plot([10-sqrt(2),10], [10,10-sqrt(2)], 'r', 'Linewidth', 2);
t = text(9,8,txt1,'Color','r','FontSize',FS,'HorizontalAlignment','left');
set(t,'Rotation',90);
text(-9,9,txt2,'Color','r','FontSize',FS,'HorizontalAlignment','left');
xticks([-10 -5 0 5 10]);
yticks([-10 -5 0 5 10]);

load('Data_Pout_vs_ori.mat')
f3 = subplot(2,2,3);
imagesc(Ru_range_y, Ru_range_z, data_planar_sumrate/1e9); hold on
title('sum-rate');
colormap(f1,parula)
colorbar
xlabel('\gamma [deg]');
ylabel('\beta [deg]');
title('(c) maximum sum-rate over BSs [Gbps]','Interpreter','latex','FontSize',12,'Units', 'normalized', 'Position', [0.5, d, 0]);
xticks([0 90 180 270 360]); 
yticks([0 90 180 270 360]);

f4 = subplot(2,2,4);
imagesc(Ru_range_y, Ru_range_z, data_planar_mtilde); hold on
title('m_tilde');
C = [0.6  0.2   1;
    0 0.6 0];
C(end+1,:) = 1;
colormap(f4, flipud(C));
colorbar('Ticks', [0 1 2])
xlabel('\gamma [deg]');
ylabel('\beta [deg]');
caxis([0 2]);
title('(d) Selection of $\tilde{m}$','Interpreter','latex','FontSize',12,'Units', 'normalized', 'Position', [0.5, d, 0]);
xticks([0 90 180 270 360]);
yticks([0 90 180 270 360]);




figure(1)
centerX = 300;
centerY = 300;
width =590;
height = 420;
set(gcf,'position',[centerX, centerY,width, height])
txt1 = 'BS 1';
txt2 = 'BS 2';
FS = 9;
d = -0.32;

load('Data_Pout_vs_pos.mat')
subplot(2,2,1)
imagesc(Pu_range_x, Pu_range_y, data_planar_pout); hold on
%contour(Pu_range_x, Pu_range_y, data_planar_pout, 'ShowText','on', 'LineColor', 'w');
title('(a) $P_{\tilde{m},k,\gamma_\mathrm{th}}^{\mathrm{out}}$ (planar) over different position','Interpreter','latex','FontSize',12,'Units', 'normalized', 'Position', [0.5, d, 0]);
colorbar
xlabel('y [m]');
ylabel('x [m]');
xticks([-10 -5 0 5 10]);
yticks([-10 -5 0 5 10]);
plot([-10,-10], [10.1,8], 'r', 'Linewidth', 2);
plot([10-sqrt(2),10], [10,10-sqrt(2)], 'r', 'Linewidth', 2);
scatter(1.8, -5.9,'rx','Linewidth', 2);
txt = 'Uncovered';
text(-5,-8.5,txt,'Color','b')
txt = 'Covered';
text(-2,0,txt,'Color','w')
t = text(9,8,txt1,'Color','r','FontSize',FS,'HorizontalAlignment','left');
set(t,'Rotation',90);
text(-9,9,txt2,'Color','r','FontSize',FS,'HorizontalAlignment','left');

subplot(2,2,2)
imagesc(Pu_range_x, Pu_range_y, data_cuboidal_pout); hold on
%contour(Pu_range_x, Pu_range_y, data_planar_pout, 'ShowText','on', 'LineColor', 'w');
title('(b) $P_{\tilde{m},k,\gamma_\mathrm{th}}^{\mathrm{out}}$ (cuboidal) over different position','Interpreter','latex','FontSize',12,'Units', 'normalized', 'Position', [0.5, d, 0]);
colorbar
xlabel('y [m]');
ylabel('x [m]');
xticks([-10 -5 0 5 10]);
yticks([-10 -5 0 5 10]);
plot([-10,-10], [10.1,8], 'r', 'Linewidth', 2);
plot([10-sqrt(2),10], [10,10-sqrt(2)], 'r', 'Linewidth', 2);
scatter(1.8, -5.9, 'rx','Linewidth', 2);
txt = 'Uncovered';
text(-5,-8.5,txt,'Color','b')
txt = 'Covered';
text(-2,0,txt,'Color','w')
t = text(9,8,txt1,'Color','r','FontSize',FS,'HorizontalAlignment','left');
set(t,'Rotation',90);
text(-9,9,txt2,'Color','r','FontSize',FS,'HorizontalAlignment','left');

load('Data_Pout_vs_ori.mat')
subplot(2,2,3)
imagesc(Ru_range_y, Ru_range_z, data_planar_pout); hold on
title('(c) $P_{\tilde{m},k,\gamma_\mathrm{th}}^{\mathrm{out}}$ (planar) over different orientation','Interpreter','latex','FontSize',12,'Units', 'normalized', 'Position', [0.5, d, 0]);
colorbar
caxis([min(min(data_planar_pout)) 1]);
xlabel('\gamma [deg]');
ylabel('\beta [deg]');
xticks([0 90 180 270 360]);
yticks([0 90 180 270 360]);


subplot(2,2,4)
imagesc(Ru_range_y, Ru_range_z, data_cuboidal_pout); hold on
title('(d) $P_{\tilde{m},k,\gamma_\mathrm{th}}^{\mathrm{out}}$ (cuboidal) over different orientation','Interpreter','latex','FontSize',12,'Units', 'normalized', 'Position', [0.5, d, 0]);
colorbar
caxis([min(min(data_planar_pout)) 1]);
xlabel('\gamma [deg]');
ylabel('\beta [deg]');
xticks([0 90 180 270 360]);
yticks([0 90 180 270 360]);





figure(2)
centerX = 300;
centerY = 300;
width =590;
height = 420;
set(gcf,'position',[centerX, centerY,width, height])

load('Data_Pout_vs_threshold.mat')
a = 1;
semilogy(Range_thereshold,data_Pout_analy(1,:,a),'b-'); hold on
scatter(Range_thereshold,data_Pout_empir(1,:,a),'k'); hold on
semilogy(Range_thereshold,data_Pout_analy(2,:,a),'r-'); hold on
scatter(Range_thereshold,data_Pout_empir(2,:,a),'ks'); hold on
a = 2;
semilogy(Range_thereshold,data_Pout_analy(1,:,a),'b-'); hold on
scatter(Range_thereshold,data_Pout_empir(1,:,a)); hold on
semilogy(Range_thereshold,data_Pout_analy(2,:,a),'r-'); hold on
scatter(Range_thereshold,data_Pout_empir(2,:,a)); hold on
a = 3;
semilogy(Range_thereshold,data_Pout_analy(1,:,a),'b-'); hold on
scatter(Range_thereshold,data_Pout_empir(1,:,a)); hold on
semilogy(Range_thereshold,data_Pout_analy(2,:,a),'r-'); hold on
scatter(Range_thereshold,data_Pout_empir(2,:,a)); hold off
ylim([1e-5, 1]);
xlim([17.5, 19]);
grid on
legend('Planar, analytical','Planar, empirical','Cuboidal, analytical','Cuboidal, empirical');
xlabel('SNR threshold [dB]');
ylabel('Outage probability');




figure(3)
centerX = 300;
centerY = 300;
width =590;
height = 420;
set(gcf,'position',[centerX, centerY,width, height])
load('Data_Coverage_Pout.mat')

[f,x] = ecdf(data_cuboidal_t1(:,1),'Bounds','on');
x = x(1:5:end); f = f(1:5:end);
plot(x,1-f,'r-','linewidth',1.5); hold on
[f,x] =ecdf(data_planar_t1(:,1),'Bounds','on'); 
x = x(1:5:end); f = f(1:5:end);
plot(x,1-f,'r--','linewidth',1.5); hold on
[f,x] = ecdf(data_cuboidal_t2(:,1),'Bounds','on'); 
x = x(1:5:end); f = f(1:5:end);
plot(x,1-f,'c-','linewidth',1.5); hold on
[f,x] =ecdf(data_planar_t2(:,1),'Bounds','on'); 
x = x(1:5:end); f = f(1:5:end);
plot(x,1-f,'c--','linewidth',1.5); hold on
[f,x] = ecdf(data_cuboidal_t3(:,1),'Bounds','on'); 
x = x(1:5:end); f = f(1:5:end);
plot(x,1-f,'b-','linewidth',1.5); hold on
[f,x] =ecdf(data_planar_t3(:,1),'Bounds','on'); 
x = x(1:5:end); f = f(1:5:end);
plot(x,1-f,'b--','linewidth',1.5); hold off

%set(gca(),'XScale','log');
set(gca(),'YScale','log');
ylim([1e-2, 1]);
xlabel('Threshold \xi [percent]');
ylabel('1-Non-outage coverage');
title('1-Non-outage coverage vs. threshold');
grid on
lgd = legend('Cuboidal, gamma=17 dB','Planar, gamma=17 dB','Cuboidal, gamma=20 dB','Planar, gamma=20 dB','Cuboidal, gamma=23 dB','Planar, gamma=23 dB');
lgd.NumColumns = 3;

