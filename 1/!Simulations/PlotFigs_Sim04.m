clc; clear; close all

figure(1)
centerX = 300;
centerY = 300;
width =590;
height = 420;
set(gcf,'position',[centerX, centerY,width, height])

load('Data_RMSE_vs_P.mat')
semilogy(Ptest,data_CCRB(:,1),'ro--'); hold on
semilogy(Ptest,data_CCRB(:,2),'go--'); hold on
semilogy(Ptest,data_Pu_LS,'r^-.'); hold on
semilogy(Ptest,data_Ru_LS,'g^-.'); hold on
semilogy(Ptest,data_Pu_ML,'rv-'); hold on
semilogy(Ptest,data_Ru_ML,'gv-'); hold off
grid on
xlabel('Average transmission power');
ylabel('RMSE');