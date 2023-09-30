clc; clear; close all
rng('default');

%% Add Paths
path(pathdef); addpath(pwd);
cd ..;
cd Channel; addpath(genpath(pwd)); cd ..;
cd System; addpath(genpath(pwd)); cd ..; 
cd Localization; addpath(genpath(pwd)); cd ..; 
cd !Simulations;


%% simulation setup
Ptest = -25:5:45;
repeatNum = 300;
data_Pu_LS = zeros(length(Ptest),1);
data_Ru_LS = zeros(length(Ptest),1);
data_Pu_ML = zeros(length(Ptest),1);
data_Ru_ML = zeros(length(Ptest),1);
data_CCRB = zeros(length(Ptest),2);
pcp.func = 'localization';


for i = 1:length(Ptest)
    P = db2pow(Ptest(i));
    disp(['P = ', num2str(Ptest(i)),'[dBm],    ', num2str(i/length(Ptest)*100), '%  start...']);

    % default setup
    sp = default_system_setup();
    cp = default_channel_setup();   

    % update setup
    sp.TypeSA = 'cuboidal';
    cp.P = P;
    cp.Kr = 1e20;
    cp.K = 2;   % # of subcarriers
    cp.G = 640;   % # of transmissions
    % ---------------------------
    % generate UE
%     sp.Pu = [1 3 2].';
%     sp.Ou_Euler = [30 40 50].';
    % ---------------------------
    sp = update_system_setup(sp);
    cp = update_channel_setup(cp);

    % generate channel
    cp = gen_channel(sp,cp);

    % generate precoders & combiners
    cp = gen_precoder_combiner(sp,cp,pcp);

    % get CRLB
    CCRB = get_CCRB(cp,sp);
    data_CCRB(i,1) = CCRB(1);
    data_CCRB(i,2) = CCRB(2);
    
    % get RMSE
    Pu = sp.Pu;
    Ru = sp.Ru;
    [RMSE_Pu_LS,RMSE_Ru_LS,RMSE_Pu_ML,RMSE_Ru_ML] = get_RMSE(repeatNum,cp,sp,Pu,Ru,data_CCRB(i,:));

    %% save data
    data_Pu_LS(i) = RMSE_Pu_LS;
    data_Ru_LS(i) = RMSE_Ru_LS;
    data_Pu_ML(i) = RMSE_Pu_ML;
    data_Ru_ML(i) = RMSE_Ru_ML;
end

figure(1)
semilogy(Ptest,data_CCRB(:,1),'ro--'); hold on
semilogy(Ptest,data_CCRB(:,2),'go--'); hold on
semilogy(Ptest,data_Pu_LS,'r^-.'); hold on
semilogy(Ptest,data_Ru_LS,'g^-.'); hold on
semilogy(Ptest,data_Pu_ML,'rv-'); hold on
semilogy(Ptest,data_Ru_ML,'gv-'); hold off


% save('Data_RMSE_vs_P.mat','Ptest','data_CCRB','data_Pu_LS','data_Pu_ML','data_Ru_LS','data_Ru_ML');


