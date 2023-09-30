clc; clear; close all

%% Add Paths
path(pathdef); addpath(pwd);
cd ..;
cd Channel; addpath(genpath(pwd)); cd ..;
cd System; addpath(genpath(pwd)); cd ..; 
cd Localization; addpath(genpath(pwd)); cd ..; 
cd !Simulations;


%% Simulation setup
repeatNum = 10000;
TesThereshold = [17,20,23];
TesThereshold = 10.^(TesThereshold/10);
data_planar_t1 = zeros(repeatNum,1);
data_cuboidal_t1 = zeros(repeatNum,1);
data_planar_t2 = zeros(repeatNum,1);
data_cuboidal_t2 = zeros(repeatNum,1);
data_planar_t3 = zeros(repeatNum,1);
data_cuboidal_t3 = zeros(repeatNum,1);
randpos = rand(3,repeatNum);
randori = rand(3,repeatNum);
pcp.func = 'communication';


parfor j = 1:repeatNum
    x = 20*randpos(1,j) - 10;
    y = 20*randpos(2,j) - 10;
    z = 5*randpos(3,j);
    pu = [x, y, z].';
    ou_euler = 360*randori(:,j);
    disp(['j = ', num2str(j), ',   --------- ', num2str(100*(j/repeatNum)),'% start...']);

    %% planar array
    % default setup
    sp = default_system_setup();
    cp = default_channel_setup();
    % update setup
    sp.TypeSA = 'planar';
    sp.Pu = pu;
    sp.Ou_Euler = ou_euler;
    cp.K = 2;   % # of subcarriers
    cp.G = 1;   % # of transmissions
    sp = update_system_setup(sp);
    cp = update_channel_setup(cp);
    % generate channel
    cp = gen_channel(sp,cp);
    % generate precoders & combiners
    cp = gen_precoder_combiner(sp,cp,pcp);
    % get outage probability of m-th BS at subcarrier fk
    k = floor(cp.K/2);
    method = "analytical";
    % t1
    threshold = TesThereshold(1);
    data_planar_t1(j) = prod(get_outage_probability(sp,cp,threshold,k,method));
    % t2
    threshold = TesThereshold(2);
    data_planar_t2(j) = prod(get_outage_probability(sp,cp,threshold,k,method));
    % t3
    threshold = TesThereshold(3);
    data_planar_t3(j) = prod(get_outage_probability(sp,cp,threshold,k,method));


    %% cuboidal array
    % default setup
    sp = default_system_setup();
    cp = default_channel_setup();
    % update setup
    sp.TypeSA = 'cuboidal';
    sp.Pu = pu;
    sp.Ou_Euler = ou_euler;
    cp.K = 2;   % # of subcarriers
    cp.G = 1;   % # of transmissions
    sp = update_system_setup(sp);
    cp = update_channel_setup(cp);
    % generate channel
    cp = gen_channel(sp,cp);
    % generate precoders & combiners
    cp = gen_precoder_combiner(sp,cp,pcp);
    % get outage probability of m-th BS at subcarrier fk
    k = floor(cp.K/2);
    method = "analytical";
    % t1
    threshold = TesThereshold(1);
    data_cuboidal_t1(j) = prod(get_outage_probability(sp,cp,threshold,k,method));
    % t2
    threshold = TesThereshold(2);
    data_cuboidal_t2(j) = prod(get_outage_probability(sp,cp,threshold,k,method));
    % t3
    threshold = TesThereshold(3);
    data_cuboidal_t3(j) = prod(get_outage_probability(sp,cp,threshold,k,method));

end

figure(1)
[f,x] =ecdf(data_planar_t1(:,1),'Bounds','on'); 
plot(x,f,'b-'); hold on
[f,x] =ecdf(data_planar_t2(:,1),'Bounds','on'); 
plot(x,f,'b-.'); hold on
[f,x] =ecdf(data_planar_t3(:,1),'Bounds','on'); 
plot(x,f,'b--'); hold on
[f,x] = ecdf(data_cuboidal_t1(:,1),'Bounds','on'); 
plot(x,f,'r-'); hold on
[f,x] = ecdf(data_cuboidal_t2(:,1),'Bounds','on'); 
plot(x,f,'r-.'); hold on
[f,x] = ecdf(data_cuboidal_t3(:,1),'Bounds','on'); 
plot(x,f,'r--'); hold off
%set(gca(),'XScale','log');
% set(gca(),'YScale','log');
xlabel('Threshold \xi [percent]');
ylabel('Non-outage coverage');
title('Non-outage coverage vs. threshold');
% xlim([1e-3, 1e4]);
grid on
lgd = legend('Cuboidal, gamma=17 dB','Planar, gamma=17 dB','Cuboidal, gamma=20 dB','Planar, gamma=20 dB','Cuboidal, gamma=23 dB','Planar, gamma=23 dB');
lgd.NumColumns = 3;

% save('Data_Coverage_Pout.mat','data_planar_t1','data_planar_t2','data_planar_t3','data_cuboidal_t1','data_cuboidal_t2','data_cuboidal_t3')
