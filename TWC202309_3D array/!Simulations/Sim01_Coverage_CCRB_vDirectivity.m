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
v = [120, 150, 180];
data_planar_v1 = zeros(repeatNum,2);
data_cuboidal_v1 = zeros(repeatNum,2);
data_planar_v2 = zeros(repeatNum,2);
data_cuboidal_v2 = zeros(repeatNum,2);
data_planar_v3 = zeros(repeatNum,2);
data_cuboidal_v3 = zeros(repeatNum,2);
randpos = rand(3,repeatNum);
randori = rand(3,repeatNum);
pcp.func = 'localization';
G = 40;
K = 4;
scalar = 2*sqrt(2);

for j = 1:repeatNum
    x = 20*randpos(1,j) - 10;
    y = 20*randpos(2,j) - 10;
    z = 5*randpos(3,j);
    pu = [x, y, z].';
    ou_euler = 360*randori(:,j);
    disp(['j = ', num2str(j), ',   --------- ', num2str(100*(j/repeatNum)),'% start...']);


    %% planar array, v1
    % default setup
    sp = default_system_setup();
    cp = default_channel_setup();
    % update setup
    cp.G = G;   % # of transmissions
    cp.K = K;   % # of subcarrier
    sp.TypeSA = 'planar';
    sp.vartheta = v(1);
    sp.Pu = pu;
    sp.Ou_Euler = ou_euler;
    sp = update_system_setup(sp);
    cp = update_channel_setup(cp);
    % generate channel
    cp = gen_channel(sp,cp);
    % generate precoders & combiners
    cp = gen_precoder_combiner(sp,cp,pcp);
    % get CRLB
    CCRB = get_CCRB(cp,sp)/scalar;
    data_planar_v1(j,1) = CCRB(1);
    data_planar_v1(j,2) = CCRB(2);

    %% planar array, v2
    % default setup
    sp = default_system_setup();
    cp = default_channel_setup();
    % update setup
    cp.G = G;   % # of transmissions
    cp.K = K;   % # of subcarrier
    sp.TypeSA = 'planar';
    sp.vartheta = v(2);
    sp.Pu = pu;
    sp.Ou_Euler = ou_euler;
    sp = update_system_setup(sp);
    cp = update_channel_setup(cp);
    % generate channel
    cp = gen_channel(sp,cp);
    % generate precoders & combiners
    cp = gen_precoder_combiner(sp,cp,pcp);
    % get CRLB
    CCRB = get_CCRB(cp,sp)/scalar;
    data_planar_v2(j,1) = CCRB(1);
    data_planar_v2(j,2) = CCRB(2);

    %% planar array, v3
    % default setup
    sp = default_system_setup();
    cp = default_channel_setup();
    % update setup
    cp.G = G;   % # of transmissions
    cp.K = K;   % # of subcarrier
    sp.TypeSA = 'planar';
    sp.vartheta = v(3);
    sp.Pu = pu;
    sp.Ou_Euler = ou_euler;
    sp = update_system_setup(sp);
    cp = update_channel_setup(cp);
    % generate channel
    cp = gen_channel(sp,cp);
    % generate precoders & combiners
    cp = gen_precoder_combiner(sp,cp,pcp);
    % get CRLB
    CCRB = get_CCRB(cp,sp)/scalar;
    data_planar_v3(j,1) = CCRB(1);
    data_planar_v3(j,2) = CCRB(2);


    %% cuboidal array, v1
    % default setup
    sp = default_system_setup();
    cp = default_channel_setup();
    % update setup
    cp.G = G;   % # of transmissions
    cp.K = K;   % # of subcarrier
    sp.TypeSA = 'cuboidal';
    sp.vartheta = v(1);
    sp.Pu = pu;
    sp.Ou_Euler = ou_euler;
    sp = update_system_setup(sp);
    cp = update_channel_setup(cp);
    % generate channel
    cp = gen_channel(sp,cp);
    % generate precoders & combiners
    cp = gen_precoder_combiner(sp,cp,pcp);
    % get CRLB
    CCRB = get_CCRB(cp,sp)/scalar;
    data_cuboidal_v1(j,1) = CCRB(1);
    data_cuboidal_v1(j,2) = CCRB(2);

    %% cuboidal array, v2
    % default setup
    sp = default_system_setup();
    cp = default_channel_setup();
    % update setup
    cp.G = G;   % # of transmissions
    cp.K = K;   % # of subcarrier
    sp.TypeSA = 'cuboidal';
    sp.vartheta = v(2);
    sp.Pu = pu;
    sp.Ou_Euler = ou_euler;
    sp = update_system_setup(sp);
    cp = update_channel_setup(cp);
    % generate channel
    cp = gen_channel(sp,cp);
    % generate precoders & combiners
    cp = gen_precoder_combiner(sp,cp,pcp);
    % get CRLB
    CCRB = get_CCRB(cp,sp)/scalar;
    data_cuboidal_v2(j,1) = CCRB(1);
    data_cuboidal_v2(j,2) = CCRB(2);

    %% cuboidal array, v3
    % default setup
    sp = default_system_setup();
    cp = default_channel_setup();
    % update setup
    cp.G = G;   % # of transmissions
    cp.K = K;   % # of subcarrier
    sp.TypeSA = 'cuboidal';
    sp.vartheta = v(3);
    sp.Pu = pu;
    sp.Ou_Euler = ou_euler;
    sp = update_system_setup(sp);
    cp = update_channel_setup(cp);
    % generate channel
    cp = gen_channel(sp,cp);
    % generate precoders & combiners
    cp = gen_precoder_combiner(sp,cp,pcp);
    % get CRLB
    CCRB = get_CCRB(cp,sp)/scalar;
    data_cuboidal_v3(j,1) = CCRB(1);
    data_cuboidal_v3(j,2) = CCRB(2);
end

% data processing
data_planar_v1(isnan(data_planar_v1)) = inf;
data_planar_v2(isnan(data_planar_v2)) = inf;
data_planar_v3(isnan(data_planar_v3)) = inf;
data_cuboidal_v1(isnan(data_cuboidal_v1)) = inf;
data_cuboidal_v2(isnan(data_cuboidal_v2)) = inf;
data_cuboidal_v3(isnan(data_cuboidal_v3)) = inf;

%% plot figures
figure(1)
subplot(2,1,1)
[f,x] =ecdf(data_planar_v1(:,1),'Bounds','on'); 
plot(x,1-f); hold on
[f,x] =ecdf(data_planar_v2(:,1),'Bounds','on'); 
plot(x,1-f); hold on
[f,x] =ecdf(data_planar_v3(:,1),'Bounds','on'); 
plot(x,1-f); hold on
[f,x] = ecdf(data_cuboidal_v1(:,1),'Bounds','on'); 
plot(x,1-f); hold on
[f,x] = ecdf(data_cuboidal_v2(:,1),'Bounds','on'); 
plot(x,1-f); hold on
[f,x] = ecdf(data_cuboidal_v3(:,1),'Bounds','on'); 
plot(x,1-f); hold off
set(gca(),'XScale','log');
set(gca(),'YScale','log');
xlabel('PEB [m]');
ylabel('CCDF');
title('Empirical CDF of PEB');
grid on
legend('2D, V1','2D, V2','2D, V3','3D, V1','3D, V2','3D, V3');



subplot(2,1,2)
[f, x] = ecdf(data_planar_v1(:,2),'Bounds','on'); 
plot(x,1-f); hold on
[f, x] = ecdf(data_planar_v2(:,2),'Bounds','on'); 
plot(x,1-f); hold on
[f, x] = ecdf(data_planar_v3(:,2),'Bounds','on'); 
plot(x,1-f); hold on
[f, x] = ecdf(data_cuboidal_v1(:,2),'Bounds','on'); 
plot(x,1-f); hold on
[f, x] = ecdf(data_cuboidal_v2(:,2),'Bounds','on'); 
plot(x,1-f); hold on
[f, x] = ecdf(data_cuboidal_v3(:,2),'Bounds','on'); 
plot(x,1-f); hold off
set(gca(),'XScale','log');
set(gca(),'YScale','log');
xlabel('OEB [deg]');
ylabel('CCDF');
title('Empirical CDF of OEB');
grid on
legend('2D, V1','2D, V2','2D, V3','3D, V1','3D, V2','3D, V3');

%save('Data_Coverage_CCRB_vDirec.mat','data_planar_v1','data_planar_v2','data_planar_v3','data_cuboidal_v1','data_cuboidal_v2','data_cuboidal_v3')
