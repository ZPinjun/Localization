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
data_planar_2BS = zeros(repeatNum,2);
data_cuboidal_2BS = zeros(repeatNum,2);
data_planar_3BS = zeros(repeatNum,2);
data_cuboidal_3BS = zeros(repeatNum,2);
data_planar_4BS = zeros(repeatNum,2);
data_cuboidal_4BS = zeros(repeatNum,2);
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


    %% planar array, 2 BSs
    % default setup
    sp = default_system_setup();
    cp = default_channel_setup();
    % update setup
    cp.G = G;   % # of transmissions
    cp.K = K;   % # of subcarrier
    sp.TypeSA = 'planar';
    sp.Pu = pu;
    sp.Ou_Euler = ou_euler;
    sp.Pb = [10.5, 10.5, 5; 10.5, -10.5, 5].';
    sp.Ob_Euler = [45 135 0; 90 0 0].';
    sp = update_system_setup(sp);
    cp = update_channel_setup(cp);
    % generate channel
    cp = gen_channel(sp,cp);
    % generate precoders & combiners
    cp = gen_precoder_combiner(sp,cp,pcp);
    % get CRLB
    CCRB = get_CCRB(cp,sp)/scalar;
    data_planar_2BS(j,1) = CCRB(1);
    data_planar_2BS(j,2) = CCRB(2);

    %% planar array, 3 BSs
    % default setup
    sp = default_system_setup();
    cp = default_channel_setup();
    % update setup
    cp.G = G;   % # of transmissions
    cp.K = K;   % # of subcarrier
    sp.TypeSA = 'planar';
    sp.Pu = pu;
    sp.Ou_Euler = ou_euler;
    sp.Pb = [10.5, 10.5, 5; 10.5, -10.5, 5; -10.5, 10.5, 5].';
    sp.Ob_Euler = [45 135 0; 90 0 0; -45 45 0].';
    sp = update_system_setup(sp);
    cp = update_channel_setup(cp);
    % generate channel
    cp = gen_channel(sp,cp);
    % generate precoders & combiners
    cp = gen_precoder_combiner(sp,cp,pcp);
    % get CRLB
    CCRB = get_CCRB(cp,sp)/scalar;
    data_planar_3BS(j,1) = CCRB(1);
    data_planar_3BS(j,2) = CCRB(2);

    %% planar array, 4 BSs
    % default setup
    sp = default_system_setup();
    cp = default_channel_setup();
    % update setup
    cp.G = G;   % # of transmissions
    cp.K = K;   % # of subcarrier
    sp.TypeSA = 'planar';
    sp.Pu = pu;
    sp.Ou_Euler = ou_euler;
    sp.Pb = [10.5, 10.5, 5; 10.5, -10.5, 5; -10.5, 10.5, 5; -10.5, -10.5, 5].';
    sp.Ob_Euler = [45 135 0; 90 0 0; -45 45 0; 45 45 0].';
    sp = update_system_setup(sp);
    cp = update_channel_setup(cp);
    % generate channel
    cp = gen_channel(sp,cp);
    % generate precoders & combiners
    cp = gen_precoder_combiner(sp,cp,pcp);
    % get CRLB
    CCRB = get_CCRB(cp,sp)/scalar;
    data_planar_4BS(j,1) = CCRB(1);
    data_planar_4BS(j,2) = CCRB(2);


    %% cuboidal array, 2 BSs
    % default setup
    sp = default_system_setup();
    cp = default_channel_setup();
    % update setup
    cp.G = G;   % # of transmissions
    cp.K = K;   % # of subcarrier
    sp.TypeSA = 'cuboidal';
    sp.Pu = pu;
    sp.Ou_Euler = ou_euler;
    sp.Pb = [10.5, 10.5, 5; 10.5, -10.5, 5].';
    sp.Ob_Euler = [45 135 0; 90 0 0].';
    sp = update_system_setup(sp);
    cp = update_channel_setup(cp);
    % generate channel
    cp = gen_channel(sp,cp);
    % generate precoders & combiners
    cp = gen_precoder_combiner(sp,cp,pcp);
    % get CRLB
    CCRB = get_CCRB(cp,sp)/scalar;
    data_cuboidal_2BS(j,1) = CCRB(1);
    data_cuboidal_2BS(j,2) = CCRB(2);

    %% cuboidal array, 3 BSs
    % default setup
    sp = default_system_setup();
    cp = default_channel_setup();
    % update setup
    cp.G = G;   % # of transmissions
    cp.K = K;   % # of subcarrier
    sp.TypeSA = 'cuboidal';
    sp.Pu = pu;
    sp.Ou_Euler = ou_euler;
    sp.Pb = [10.5, 10.5, 5; 10.5, -10.5, 5; -10.5, 10.5, 5].';
    sp.Ob_Euler = [45 135 0; 90 0 0; -45 45 0].';
    sp = update_system_setup(sp);
    cp = update_channel_setup(cp);
    % generate channel
    cp = gen_channel(sp,cp);
    % generate precoders & combiners
    cp = gen_precoder_combiner(sp,cp,pcp);
    % get CRLB
    CCRB = get_CCRB(cp,sp)/scalar;
    data_cuboidal_3BS(j,1) = CCRB(1);
    data_cuboidal_3BS(j,2) = CCRB(2);

    %% cuboidal array, 4 BSs
    % default setup
    sp = default_system_setup();
    cp = default_channel_setup();
    % update setup
    cp.G = G;   % # of transmissions
    cp.K = K;   % # of subcarrier
    sp.TypeSA = 'cuboidal';
    sp.Pu = pu;
    sp.Ou_Euler = ou_euler;
    sp.Pb = [10.5, 10.5, 5; 10.5, -10.5, 5; -10.5, 10.5, 5; -10.5, -10.5, 5].';
    sp.Ob_Euler = [45 135 0; 90 0 0; -45 45 0; 45 45 0].';
    sp = update_system_setup(sp);
    cp = update_channel_setup(cp);
    % generate channel
    cp = gen_channel(sp,cp);
    % generate precoders & combiners
    cp = gen_precoder_combiner(sp,cp,pcp);
    % get CRLB
    CCRB = get_CCRB(cp,sp)/scalar;
    data_cuboidal_4BS(j,1) = CCRB(1);
    data_cuboidal_4BS(j,2) = CCRB(2);
end

% data processing
data_planar_2BS(isnan(data_planar_2BS)) = inf;
data_planar_3BS(isnan(data_planar_3BS)) = inf;
data_planar_4BS(isnan(data_planar_4BS)) = inf;
data_cuboidal_2BS(isnan(data_cuboidal_2BS)) = inf;
data_cuboidal_3BS(isnan(data_cuboidal_3BS)) = inf;
data_cuboidal_4BS(isnan(data_cuboidal_4BS)) = inf;

% plot figures
figure(1)
subplot(2,1,1)
[f,x] =ecdf(data_planar_2BS(:,1),'Bounds','on'); 
plot(x,1-f); hold on
[f,x] =ecdf(data_planar_3BS(:,1),'Bounds','on'); 
plot(x,1-f); hold on
[f,x] =ecdf(data_planar_4BS(:,1),'Bounds','on'); 
plot(x,1-f); hold on
[f,x] = ecdf(data_cuboidal_2BS(:,1),'Bounds','on'); 
plot(x,1-f); hold on
[f,x] = ecdf(data_cuboidal_3BS(:,1),'Bounds','on'); 
plot(x,1-f); hold on
[f,x] = ecdf(data_cuboidal_4BS(:,1),'Bounds','on'); 
plot(x,1-f); hold off
set(gca(),'XScale','log');
set(gca(),'YScale','log');
xlabel('PEB [m]');
ylabel('CCDF');
title('Empirical CDF of PEB');
xlim([1e-3, 1e4]);
grid on
legend('2D, 2BS','2D, 3BS','2D, 4BS','3D, 2BS','3D, 3BS','3D, 4BS');
yticks([1e-4,1e-3,1e-2,1e-1,1]);

subplot(2,1,2)
[f, x] = ecdf(data_planar_2BS(:,2),'Bounds','on'); 
plot(x,1-f); hold on
[f, x] = ecdf(data_planar_3BS(:,2),'Bounds','on'); 
plot(x,1-f); hold on
[f, x] = ecdf(data_planar_4BS(:,2),'Bounds','on'); 
plot(x,1-f); hold on
[f, x] = ecdf(data_cuboidal_2BS(:,2),'Bounds','on'); 
plot(x,1-f); hold on
[f, x] = ecdf(data_cuboidal_3BS(:,2),'Bounds','on'); 
plot(x,1-f); hold on
[f, x] = ecdf(data_cuboidal_4BS(:,2),'Bounds','on'); 
plot(x,1-f); hold off
set(gca(),'XScale','log');
set(gca(),'YScale','log');
xlabel('OEB [deg]');
ylabel('CCDF');
title('Empirical CDF of OEB');
xlim([1e-3, 1e4]);
yticks([1e-4,1e-3,1e-2,1e-1,1]);
grid on
legend('2D, 2BS','2D, 3BS','2D, 4BS','3D, 2BS','3D, 3BS','3D, 4BS');

%save('Data_Coverage_CCRB_vBS.mat','data_planar_2BS','data_planar_3BS','data_planar_4BS','data_cuboidal_2BS','data_cuboidal_3BS','data_cuboidal_4BS')

