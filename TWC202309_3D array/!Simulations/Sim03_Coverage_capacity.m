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
KrRange = [1,4,16];
data_planar_v1 = zeros(repeatNum,1);
data_cuboidal_v1 = zeros(repeatNum,1);
data_planar_v2 = zeros(repeatNum,1);
data_cuboidal_v2 = zeros(repeatNum,1);
data_planar_v3 = zeros(repeatNum,1);
data_cuboidal_v3 = zeros(repeatNum,1);
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


    %% planar array, v1
    % default setup
    sp = default_system_setup();
    cp = default_channel_setup();
    % update setup
    cp.G = 1;   % # of transmissions
    cp.K = 2;   % to reduce running time
    sp.TypeSA = 'planar';
    sp.Pu = pu;
    sp.Ou_Euler = ou_euler;
    %     sp.vartheta = vartheta(1);
    cp.Kr = KrRange(1);
    sp = update_system_setup(sp);
    cp = update_channel_setup(cp);
    % generate channel
    cp = gen_channel(sp,cp);
    % generate precoders & combiners
    cp = gen_precoder_combiner(sp,cp,pcp);
    % get ergodic capacity
    C = get_ergodic_capacity(sp,cp,20);
    data_planar_v1(j) = C;

    %% planar array, v2
    % default setup
    sp = default_system_setup();
    cp = default_channel_setup();
    % update setup
    cp.G = 1;   % # of transmissions
    cp.K = 2;   % to reduce running time
    sp.TypeSA = 'planar';
    sp.Pu = pu;
    sp.Ou_Euler = ou_euler;
    %     sp.vartheta = vartheta(2);
    cp.Kr = KrRange(2);
    sp = update_system_setup(sp);
    cp = update_channel_setup(cp);
    % generate channel
    cp = gen_channel(sp,cp);
    % generate precoders & combiners
    cp = gen_precoder_combiner(sp,cp,pcp);
    % get ergodic capacity
    C = get_ergodic_capacity(sp,cp,20)*128/2;
    data_planar_v2(j) = C;

    %% planar array, v3
    % default setup
    sp = default_system_setup();
    cp = default_channel_setup();
    % update setup
    cp.G = 1;   % # of transmissions
    cp.K = 2;   % to reduce running time
    sp.TypeSA = 'planar';
    sp.Pu = pu;
    sp.Ou_Euler = ou_euler;
    %     sp.vartheta = vartheta(3);
    cp.Kr = KrRange(3);
    sp = update_system_setup(sp);
    cp = update_channel_setup(cp);
    % generate channel
    cp = gen_channel(sp,cp);
    % generate precoders & combiners
    cp = gen_precoder_combiner(sp,cp,pcp);
    % get ergodic capacity
    C = get_ergodic_capacity(sp,cp,20)*128/2;
    data_planar_v3(j) = C;



    %% cuboidal array, v1
    % default setup
    sp = default_system_setup();
    cp = default_channel_setup();
    % update setup
    cp.G = 1;   % # of transmissions
    cp.K = 2;   % to reduce running time
    sp.TypeSA = 'cuboidal';
    sp.Pu = pu;
    sp.Ou_Euler = ou_euler;
    %     sp.vartheta = vartheta(1);
    cp.Kr = KrRange(1);
    sp = update_system_setup(sp);
    cp = update_channel_setup(cp);
    % generate channel
    cp = gen_channel(sp,cp);
    % generate precoders & combiners
    cp = gen_precoder_combiner(sp,cp,pcp);
    % get ergodic capacity
    C = get_ergodic_capacity(sp,cp,20)*128/2;
    data_cuboidal_v1(j) = C;

    %% cuboidal array, v2
    % default setup
    sp = default_system_setup();
    cp = default_channel_setup();
    % update setup
    cp.G = 1;   % # of transmissions
    cp.K = 2;   % to reduce running time
    sp.TypeSA = 'cuboidal';
    sp.Pu = pu;
    sp.Ou_Euler = ou_euler;
    %     sp.vartheta = vartheta(2);
    cp.Kr = KrRange(2);
    sp = update_system_setup(sp);
    cp = update_channel_setup(cp);
    % generate channel
    cp = gen_channel(sp,cp);
    % generate precoders & combiners
    cp = gen_precoder_combiner(sp,cp,pcp);
    % get ergodic capacity
    C = get_ergodic_capacity(sp,cp,20)*128/2;
    data_cuboidal_v2(j) = C;

    %% cuboidal array, v3
    % default setup
    sp = default_system_setup();
    cp = default_channel_setup();
    % update setup
    cp.G = 1;   % # of transmissions
    cp.K = 2;   % to reduce running time
    sp.TypeSA = 'cuboidal';
    sp.Pu = pu;
    sp.Ou_Euler = ou_euler;
    %     sp.vartheta = vartheta(3);
    cp.Kr = KrRange(3);
    sp = update_system_setup(sp);
    cp = update_channel_setup(cp);
    % generate channel
    cp = gen_channel(sp,cp);
    % generate precoders & combiners
    cp = gen_precoder_combiner(sp,cp,pcp);
    % get ergodic capacity
    C = get_ergodic_capacity(sp,cp,20)*128/2;
    data_cuboidal_v3(j) = C;


end

figure(1)
[f,x] =ecdf(data_planar_v1(:,1),'Bounds','on'); 
plot(x,f,'b--'); hold on
[f,x] =ecdf(data_planar_v2(:,1),'Bounds','on'); 
plot(x,f,'b--'); hold on
[f,x] =ecdf(data_planar_v3(:,1),'Bounds','on'); 
plot(x,f,'b--'); hold on
[f,x] = ecdf(data_cuboidal_v1(:,1),'Bounds','on'); 
plot(x,f); hold on
[f,x] = ecdf(data_cuboidal_v2(:,1),'Bounds','on'); 
plot(x,f); hold on
[f,x] = ecdf(data_cuboidal_v3(:,1),'Bounds','on'); 
plot(x,f); hold off
set(gca(),'XScale','log');
xlabel('Ergodic capacity threshold [bps/Hz]');
ylabel('Ergodic capacity coverage');
title('Ergodic capacity coverage vs. threshold');
%xlim([1, 100]);
ylim([1e-2, 1]);
grid on
legend('Planar,K1','Planar,K2','Planar,K3','Cuboidal,K1','Cuboidal,K2','Cuboidal,K3');



% save('Data_Coverage_capacity.mat','data_planar_v1','data_planar_v2','data_planar_v3','data_cuboidal_v1','data_cuboidal_v2','data_cuboidal_v3')
