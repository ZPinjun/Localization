clc; clear; close all

%% Add Paths
path(pathdef); addpath(pwd);
cd ..;
cd Channel; addpath(genpath(pwd)); cd ..;
cd System; addpath(genpath(pwd)); cd ..; 
cd Localization; addpath(genpath(pwd)); cd ..; 
cd !Simulations;



%% simulation setup
Pu_range_x = -10:0.05:10;
Pu_range_y = -10:0.05:10;
Pu_z = 0;
pcp.func = 'communication';
threshold = 100;

%% simulation for planar array
data_planar_pout = zeros(length(Pu_range_x),length(Pu_range_y));
data_planar_sumrate = zeros(length(Pu_range_x),length(Pu_range_y));
data_planar_mtilde = zeros(length(Pu_range_x),length(Pu_range_y));
for i = 1:length(Pu_range_x)
    tic
    parfor j = 1:length(Pu_range_y)

        disp(['**Planar**  ', 'x = ', num2str(Pu_range_x(i)),', y = ', num2str(Pu_range_y(j)), '   start...']);

        % default setup
        sp = default_system_setup();
        cp = default_channel_setup();
        % update setup
        sp.TypeSA = 'planar';
        sp.Pu = [Pu_range_x(i), Pu_range_y(j), Pu_z].';
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
        data_planar_pout(i,j) = prod(get_outage_probability(sp,cp,threshold,k,method));
        data_planar_sumrate(i,j) = cp.m_sumrate;
        data_planar_mtilde(i,j) = cp.m_tilde;
    end
    toc
end


%% simulation for cuboidal array
data_cuboidal_pout = zeros(length(Pu_range_x),length(Pu_range_y));
data_cuboidal_sumrate = zeros(length(Pu_range_x),length(Pu_range_y));
data_cuboidal_mtilde = zeros(length(Pu_range_x),length(Pu_range_y));
for i = 1:length(Pu_range_x)
    parfor j = 1:length(Pu_range_y)

        disp(['**Cuboidal**  ', 'x = ', num2str(Pu_range_x(i)),', y = ', num2str(Pu_range_y(j)), '   start...']);

        % default setup
        sp = default_system_setup();
        cp = default_channel_setup();
        % update setup
        sp.TypeSA = 'cuboidal';
        sp.Pu = [Pu_range_x(i), Pu_range_y(j), Pu_z].';
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
        data_cuboidal_pout(i,j) = prod(get_outage_probability(sp,cp,threshold,k,method));
        data_cuboidal_sumrate(i,j) = cp.m_sumrate;
        data_cuboidal_mtilde(i,j) = cp.m_tilde;
    end
end

%% plot figs
figure(1)
subplot(1,2,1)
imagesc(Pu_range_x, Pu_range_y, data_planar_sumrate); hold on
title('sum-rate');
colorbar
xlabel('y [m]');
ylabel('x [m]');

subplot(1,2,2)
imagesc(Pu_range_x, Pu_range_y, data_planar_mtilde); hold on
title('m_tilde');
colorbar
xlabel('y [m]');
ylabel('x [m]');


figure(2)
subplot(1,2,1)
imagesc(Pu_range_x, Pu_range_y, data_cuboidal_sumrate); hold on
title('sum-rate');
colorbar
xlabel('y [m]');
ylabel('x [m]');

subplot(1,2,2)
imagesc(Pu_range_x, Pu_range_y, data_cuboidal_mtilde); hold on
title('m_tilde');
colorbar
xlabel('y [m]');
ylabel('x [m]');


figure(3)
subplot(1,2,1)
imagesc(Pu_range_x, Pu_range_y, data_planar_pout); hold on
title('Outage probability of planar array [m]');
colorbar
xlabel('y [m]');
ylabel('x [m]');

subplot(1,2,2)
imagesc(Pu_range_x, Pu_range_y, data_cuboidal_pout); hold on
title('Outage probability of cuboidal array [m]');
colorbar
xlabel('y [m]');
ylabel('x [m]');


% save Data_Pout_vs_pos.mat
