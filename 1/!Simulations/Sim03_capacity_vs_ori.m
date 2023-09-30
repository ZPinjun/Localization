clc; clear; close all

%% Add Paths
path(pathdef); addpath(pwd);
cd ..;
cd Channel; addpath(genpath(pwd)); cd ..;
cd System; addpath(genpath(pwd)); cd ..; 
cd Localization; addpath(genpath(pwd)); cd ..; 
cd !Simulations;



%% simulation setup
Ru_x = 0;
Ru_range_y = 0:3:360;
Ru_range_z = 0:3:360;
pcp.func = 'communication';
repeatNum = 20;


%% simulation for planar array
data_planar_capacity = zeros(length(Ru_range_y),length(Ru_range_z));
for i = 1:length(Ru_range_y)
    parfor j = 1:length(Ru_range_z)

        disp(['**Planar**  ', 'beta = ', num2str(Ru_range_y(i)),', gamma = ', num2str(Ru_range_z(j)), '   start...']);
    
        % default setup
        sp = default_system_setup();
        cp = default_channel_setup();
        % update setup
        cp.G = 1;   % # of transmissions
        cp.K = 2;   % to reduce running time
        sp.TypeSA = 'planar';
        sp.Ou_Euler = [ Ru_range_z(j) Ru_range_y(i) Ru_x ].';
        sp = update_system_setup(sp);
        cp = update_channel_setup(cp);
        % generate channel
        cp = gen_channel(sp,cp);
        % generate precoders & combiners
        cp = gen_precoder_combiner(sp,cp,pcp);
        % get ergodic capacity
        C = get_ergodic_capacity(sp,cp,repeatNum);
        data_planar_capacity(i,j) = C;
    end
end


%% simulation for cuboidal array
data_cuboidal_capacity = zeros(length(Ru_range_y),length(Ru_range_z));
for i = 1:length(Ru_range_y)
    parfor j = 1:length(Ru_range_z)

        disp(['**Cuboidal**  ', 'beta = ', num2str(Ru_range_y(i)),', gamma = ', num2str(Ru_range_z(j)), '   start...']);

        % default setup
        sp = default_system_setup();
        cp = default_channel_setup();
        % update setup
        cp.G = 1;   % # of transmissions
        cp.K = 2;   % to reduce running time
        sp.TypeSA = 'cuboidal';
        sp.Ou_Euler = [ Ru_range_z(j) Ru_range_y(i) Ru_x ].';
        sp = update_system_setup(sp);
        cp = update_channel_setup(cp);
        % generate channel
        cp = gen_channel(sp,cp);
        % generate precoders & combiners
        cp = gen_precoder_combiner(sp,cp,pcp);
        % get ergodic capacity
        C = get_ergodic_capacity(sp,cp,repeatNum);
        data_cuboidal_capacity(i,j) = C;
    end
end




%% plot figs
figure(2)
subplot(1,2,1)
imagesc(Ru_range_y, Ru_range_z, data_planar_capacity); hold on
title('Ergodic capacity of planar array [m]');
colorbar
xlabel('\gamma [deg]');
ylabel('\beta [deg]');

subplot(1,2,2)
imagesc(Ru_range_y, Ru_range_z, data_cuboidal_capacity); hold on
title('Ergodic capacity of cuboidal array [m]');
colorbar
xlabel('\gamma [deg]');
ylabel('\beta [deg]');

% save('Data_capacity_vs_ori.mat','Ru_range_y','Ru_range_z','data_planar_capacity','data_cuboidal_capacity');

