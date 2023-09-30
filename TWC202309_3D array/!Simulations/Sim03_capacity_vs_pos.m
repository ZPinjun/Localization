clc; clear; close all

%% Add Paths
path(pathdef); addpath(pwd);
cd ..;
cd Channel; addpath(genpath(pwd)); cd ..;
cd System; addpath(genpath(pwd)); cd ..; 
cd Localization; addpath(genpath(pwd)); cd ..; 
cd !Simulations;



%% simulation setup
Pu_range_x = -10:0.2:10;
Pu_range_y = -10:0.2:10;
Pu_z = 0;
pcp.func = 'communication';
repeatNum = 20;


%% simulation for planar array
data_planar_capacity = zeros(length(Pu_range_x),length(Pu_range_y));
for i = 1:length(Pu_range_x)
    parfor j = 1:length(Pu_range_y)

        disp(['**Planar**  ', 'x = ', num2str(Pu_range_x(i)),', y = ', num2str(Pu_range_y(j)), '   start...']);

        % default setup
        sp = default_system_setup();
        cp = default_channel_setup();
        % update setup
        cp.G = 1;   % # of transmissions
        cp.K = 2;   % to reduce running time
        sp.TypeSA = 'planar';
        sp.Pu = [Pu_range_x(i), Pu_range_y(j), Pu_z].';
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
data_cuboidal_capacity = zeros(length(Pu_range_x),length(Pu_range_y));
for i = 1:length(Pu_range_x)
    parfor j = 1:length(Pu_range_y)

        disp(['**Cuboidal**  ', 'x = ', num2str(Pu_range_x(i)),', y = ', num2str(Pu_range_y(j)), '   start...']);

        % default setup
        sp = default_system_setup();
        cp = default_channel_setup();
        % update setup
        cp.G = 1;   % # of transmissions
        cp.K = 2;   % to reduce running time
        sp.TypeSA = 'cuboidal';
        sp.Pu = [Pu_range_x(i), Pu_range_y(j), Pu_z].';
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
figure(1)
subplot(1,2,1)
imagesc(Pu_range_x, Pu_range_y, data_planar_capacity); hold on
%contour(Pu_range_x, Pu_range_y, data_planar_pout, 'ShowText','on', 'LineColor', 'w');
title('Ergodic capacity of planar array [m]');
colorbar
xlabel('y [m]');
ylabel('x [m]');

subplot(1,2,2)
imagesc(Pu_range_x, Pu_range_y, data_cuboidal_capacity); hold on
%contour(Pu_range_x, Pu_range_y, data_planar_pout, 'ShowText','on', 'LineColor', 'w');
title('Ergodic capacity of cuboidal array [m]');
colorbar
xlabel('y [m]');
ylabel('x [m]');

% save('Data_capacity_vs_pos.mat','Pu_range_x','Pu_range_y','data_planar_capacity','data_cuboidal_capacity');

