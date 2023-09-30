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
pcp.func = 'localization';


%% simulation for planar array
data_planar_peb = zeros(length(Pu_range_x),length(Pu_range_y));
data_planar_oeb = zeros(length(Pu_range_x),length(Pu_range_y));
for i = 1:length(Pu_range_x)
    parfor j = 1:length(Pu_range_y)
        disp(['**Planar**  ', 'x = ', num2str(Pu_range_x(i)),', y = ', num2str(Pu_range_y(j)), '   start...']);

        % default setup
        sp = default_system_setup();
        cp = default_channel_setup();
        % update setup
        cp.G = 320;   % # of transmissions
        cp.K = 4;   % # of subcarrier
        sp.TypeSA = 'planar';
        sp.Pu = [Pu_range_x(i), Pu_range_y(j), Pu_z].';
        sp = update_system_setup(sp);
        cp = update_channel_setup(cp);
        % generate channel
        cp = gen_channel(sp,cp);
        % generate precoders & combiners
        cp = gen_precoder_combiner(sp,cp,pcp);
        % get CRLB
        CCRB = get_CCRB(cp,sp);

        data_planar_peb(i,j) = CCRB(1);
        data_planar_oeb(i,j) = CCRB(2);
    end
end



%% simulation for cuboidal array
data_cuboidal_peb = zeros(length(Pu_range_x),length(Pu_range_y));
data_cuboidal_oeb = zeros(length(Pu_range_x),length(Pu_range_y));
for i = 1:length(Pu_range_x)
    parfor j = 1:length(Pu_range_y)
        disp(['**Cuboidal**', 'x = ', num2str(Pu_range_x(i)),', y = ', num2str(Pu_range_y(j)), '   start...']);
        
        % default setup
        sp = default_system_setup();
        cp = default_channel_setup();
        % update setup
        cp.G = 320;   % # of transmissions
        cp.K = 4;   % # of subcarrier
        sp.TypeSA = 'cuboidal';
        sp.Pu = [Pu_range_x(i), Pu_range_y(j), Pu_z].';
        sp = update_system_setup(sp);
        cp = update_channel_setup(cp);
        % generate channel
        cp = gen_channel(sp,cp);
        % generate precoders & combiners
        cp = gen_precoder_combiner(sp,cp,pcp);
        % get CRLB
        CCRB = get_CCRB(cp,sp);

        data_cuboidal_peb(i,j) = CCRB(1);
        data_cuboidal_oeb(i,j) = CCRB(2);
    end
end

%% plot figs
figure(1)
subplot(2,2,1)
imagesc(Pu_range_x, Pu_range_y, data_planar_peb); hold on
title('PEB of planar array [m]');
colorbar
xlabel('y [m]');
ylabel('x [m]');
subplot(2,2,2)
imagesc(Pu_range_x, Pu_range_y, data_planar_oeb); hold on
title('OEB of planar array [deg]');
colorbar
% caxis([0 0.04]);
xlabel('y [m]');
ylabel('x [m]');
subplot(2,2,3)
imagesc(Pu_range_x, Pu_range_y, data_cuboidal_peb); hold on
title('PEB of cuboidal array [m]');
colorbar
xlabel('y [m]');
ylabel('x [m]');
subplot(2,2,4)
imagesc(Pu_range_x, Pu_range_y, data_cuboidal_oeb); hold on
title('OEB of cuboidal array [deg]');
colorbar
xlabel('y [m]');
ylabel('x [m]');



%save('Data_CCRB_vs_pos.mat','Pu_range_x','Pu_range_y','data_planar_peb','data_planar_oeb','data_cuboidal_peb','data_cuboidal_oeb');


