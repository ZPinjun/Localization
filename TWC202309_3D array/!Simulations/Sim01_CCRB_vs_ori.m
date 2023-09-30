clc;clear;close all

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
pcp.func = 'localization';



data_planar_peb = zeros(length(Ru_range_y),length(Ru_range_z));
data_planar_oeb = zeros(length(Ru_range_y),length(Ru_range_z));

%% simulation for planar array
for i = 1:length(Ru_range_y)
    parfor j = 1:length(Ru_range_z)
        disp(['**Planar**  ', 'beta = ', num2str(Ru_range_y(i)),', gamma = ', num2str(Ru_range_z(j)), '   start...']);
        
        % default setup
        sp = default_system_setup();
        cp = default_channel_setup();
        % update setup
        cp.G = 320;   % # of transmissions
        cp.K = 4;   % # of subcarrier
        sp.TypeSA = 'planar';
        sp.Ou_Euler = [ Ru_range_z(j) Ru_range_y(i) Ru_x ].';
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



data_cuboidal_peb = zeros(length(Ru_range_y),length(Ru_range_z));
data_cuboidal_oeb = zeros(length(Ru_range_y),length(Ru_range_z));

%% simulation for cuboidal array
for i = 1:length(Ru_range_y)
    parfor j = 1:length(Ru_range_z)
        disp(['**Cuboidal**  ', 'beta = ', num2str(Ru_range_y(i)),', gamma = ', num2str(Ru_range_z(j)), '   start...']);
        
        % default setup
        sp = default_system_setup();
        cp = default_channel_setup();
        % update setup
        cp.G = 320;   % # of transmissions
        cp.K = 4;   % # of subcarrier
        sp.TypeSA = 'cuboidal';
        sp.Ou_Euler = [ Ru_range_z(j) Ru_range_y(i) Ru_x ].';
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
imagesc(Ru_range_y, Ru_range_z, data_planar_peb); hold on
title('PEB of planar array [m]');
colorbar
xlabel('\gamma [deg]');
ylabel('\beta [deg]');
subplot(2,2,2)
imagesc(Ru_range_y, Ru_range_z, data_planar_oeb); hold on
title('OEB of planar array [deg]');
colorbar
xlabel('\gamma [deg]');
ylabel('\beta [deg]');
subplot(2,2,3)
imagesc(Ru_range_y, Ru_range_z, data_cuboidal_peb); hold on
title('PEB of cuboidal array [m]');
colorbar
xlabel('\gamma [deg]');
ylabel('\beta [deg]');
subplot(2,2,4)
imagesc(Ru_range_y, Ru_range_z, data_cuboidal_oeb); hold on
title('OEB of cuboidal array [deg]');
colorbar
xlabel('\gamma [deg]');
ylabel('\beta [deg]');



%save('Data_CCRB_vs_ori.mat','Ru_range_y','Ru_range_z','data_planar_peb','data_planar_oeb','data_cuboidal_peb','data_cuboidal_oeb');



