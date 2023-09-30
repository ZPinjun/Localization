clc; clear; close all

%% Add Paths
path(pathdef); addpath(pwd);
cd ..;
cd Channel; addpath(genpath(pwd)); cd ..;
cd System; addpath(genpath(pwd)); cd ..; 
cd Localization; addpath(genpath(pwd)); cd ..; 
cd !Simulations;


%% simulation setup
Kr_test = 10.^((-40:10:20)/10);
repeatNum = 10;
pcp.func = 'localization';


%% simulation for planar array
data_planar_peb = zeros(repeatNum,length(Kr_test));
data_planar_oeb = zeros(repeatNum,length(Kr_test));
for i = 1:length(Kr_test)
    parfor j = 1:repeatNum
        disp(['**Planar**  ', 'Kr = ', num2str(Kr_test(i)), '   start...']);
    
        % default setup
        sp = default_system_setup();
        cp = default_channel_setup();
        % update setup
        sp.TypeSA = 'planar';
        cp.Kr = Kr_test(i);
        cp.G = 640;   % # of transmissions
        cp.K = 2;   % # of subcarrier
        sp = update_system_setup(sp);
        cp = update_channel_setup(cp);
        % generate channel
        cp = gen_channel(sp,cp);
        % generate precoders & combiners
        cp = gen_precoder_combiner(sp,cp,pcp);
        % get CRLB
        CCRB = get_CCRB(cp,sp);
    
        data_planar_peb(j,i) = CCRB(1);
        data_planar_oeb(j,i) = CCRB(2);
    end
end


%% simulation for planar array
data_cuboidal_peb = zeros(repeatNum,length(Kr_test));
data_cuboidal_oeb = zeros(repeatNum,length(Kr_test));
for i = 1:length(Kr_test)
    parfor j = 1:repeatNum
        disp(['**Cuboidal**  ', 'Kr = ', num2str(Kr_test(i)),  '   start...']);
        % default setup
        sp = default_system_setup();
        cp = default_channel_setup();
        % update setup
        sp.TypeSA = 'cuboidal';
        cp.Kr = Kr_test(i);
        cp.G = 640;   % # of transmissions
        cp.K = 2;   % # of subcarrier
        sp = update_system_setup(sp);
        cp = update_channel_setup(cp);
        % generate channel
        cp = gen_channel(sp,cp);
        % generate precoders & combiners
        cp = gen_precoder_combiner(sp,cp,pcp);
        % get CRLB
        CCRB = get_CCRB(cp,sp);

        data_cuboidal_peb(j,i) = CCRB(1);
        data_cuboidal_oeb(j,i) = CCRB(2);
    end
end


% data average
data_planar_peb = mean(data_planar_peb);
data_planar_oeb = mean(data_planar_oeb);
data_cuboidal_peb = mean(data_cuboidal_peb);
data_cuboidal_oeb = mean(data_cuboidal_oeb);


%% plot figs
figure(1)
subplot(2,1,1)
plot(Kr_test, data_planar_peb(1,:),'bo-'); hold on
plot(Kr_test, data_cuboidal_peb(1,:),'rs-'); hold off
set(gca(),'XScale','log');
set(gca(),'YScale','log');
title('PEB vs. Kr');
xlabel('Kr');
ylabel('PEB');
legend('Planar','Cuboidal')
grid on

subplot(2,1,2)
plot(Kr_test, data_planar_oeb(1,:),'bo-'); hold on
plot(Kr_test, data_cuboidal_oeb(1,:),'rs-'); hold off
set(gca(),'XScale','log');
set(gca(),'YScale','log');
title('OEB vs. Kr');
xlabel('Kr');
ylabel('OEB');
legend('Planar','Cuboidal')
grid on

%save('Data_CCRB_vs_Kr.mat','Kr_test','data_planar_peb','data_cuboidal_peb','data_planar_oeb','data_cuboidal_oeb')

