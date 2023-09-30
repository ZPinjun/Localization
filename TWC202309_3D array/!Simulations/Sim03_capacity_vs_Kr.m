clc; clear; close all

%% Add Paths
path(pathdef); addpath(pwd);
cd ..;
cd Channel; addpath(genpath(pwd)); cd ..;
cd System; addpath(genpath(pwd)); cd ..; 
cd Localization; addpath(genpath(pwd)); cd ..; 
cd !Simulations;


%% simulation setup
RangeKrdB = -20:6:31;   % [dB]
RangeKr = 10.^(RangeKrdB/10); 
RangeP = [5, 10, 15];
data_capacity_planar = zeros(3,length(RangeKr));
data_capacity_cuboidal = zeros(3,length(RangeKr));
pcp.func = 'communication';
repeatNum = 50;

for i = 1:length(RangeKr)
    disp(['Kr = ', num2str(RangeKrdB(i)), ' [dB],   --------- ', num2str(100*(i/length(RangeKr))),'% start...']);
    %% planar array
    parfor j = 1:length(RangeP)
        % default setup
        sp = default_system_setup();
        cp = default_channel_setup();
        % update setup
        cp.G = 1;   % # of transmissions
        sp.TypeSA = 'planar';
        cp.Kr = RangeKr(i);
        cp.P = RangeP(j);
        sp = update_system_setup(sp);
        cp = update_channel_setup(cp);
        % generate channel
        cp = gen_channel(sp,cp);
        % generate precoders & combiners
        cp = gen_precoder_combiner(sp,cp,pcp);
        % get ergodic capacity
        C = get_ergodic_capacity(sp,cp,repeatNum);
        data_capacity_planar(j,i) = C;
    end

    %% cuboidal array
    parfor j = 1:length(RangeP)
        % default setup
        sp = default_system_setup();
        cp = default_channel_setup();
        % update setup
        cp.G = 1;   % # of transmissions
        sp.TypeSA = 'cuboidal';
        cp.Kr = RangeKr(i);
        cp.P = RangeP(j);
        sp = update_system_setup(sp);
        cp = update_channel_setup(cp);
        % generate channel
        cp = gen_channel(sp,cp);
        % generate precoders & combiners
        cp = gen_precoder_combiner(sp,cp,pcp);
        % get ergodic capacity
        C = get_ergodic_capacity(sp,cp,repeatNum);
        data_capacity_cuboidal(j,i) = C
    end
end


%% plot figs
figure(2)
semilogx(RangeKr,data_capacity_planar(3,:),'b^--'); hold on
semilogx(RangeKr,data_capacity_planar(2,:),'bs--'); hold on
semilogx(RangeKr,data_capacity_planar(1,:),'bv--'); hold on
semilogx(RangeKr,data_capacity_cuboidal(3,:),'r^-'); hold on
semilogx(RangeKr,data_capacity_cuboidal(2,:),'rs-'); hold on
semilogx(RangeKr,data_capacity_cuboidal(1,:),'rv-'); hold off
legend('Planar,P=15mW','Planar,P=10mW','Planar,P=5mW','Cuboidal,P=15mW','Cuboidal,P=10mW','Cuboidal,P=5mW');
xlabel('Kr');
ylabel('Ergodic capacity [bps/Hz]');

% save('Data_capacity_vs_Kr.mat','RangeKr','RangeP','data_capacity_planar','data_capacity_cuboidal');