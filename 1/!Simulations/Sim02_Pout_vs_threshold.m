clc;clear;close all

%% Add Paths
path(pathdef); addpath(pwd);
cd ..;
cd Channel; addpath(genpath(pwd)); cd ..;
cd System; addpath(genpath(pwd)); cd ..; 
cd Localization; addpath(genpath(pwd)); cd ..; 
cd !Simulations;



%% simulation setup
Range_thereshold = 17.5:0.02:18.8;   % [dB]
Ptest = [9,10,11];
data_Pout_analy = zeros(2,length(Range_thereshold),length(Ptest));
data_Pout_empir = zeros(2,length(Range_thereshold),length(Ptest));
pcp.func = 'communication';


for j = 1:length(Ptest)

    P = Ptest(j);

    %% simulation for planar array, Ou_Euler = [0 -90 0].';
    % default setup
    sp = default_system_setup();
    cp = default_channel_setup();
    % update setup
    cp.P = P;
    sp.TypeSA = 'planar';
    sp.Pu = [-10 0 0].';
    sp.Ou_Euler = [0 0 0].';
    cp.K = 2;   % # of subcarriers
    cp.G = 1;   % # of transmissions
    sp = update_system_setup(sp);
    cp = update_channel_setup(cp);
    % generate channel
    cp = gen_channel(sp,cp);
    % generate precoders & combiners
    cp = gen_precoder_combiner(sp,cp,pcp);
    % calculate outage probability
    k = floor(cp.K/2);
    parfor i = 1:length(Range_thereshold)
        disp(['P =', num2str(Ptest(j)), '  **Planar**  ', 'thereshold = ', num2str(Range_thereshold(i)) ,'[dB]   start...']);
        threshold = 10^(Range_thereshold(i)/10);
        method = "analytical";
        Pout = get_outage_probability(sp,cp,threshold,k,method)
        data_Pout_analy(1,i,j) = prod(Pout);
        method = "empirical";
        Pout = get_outage_probability(sp,cp,threshold,k,method);
        data_Pout_empir(1,i,j) = prod(Pout);    
    end
    
    
    
    %% simulation for cuboidal array, Ou_Euler = [0 -90 0].'; 
    % default setup
    sp = default_system_setup();
    cp = default_channel_setup();
    % update setup
    cp.P = P;
    sp.TypeSA = 'cuboidal';
    sp.Pu = [-10 0 0].';
    sp.Ou_Euler = [0 0 0].';
    cp.K = 2;   % # of subcarriers
    cp.G = 1;   % # of transmissions
    sp = update_system_setup(sp);
    cp = update_channel_setup(cp);
    % generate channel
    cp = gen_channel(sp,cp);
    % generate precoders & combiners
    cp = gen_precoder_combiner(sp,cp,pcp);
    % calculate outage probability
    k = floor(cp.K/2);
    parfor i = 1:length(Range_thereshold)
        disp(['P =', num2str(Ptest(j)), '  **Cuboidal**  ', 'thereshold = ', num2str(Range_thereshold(i)) ,'[dB]   start...']);
        threshold = 10^(Range_thereshold(i)/10);
        method = "analytical";
        Pout = get_outage_probability(sp,cp,threshold,k,method);
        data_Pout_analy(2,i,j) = prod(Pout);
        method = "empirical";
        Pout = get_outage_probability(sp,cp,threshold,k,method);
        data_Pout_empir(2,i,j) = prod(Pout);
    end


end

figure(1)
a = 1;
semilogy(Range_thereshold,data_Pout_analy(1,:,a),'r-'); hold on
scatter(Range_thereshold,data_Pout_empir(1,:,a)); hold on
semilogy(Range_thereshold,data_Pout_analy(2,:,a),'b-'); hold on
scatter(Range_thereshold,data_Pout_empir(2,:,a)); hold on
a = 2;
semilogy(Range_thereshold,data_Pout_analy(1,:,a),'r-'); hold on
scatter(Range_thereshold,data_Pout_empir(1,:,a)); hold on
semilogy(Range_thereshold,data_Pout_analy(2,:,a),'b-'); hold on
scatter(Range_thereshold,data_Pout_empir(2,:,a)); hold on
a = 3;
semilogy(Range_thereshold,data_Pout_analy(1,:,a),'r-'); hold on
scatter(Range_thereshold,data_Pout_empir(1,:,a)); hold on
semilogy(Range_thereshold,data_Pout_analy(2,:,a),'b-'); hold on
scatter(Range_thereshold,data_Pout_empir(2,:,a)); hold off
ylim([1e-5, 1]);


% save('Data_Pout_vs_threshold.mat','Range_thereshold','Ptest','data_Pout_analy','data_Pout_empir')
