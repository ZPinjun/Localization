% This is a function from TeraMIMO:
% S. Tarboush, H. Sarieddeen, H. Chen, M.-H. Loukil, H. Jemaa, M.-S. Alouini, and T. Y. Al-Naffouri, 
%  "TeraMIMO: A Channel Simulator for Wideband Ultra-Massive MIMO Terahertz Communications,"
%  in IEEE Transactions on Vehicular Technology, vol. 70, no. 12, pp. 12325-12341, Dec. 2021, doi: 10.1109/TVT.2021.3123131.


function K_abs = Abs_Coef_Hitran(p)
% =========================================================================
% -- Function to compute exact molecular absorption coefficients in the THz band
%     using the "HITRAN" database, valid in the frequency range: [0.1-10] THz
% =========================================================================

% -- Function: K_abs = Abs_Coef_Hitran(p)

% -- Input Arguments:
%       p: Channel struct that contains the channel parameters

% -- Output Arguments:
%       K_abs: Molecular absorption coefficient, a matrix of size(p.Nsub_b, p.Nsub_c),
%              (number of subbands @ each subcarrier, number of subcarriers)

%=================================================

% -- (c) 2021 Simon Tarboush, Hadi Sarieddeen, Hui Chen, 
%             Mohamed Habib Loukil, Hakim Jemaa, 
%             Mohamed-Slim Alouini, Tareq Y. Al-Naffouri

% -- e-mail: simon.w.tarboush@gmail.com; hadi.sarieddeen@kaust.edu.sa; hui.chen@kaust.edu.sa;
%            mohamedhabib.loukil@kaust.edu.sa; hakim.jemaa@kaust.edu.sa;
%            slim.alouini@kaust.edu.sa; tareq.alnaffouri@kaust.edu.sa

% =========================================================================

% -- References: 
%       Ref [1]: J. M. Jornet and I. F. Akyildiz, "Channel Modeling and Capacity Analysis for Electromagnetic Wireless Nanonetworks in the Terahertz Band,"
%                in IEEE Transactions on Wireless Communications, vol. 10, no. 10, pp. 3211-3221, October 2011, doi: 10.1109/TWC.2011.081011.100545.

% =========================================================================

%% 
K_abs = zeros(p.K,1);


for indx_molc= 1:length(p.molecules)
    % Absorption coefficient for a specific gas mixture

    % load data for each molecule of the transmissionmedium, with a specific mixing ratio,
    % from HITRAN database.
    % For more Molecules, check "Molecular_Absorption" folder --->
    % "Data" subfolders, and edit the inputs from function
    % "generate_channel_param_TIV" or "generate_channel_param_TV"
    Gas_Data = load_Gas_Data(p.molecules{indx_molc}, p.moleculesRatio(indx_molc), p.c);

    % Compute absorption coefficient for a specific gas
    K_Gas = Gas_Abs_Coef(p, Gas_Data);

    % Accumulate all gases to get the total absorption coefficient
    K_abs = K_abs + K_Gas;

end


end