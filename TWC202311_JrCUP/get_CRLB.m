function [CRLB_eta, CRLB_r] = get_CRLB(sp, scale, CovMat, o3_clkknown)

%% get parameters
K = sp.K;
G = sp.G;

%% get FIM of channel parameters
FIM_etach = zeros(12,12);
for k=1:K
    for g=1:G
        D_eta = CE_get_JM_gk(k,g,sp);
        D_eta2 = [ real(D_eta); imag(D_eta) ];
        Sigma = CovMat(:,:,k,g)*scale^2;
        FIM_etach = FIM_etach + D_eta2.'* Sigma^(-1) *D_eta2;
    end
end

%% compute the CRLB of channel parameters
EFIM_eta = FIM_etach(1:8,1:8) - FIM_etach(1:8,9:12)*FIM_etach(9:12,9:12)^(-1)*FIM_etach(9:12,1:8);
B = EFIM_eta^(-1);
b1 = (180/pi)*sqrt(trace(B(1:2,1:2)));      % in units of degree
b2 = (180/pi)*sqrt(trace(B(3:4,3:4)));
b3 = sqrt(trace(B(5,5)));                   % in units of meter
b4 = sqrt(trace(B(6,6)));
b5 = sqrt(trace(B(7:8,7:8)));               % no unit
CRLB_eta = [b1,b2,b3,b4,b5];

%% transfer FIM to localization parameters
T = Loc_get_JM(sp);
if exist('o3_clkknown','var') 
    if o3_clkknown==1
        T(:,7) = [];
    elseif o3_clkknown==2
        T(:,8) = [];
    end
end
FIM_r = T.'*EFIM_eta*T;
C = FIM_r^(-1);
c1 = sqrt(trace(C(1:3,1:3)));   % in units of meter
c2 = sqrt(trace(C(4:6,4:6)));
c3 = 180/pi*sqrt(C(7,7));       % in units of degree
c4 = sqrt(C(end,end));          % in units of meter
CRLB_r = [c1,c2,c3,c4];


end












