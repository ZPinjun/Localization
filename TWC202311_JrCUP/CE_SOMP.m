function eta_OMP = CE_SOMP(y, sp, dicsize)

%% get parameters
K = sp.K;
G = sp.G;
df = sp.df;
c = sp.c;
s = sp.s;
W = sp.W;
Gamma = sp.Gamma;
N1 = sp.RFC_dim(1); 
N2 = sp.RFC_dim(2);
Nb1 = sp.BS_dim(1);
Nb2 = sp.BS_dim(2);
Nr1 = sp.RIS_dim(1);
Nr2 = sp.RIS_dim(2);
lambdac = sp.lambda;
pAE_BS = sp.pAE_BS;
pAE_RIS = sp.pAE_RIS;


%% step 1: estimate thetaL & thetaR using SOMP
% --- generate dictionary
theta_az_dic = linspace(-90+rand(1),90+rand(1), dicsize);
theta_el_dic = linspace(-90+rand(1),90+rand(1), dicsize);
theta_dic = [kron(theta_az_dic,ones(1,dicsize)); kron(ones(1,dicsize),theta_el_dic)];
A = zeros(Nb1*Nb2,dicsize^2);
for i = 1:dicsize^2
    A(:,i) = get_array_response(lambdac, pAE_BS, theta_dic(:,i));
end
A = W'*A;
% --- call SOMP algorithm
S = reshape(y,[N1*N2,K*G]);
[~,supp] = func_SOMP(A,S,2);
theta1 = theta_dic(:,supp(1));
theta2 = theta_dic(:,supp(2));


%% step 2: estimate tauL & tauR
A_sub = A(:,supp);
Mat = zeros(G,K,2);
for k = 1:K
    for g = 1:G
        xk = (A_sub'*A_sub)^(-1)*A_sub'*y(:,k,g);
        Mat(g,k,1) = xk(1);
        Mat(g,k,2) = xk(2);
    end
end
Mat(:,:,1) = Mat(:,:,1)./s.';
Mat(:,:,2) = Mat(:,:,2)./s.';
% --- identify two paths based on the variance of (column of) Mat1 & Mat2
cov1 = cov(Mat(:,1,1));
cov2 = cov(Mat(:,1,2));
order = [1, 2]; % 1: LOS path; 2: RIS path
thetaL = theta1;
thetaR = theta2;
if cov1 > cov2
    order = [2, 1];
    thetaL = theta2;
    thetaR = theta1;
end
% --- perform 1D search
tauTest = linspace(30,50,0.1*dicsize^2);
vec1 = zeros(length(tauTest),1);
vec2 = zeros(length(tauTest),1);
for i = 1:length(tauTest)
    tau = tauTest(i);
    d = get_frequency_response(tau, K, df, c);
    vec1(i) = norm(Mat(:,:,1)*d,2);
    vec2(i) = norm(Mat(:,:,2)*d,2);
end
[~,ind1] = max(vec1);
[~,ind2] = max(vec2);
ind = [ind1; ind2];
tauL = tauTest(ind(order(1)));
tauR = tauTest(ind(order(2)));


%% step 3: estimate vartheta
M1 = Mat(:,:,order(2));
dR_est = get_frequency_response(tauR, K, df, c);
b1 = M1*dR_est;
% --- generate dictionary
vartheta2_dic = linspace(-2, 2, dicsize);
vartheta3_dic = linspace(-2, 2, dicsize);
vartheta_dic = [zeros(1,dicsize^2); kron(vartheta2_dic,ones(1,dicsize)); kron(ones(1,dicsize),vartheta3_dic)];
A1 = zeros(Nr1*Nr2,dicsize^2);
for i = 1:dicsize^2
    tR_tilde = vartheta_dic(:,i);
    A1(:,i) = exp((1j*2*pi/lambdac)*tR_tilde.'*pAE_RIS).';
end
A1 = Gamma'*A1;
% --- call OMP algorithm
[~,supp] = func_OMP(A1,b1,1);
vartheta = vartheta_dic(2:3,supp);


%% form outputs
eta_OMP = [thetaL(1),thetaL(2),tauL,thetaR(1),thetaR(2),tauR,vartheta(1),vartheta(2)];


%% ================================ functions
    function a = get_array_response(lambdac, pE, phi)
        t = get_dir_from_angle(phi);
        a = exp((1j*2*pi/lambdac)*t.'*pE).';
    end

    function d = get_frequency_response(tau, K, df, c)
        kk = (1:K).';
        d = exp(1j*2*pi*(kk-1)*df*tau/c);
    end

end

