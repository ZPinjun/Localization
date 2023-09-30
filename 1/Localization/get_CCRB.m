function CCRB = get_CCRB(cp,sp)

% get parameters
visi_mat = sp.visi_mat;
D = sum(sum(visi_mat));

% connot do localization if there is only one BS can be seen
if sum(sum(visi_mat.') > 0) <= 1
    CCRB = nan(2,1);
else
    % get FIM to channel parameters
    FIM_eta_ch = get_FIM_channel_params(cp,sp);
    
    if sum(sum(abs(FIM_eta_ch))) == 0
        CCRB = nan(2,1);
    else
        % get EFIM to localization-related parameters
        A = FIM_eta_ch(1:5*D,1:5*D);
        B = FIM_eta_ch(1:5*D, (5*D+1):7*D);
        C = FIM_eta_ch((5*D+1):7*D,(5*D+1):7*D);
        FIM_eta = A - B*C^(-1)*B.';
        % get FIM to state parameters
        FIM = get_FIM_state_params(FIM_eta,sp);
        % add rotation matrix constraint
        R = sp.Ru;
        r1 = R(:,1);
        r2 = R(:,2);
        r3 = R(:,3);
        M1 = [-r3 zeros(3,1) r2; zeros(3,1) -r3 -r1; r1 r2 zeros(3,1)];
        M1 = M1/sqrt(2);
        M2 = [zeros(3,3), eye(3), zeros(3,1); M1, zeros(size(M1,1),4); zeros(1,6), 1];
        MM = M2*(M2'*FIM*M2)^(-1)*M2';
        CRLB_pos = abs(sqrt(trace(MM(1:3,1:3))));
        CRLB_ori = abs(40.5142*sqrt(trace(MM(4:12,4:12))));
    
        CCRB = [CRLB_pos; CRLB_ori];
    end
end

end

