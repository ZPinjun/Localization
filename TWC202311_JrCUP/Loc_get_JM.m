function D = Loc_get_JM(sp)

%% get parameters
pb = sp.Pb;
pr = sp.Pr;
pu = sp.Pu;
u1 = [1 0 0].';
u2 = [0 1 0].';
u3 = [0 0 1].';
Rb = sp.Rb;
o3 = sp.Or_Euler(1);
o2 = sp.Or_Euler(2);
o1 = sp.Or_Euler(3);
Rr = eul2rotm(reshape(deg2rad(sp.Or_Euler), [1,3]), 'ZYX');

%% compute Jacobian matrix
tub = pu - pb;
dub = norm(tub,2);
tub_b = Rb.'*tub;
trb = pr - pb;
drb = norm(trb,2);
trb_b = Rb.'*trb;
tur = pu - pr;
dur = norm(tur,2);


% row 1
term1 = ( 1 + (tub_b(2)/tub_b(1))^2 )^(-1);
dazL_dpu = term1 * ( tub_b(1)*Rb*u2 - tub_b(2)*Rb*u1 ) / tub_b(1)^2;
row1 = [dazL_dpu.', zeros(1,5)];

% row 2
term1 = ( 1 - (tub_b(3)/dub)^2 )^(-1/2);
delL_dpu = term1 * ( Rb*u3/dub - tub_b(3)*tub/dub^3 );
row2 = [delL_dpu.', zeros(1,5)];

% row 3
term1 = ( 1 + (trb_b(2)/trb_b(1))^2 )^(-1);
dazR_dpr = term1 * ( trb_b(1)*Rb*u2 - trb_b(2)*Rb*u1 ) / trb_b(1)^2;
row3 = [zeros(1,3), dazR_dpr.', zeros(1,2)];

% row 4
term1 = ( 1 - (trb_b(3)/drb)^2 )^(-1/2);
delR_dpr = term1 * ( Rb*u3/drb - trb_b(3)*trb/drb^3 );
row4 = [zeros(1,3), delR_dpr.', zeros(1,2)];

% row 5
dtauL_dpu = tub/dub;
row5 = [dtauL_dpu.', zeros(1,4), 1];

% row 6
dtauR_dpu = tur/dur;
dtauR_dpr = -tur/dur + trb/drb;
row6 = [dtauR_dpu.', dtauR_dpr.', 0, 1];

% row 7
dvar2_dpu = Rr*u2/dur - (u2.'*Rr.'*tur)*tur/dur^3;
dvar2_dpr = -Rr*u2/dur + (u2.'*Rr.'*tur)*tur/dur^3 - Rr*u2/drb + (u2.'*Rr.'*trb)*trb/drb^3;
v1 = [-sind(o2)*sind(o1)*sind(o3)-cosd(o3)*cosd(o1), sind(o2)*sind(o1)*cosd(o3)+sind(o3)*cosd(o1), 0];
dvar2_do3 = v1 * ( tur/dur - trb/drb );
row7 = [dvar2_dpu.', dvar2_dpr.', dvar2_do3, 0];

% row 8
dvar3_dpu = Rr*u3/dur - (u3.'*Rr.'*tur)*tur/dur^3;
dvar3_dpr = -Rr*u3/dur + (u3.'*Rr.'*tur)*tur/dur^3 - Rr*u3/drb + (u3.'*Rr.'*trb)*trb/drb^3;
v2 =  [-sind(o2)*cosd(o1)*sind(o3)+cosd(o3)*sind(o1), sind(o2)*cosd(o1)*cosd(o3)+sind(o3)*sind(o1), 0];
dvar3_do3 = v2 * ( tur/dur - trb/drb );
row8 = [dvar3_dpu.', dvar3_dpr.', dvar3_do3, 0];


D = [row1; row2; row3; row4; row5; row6; row7; row8];


end

