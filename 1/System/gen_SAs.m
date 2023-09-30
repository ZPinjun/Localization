function [Dsubarray,Rsubarray] = gen_SAs(TypeSubarray)

r = 0.05; % spatial distance of the SAs

if TypeSubarray == "planar"
    % position of SAs (in UE's coordinate system)
    Dsubarray = [   
        0 -r -r;
        0 -r r;
        0 0 -r;
        0 0 r;
        0 r -r;
        0 r r;
        ].';
    % orientation of SAs (in UE's coordinate system)
    R_SA = zeros(3,3,6);
    R_SA(:,:,1) = eye(3);   
    R_SA(:,:,2) = eye(3);
    R_SA(:,:,3) = eye(3);
    R_SA(:,:,4) = eye(3);
    R_SA(:,:,5) = eye(3);
    R_SA(:,:,6) = eye(3);
    Rsubarray = R_SA;
elseif TypeSubarray == "cuboidal"
    % position of SAs (in UE's coordinate system)
    Dsubarray = [   
        r 0 0;
        -r 0 0;
        0 -r 0;
        0 r 0;
        0 0 -r;
        0 0 r
        ].';
    % orientation of SAs (in UE's coordinate system)
    R_SA = zeros(3,3,6);
    R_SA(:,:,1) = eye(3);   % front
    R_SA(:,:,2) = eul2rotm([0 pi 0], 'ZYX');   % back
    R_SA(:,:,3) = eul2rotm([-0.5*pi 0 0], 'ZYX');   % left
    R_SA(:,:,4) = eul2rotm([0.5*pi 0 0], 'ZYX');   % right
    R_SA(:,:,5) = eul2rotm([0 0.5*pi 0], 'ZYX');   % bottom
    R_SA(:,:,6) = eul2rotm([0 -0.5*pi 0], 'ZYX');   % top
    Rsubarray = R_SA;
end


end

