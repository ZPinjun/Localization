% convert direction to angle

function output = get_angle_from_dir(tv)
    %tv = normalize(tv, 'norm', 2);
    tv = tv/norm(tv,2);
    azimuth = atan2d(tv(2,:), tv(1,:));
    elevation = asind(tv(3,:));
    output = [azimuth; elevation];
end