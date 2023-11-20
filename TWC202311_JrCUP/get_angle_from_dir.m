% convert direction to angle
function [azimuth, elevation] = get_angle_from_dir(tv)
    tv = tv ./ norm(tv,2);
    azimuth = atan2d(tv(2,:), tv(1,:));
    elevation = asind(tv(3,:));
end