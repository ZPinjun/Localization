% convert angle to direction
function tv = get_dir_from_angle(phi)
    % phi = [az, el] angles in [deg]
    az = phi(1);
    el = phi(2);
    tv = [cosd(az)*cosd(el); sind(az)*cosd(el); sind(el)];
end
