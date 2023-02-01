function TA = transfer_angle(r1v,r2v)
    % calculated the transfer angle for a space triangle, or simply the
    % angle between two vectors.

    % assuming prograde orbit

    if asind((r1v(1)*r2v(2) - r1v(2)*r2v(1))/(norm(r1v)*norm(r2v))) < 1
        TA = 360-acosd(dot(r1v,r2v)/(norm(r1v)*norm(r2v)));
    else
        TA = acosd(dot(r1v,r2v)/(norm(r1v)*norm(r2v)));
    end
end

