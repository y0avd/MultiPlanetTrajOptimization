function [x,y,z,vx,vy,vz] = orb2cart(a,eccx,eccy,inc,raan,theta,mu)
    % orbit around earth

    % this calculation is done usign simple Keplerian math and a DCM on the
    % basic 3-1-3 Euler Angle rotation to turn distance and velocity into
    % x y z coordinates

    inc = deg2rad(inc);
    raan = deg2rad(raan);
    theta = deg2rad(theta);

    ecc = norm([eccx; eccy]);

    argp = atan2(eccy,eccx);

    f = theta - argp;

    p = a*(1 - ecc^2);

    DCM = [cos(raan)*cos(argp) - sin(raan)*sin(argp)*cos(inc),...
        -cos(raan)*sin(argp) - sin(raan)*cos(argp)*cos(inc),...
        sin(raan)*sin(inc);...
        sin(raan)*cos(argp) + cos(raan)*sin(argp)*cos(inc), ...
        -sin(raan)*sin(argp) + cos(raan)*cos(argp)*cos(inc),...
        -cos(raan)*sin(inc);...
        sin(argp)*sin(inc), cos(argp)*sin(inc), cos(inc)];

    rv = [cos(f);sin(f);0];
    rv = p*rv/(1 + ecc*cos(f));
    rv = DCM*rv;

    vv = [-sin(f); (ecc + cos(f)); 0];
    vv = sqrt(mu/p)*vv;
    vv = DCM*vv;

    x = rv(1);
    y = rv(2);
    z = rv(3);
    vx = vv(1);
    vy = vv(2);
    vz = vv(3);
end