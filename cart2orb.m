function [a,eccx,eccy,inc,raan,theta] = cart2orb(x,y,z,vx,vy,vz,mu)
    % orbit around earth

    % these calculations are done with the provided equations in lecture

    r = [x,y,z];
    v = [vx,vy,vz];

    h = cross(r,v);

    ecc_vec = cross(v,h)/mu - r/norm(r);
    ecc = norm(ecc_vec);

    if ecc ~= 0
        uecc = ecc_vec/ecc;
    else
        uecc = [1 0 0];
    end

    a = (norm(h)^2/mu)/(1-ecc^2);

    inc = rad2deg(acos(h(3)/norm(h)));

    uz = [0,0,1];

    uh = h/norm(h);
    n = cross(uz,uh);
    if norm(n) ~= 0
        un = n/norm(n);
    else
        un = [1 0 0];
    end

    raan = rad2deg(atan2(uh(1),-uh(2)));

    cosw = uecc(1)*un(1) + uecc(2)*un(2);
    sinw = dot(cross(un,uecc),uh);
    argp = atan2(sinw,cosw);

    cosf = dot(uecc,r)/norm(r);
    sinf = dot(cross(uecc,r),uh)/norm(r);
    f = atan2(sinf,cosf);

    eccx = ecc*cos(argp);
    eccy = ecc*sin(argp);

    theta = mod(rad2deg(f + argp),360);

    if theta < 0
        theta = theta + 360;
    end
end