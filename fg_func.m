function output = fg_func(r1v,v1v,TA1,r2,TA2,mu,p,output_type)

    dTA = TA2 - TA1;

    f = 1 - (r2/p)*(1 - cosd(dTA));
    g = (r2*norm(r1v)/sqrt(mu*p))*sind(dTA);
    df = (dot(r1v,v1v)/(p*norm(r1v)))*(1 - cosd(dTA)) - (1/norm(r1v))*sqrt(mu/p)*sind(dTA);
    dg = 1 - (norm(r1v)/p)*(1-cosd(dTA));

    switch output_type
        case "rv"
            output = f*r1v + g*v1v;
        case "vv"
            output = df*r1v + dg*v1v;
    end
end