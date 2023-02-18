function [dvGA, rp] = getDvGA(v_rel,delta,planetID)
    mu = getPlanetMu(planetID);

    v_inf = norm(v_rel);
    specE = v_inf^2/2;
    a = -mu/(2*specE);
    rp = a*(1-1/sin(delta/2));

    dvGA = sqrt(2*v_inf^2 - 2*v_inf^2*cos(delta));
end

