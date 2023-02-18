function plotGA(v_rel,delta,planetID)
    mu = getPlanetMu(planetID);

    rSOI = getPlanetSOI(planetID);

    v_inf = norm(v_rel);
    specE = v_inf^2/2;
    a = -mu/(2*specE);
    e = 1/sin(delta/2);
    p = a*(1-e^2);

    TA = acos((1 - p/rSOI)/e);
end

