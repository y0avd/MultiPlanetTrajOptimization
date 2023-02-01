function TOF = TOF_calc(mu,a,e,E1,E2)
    t1p = (deg2rad(E1) - e*sind(E1))/sqrt(mu/a^3);
    t2p = (deg2rad(E2) - e*sind(E2))/sqrt(mu/a^3);
    TOF = t2p - t1p;
end

