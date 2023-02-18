function TOF = getHohmannTOF(planet1,planet2)
    mu = 1.32712440018E11;

    n1 = 2*pi/getPlanetT(planet1);
    n2 = 2*pi/getPlanetT(planet2);

    a1 = (mu/n1^2)^(1/3);
    a2 = (mu/n2^2)^(1/3);

    a = (a1+a2)/2;

    TOF = (pi/sqrt(mu/a^3))/(3600*24);
end