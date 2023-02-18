function T = getPlanetT(planetID)
    mu = 1.32712440018E11;

    AU = 149597870.7;

    switch planetID
        case 1
            a = 0.38709893*AU;
        case 2
            a = 0.72333199*AU;
        case 3
            a = AU;
        case 4
            a = 1.52366231*AU;
        case 5
            a = 5.20336301*AU;
        case 6
            a = 9.53707032*AU;
        case 7
            a = 19.19126393*AU;
    end

    T = 2*pi*sqrt(a^3/mu);
end

