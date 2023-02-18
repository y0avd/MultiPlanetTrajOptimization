function mu = getPlanetMu(planetID)
    switch planetID
        case 1
            mu = 2.2032e4;
        case 2
            mu = 3.24859e5;
        case 3
            mu = 3.986004418e5;
        case 4
            mu = 4.282837e4;
        case 5
            mu = 1.26686534e8;
        case 6
            mu = 3.7931187e7;
        case 7
            mu = 5.793939e6;
    end
end

