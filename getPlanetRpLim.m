function rp = getPlanetRpLim(planetID)
    switch planetID
        case 1
            alt = 1e3;
            rp = alt + getPlanetRad(planetID);
        case 2
            alt = 1e3;
            rp = alt + getPlanetRad(planetID);
        case 3
            alt = 1e3;
            rp = alt + getPlanetRad(planetID);
        case 4
            alt = 1e3;
            rp = alt + getPlanetRad(planetID);
        case 5
            alt = 5e3;
            rp = alt + getPlanetRad(planetID);
        case 6
            alt = 1e3;
            rp = alt + getPlanetRad(planetID);
        case 7
            alt = 1e3;
            rp = alt + getPlanetRad(planetID);
    end
end

