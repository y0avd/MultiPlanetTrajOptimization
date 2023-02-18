% Data from
% https://nssdc.gsfc.nasa.gov/planetary/factsheet/

function R = getPlanetRad(planetID)
    switch planetID
        case 1
            R = 4879/2;
        case 2
            R = 12104/2;
        case 3
            R = 12756/2;
        case 4
            R = 6792/2;
        case 5
            R = 142984/2;
        case 6
            R = 120536/2;
        case 7
            R = 51118/2;
    end
end

