function printOrbitInfo(input,launchDateRange,sequence)
    mu = 1.32712440018E11;
    
    K_launchDate = input(1); %0 to 1
    K_TOFleg1 = input(2); %0 to 1
    K_TOFleg2 = input(3); %0 to 1
    K_TOFleg3 = input(4); %0 to 1
    
    launchDate = launchDateRange(1) +...
        K_launchDate*(launchDateRange(2)-launchDateRange(1)); %days
    TOFleg1 = K_TOFleg1*getHohmannTOF(sequence(1),sequence(2)); %days
    TOFleg2 = K_TOFleg2*getHohmannTOF(sequence(2),sequence(3)); %days
    TOFleg3 = K_TOFleg3*getHohmannTOF(sequence(3),sequence(4)); %days

    %% Calculating Trajectory
    [r1, vp1] = extractEphem(launchDate,sequence(1),true);
    [r2, vp2] = extractEphem(launchDate+TOFleg1,sequence(2),true);
    [r3, vp3] = extractEphem(launchDate+TOFleg1+TOFleg2,sequence(3),true);
    [r4, vp4] = extractEphem(launchDate+TOFleg1+TOFleg2+TOFleg3,sequence(4),true);
    
    [v1,v2i,~,~] = lambert(r1,r2,TOFleg1,0,mu);
    [v2f,v3i,~,~] = lambert(r2,r3,TOFleg2,0,mu);
    [v3f,v4,~,~] = lambert(r3,r4,TOFleg3,0,mu);

    %% Calculating GA information
    dv_reqGA1 = v2f - v2i;
    v_relGA1 = v2i - vp2;
    alphaGA1 = pi - angleBetween(v_relGA1,dv_reqGA1);
    deltaGA1 = pi - 2*alphaGA1;
    
    [dvGA1, passbyr1] = getDvGA(v_relGA1,deltaGA1,sequence(2));
    
    dv_reqGA2 = v3f - v3i;
    v_relGA2 = v3i - vp3;
    alphaGA2 = pi - angleBetween(v_relGA2,dv_reqGA2);
    deltaGA2 = pi - 2*alphaGA2;
    
    [dvGA2, passbyr2] = getDvGA(v_relGA2,deltaGA2,sequence(3));

    earth_vinf = norm(v1 - vp1);
    uranus_vinf = norm(v4 - vp4);

    [dv_tot, dv_earth, dv_uranus] = missiondV(earth_vinf,uranus_vinf);

    fprintf("\nLaunch Date:")
    disp(datetime(launchDate,'convertfrom','juliandate'))
    fprintf("V infinity from Earth: %.4f km/s\n",earth_vinf)
    fprintf("V infinity at Uranus: %.4f km/s\n",uranus_vinf)
    fprintf("Delta V at Earth: %.4f km/s\n",dv_earth)
    fprintf("Delta V at Uranus: %.4f km/s\n",dv_uranus)
    fprintf("Total Mission Delta V: %.4f km/s\n",dv_tot)
    fprintf("\nTotal TOF: %.4f years\n", day2yr(TOFleg1+TOFleg2+TOFleg3))
    if TOFleg1 > yr2day(1)
        fprintf("Leg 1 TOF: %.4f years (K-%.4f)\n", day2yr(TOFleg1),K_TOFleg1)
    else
        fprintf("Leg 1 TOF: %.4f days (K-%.4f)\n", TOFleg1,K_TOFleg1)
    end
    if TOFleg2 > yr2day(1)
        fprintf("Leg 2 TOF: %.4f years (K-%.4f)\n", day2yr(TOFleg2),K_TOFleg2)
    else
        fprintf("Leg 2 TOF: %.4f days (K-%.4f)\n", TOFleg2,K_TOFleg2)
    end
    if TOFleg3 > yr2day(1)
        fprintf("Leg 3 TOF: %.4f years (K-%.4f)\n", day2yr(TOFleg3),K_TOFleg3)
    else
        fprintf("Leg 3 TOF: %.4f days (K-%.4f)\n", TOFleg3,K_TOFleg3)
    end
    fprintf("\nPassby altitude at GA1: %.4f km (%.2f R)\n",...
        passbyr1-getPlanetRad(sequence(2)),passbyr1/getPlanetRad(sequence(2)))
    fprintf("Vinf at GA1: %.4f km/s\n",norm(v_relGA1))
    fprintf("Passby angle at GA1: %.4f degrees\n",rad2deg(deltaGA1))
    fprintf("Magnitude of Delta V from GA1: %.4f km/s\n",dvGA1)
    fprintf("Passby altitude at GA2: %.4f km (%.2f R)\n",...
        passbyr2-getPlanetRad(sequence(3)),passbyr2/getPlanetRad(sequence(3)))
    fprintf("Vinf at GA2: %.4f km/s\n",norm(v_relGA2))
    fprintf("Passby angle at GA1: %.4f degrees\n",rad2deg(deltaGA2))
    fprintf("Magnitude of Delta V from GA2: %.4f km/s\n",dvGA2)
end

function [dVtot,dVe,dVu] = missiondV(vInfe,vInfu)

    %% CONSTANTS
    muE = getPlanetMu(3);
    muU = getPlanetMu(7);
    
    %% ORBIT PARAMETERS
    rpE = 6578;
    raE = 42164;
    
    rpU = 33500;
    raU = 1535860;
    
    %% CALCULATIONS
    aEe = 0.5*(rpE+raE);
    aHe = -(muE)/vInfe^2;
    vpHe = sqrt(muE*((2/rpE)-(1/aHe)));
    vpEe = sqrt(muE*((2/rpE)-(1/aEe)));
    dVe = vpHe - vpEe;
    
    aEu = 0.5*(rpU+raU);
    aHu = -(muU)/vInfu^2;
    vpHu = sqrt(muU*((2/rpU)-(1/aHu)));
    vpEu = sqrt(muU*((2/rpU)-(1/aEu)));
    dVu = vpHu - vpEu;
    
    dVtot = dVe+dVu;
end