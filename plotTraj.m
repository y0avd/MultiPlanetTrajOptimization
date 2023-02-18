function plotTraj(input,launchDateRange,sequence)
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
    [r1, ~] = extractEphem(launchDate,sequence(1),true);
    [r2, vp2] = extractEphem(launchDate+TOFleg1,sequence(2),true);
    [r3, vp3] = extractEphem(launchDate+TOFleg1+TOFleg2,sequence(3),true);
    [r4, ~] = extractEphem(launchDate+TOFleg1+TOFleg2+TOFleg3,sequence(4),true);
    
    [~,v2i,~,~] = lambert(r1,r2,TOFleg1,0,mu);
    [v2f,v3i,~,~] = lambert(r2,r3,TOFleg2,0,mu);
    [v3f,~,~,~] = lambert(r3,r4,TOFleg3,0,mu);

    %% Calculating GA information
    dv_reqGA1 = v2f - v2i;
    v_relGA1 = v2i - vp2;
    alphaGA1 = pi - angleBetween(v_relGA1,dv_reqGA1);
    deltaGA1 = pi - 2*alphaGA1;
    
    dv_reqGA2 = v3f - v3i;
    v_relGA2 = v3i - vp3;
    alphaGA2 = pi - angleBetween(v_relGA2,dv_reqGA2);
    deltaGA2 = pi - 2*alphaGA2;

    % Plotting
    figure;
    title("Interplanetary Trajectory")
    hold on; grid on; axis equal
    plotPlanetaryOrbits(launchDate,sequence)

    plotLambertArc(r1,r2,TOFleg1,mu)
    plotLambertArc(r2,r3,TOFleg2,mu)
    plotLambertArc(r3,r4,TOFleg3,mu)

    figure;
    title("GA1 Passby Trajectory")
    hold on; grid on; axis equal
    %plotGA(v_relGA1,deltaGA1,sequence(2))

    figure;
    title("GA2 Passby Trajectory")
    hold on; grid on; axis equal
    %plotGA(v_relGA2,deltaGA2,sequence(3))
end

%% SUBFUNCTIONS
function plotVector(K,r,v,color)
    plot3([r(1),r(1) + K*v(1)],...
        [r(2),r(2) + K*v(2)],...
        [r(3),r(3) + K*v(3)],color)
end

function plotPlanetaryOrbits(launchDate, sequence)
    nump = 1e3;

    r_sun = 696000;
    [X,Y,Z] = sphere;
    surf(X*r_sun,Y*r_sun,Z*r_sun)

    for i = 1:length(sequence)
        T = getPlanetT(sequence(i))/(3600*24);

        tspan = linspace(0,T,nump);

        for k = 1:nump
            r = extractEphem(launchDate+tspan(k),sequence(i),false);
            x(k) = r(1); %#ok<AGROW> 
            y(k) = r(2); %#ok<AGROW> 
            z(k) = r(3); %#ok<AGROW> 
        end

        plot3(x,y,z,'k')
    end
end

function plotLambertArc(r1,r2,TOF,mu)
    nump = 1e3;

    [v1,v2] = lambert(r1,r2,TOF,0,mu);

    [a,eccx,eccy,inc,raan,theta1] = cart2orb(...
        r1(1),r1(2),r1(3),v1(1),v1(2),v1(3),mu);
    [~,~,~,~,~,theta2] = cart2orb(...
        r2(1),r2(2),r2(3),v2(1),v2(2),v2(3),mu);

    theta = linspace(theta1,theta2,nump);

    if theta2 > theta1
        theta = linspace(theta1,theta2,nump);
    else
        theta = linspace(theta1,360 + theta2,nump);
    end

    for i = 1:nump
        [x(i),y(i),z(i)] = orb2cart(a,eccx,eccy,inc,raan,theta(i),mu);
    end

    plot3(x,y,z);
end

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