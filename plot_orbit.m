function plot_orbit(orb,theta_range,AOP,color)
    theta = linspace(theta_range(1),theta_range(2),1e4);
    orb.p = orb.a*(1-orb.e^2);
    r = orb.p./(1 + orb.e*cosd(theta));

    au = 149597870;
    
    rx = r.*cosd(theta+AOP)/au;
    ry = r.*sind(theta+AOP)/au;
    
    plot(rx,ry,color)
end