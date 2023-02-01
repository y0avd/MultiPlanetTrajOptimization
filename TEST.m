clc, clear, close all;

au = 149597870; %km

global sun
global earth
global jupiter
global uranus

sun = struct('mu',1.32712440018e11);

earth = struct('R',6378.1363,'mu',398600.435507,'a',au,'e',0);
earth.vc = sqrt(sun.mu/earth.a);

jupiter = struct('R',69911,'mu',1.26686534e8,'a',5.2038*au,'e',0);
jupiter.vc = sqrt(sun.mu/jupiter.a);

uranus = struct('R',25559,'mu',5794556.4,'a',19.19126*au,'e',0);
uranus.vc = sqrt(sun.mu/uranus.a);

global rp_E1
global rp_E2
global rp_J
global rp_U

rp_E1 = 500 + earth.R; %km
rp_E2 = 1000 + earth.R; %km
rp_J = 30000 + jupiter.R; %km
rp_U = 10000 + uranus.R; %km

x0 = [10 1 yr2sec(0.9) 4.2 0];

output = EJGAtraj(x0, "dv");