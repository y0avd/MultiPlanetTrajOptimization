clc, clear, close all;

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
global rp_U

rp_E1 = 500 + earth.R; %km
rp_U = 3000 + uranus.R; %km

x0 = [1 15 yr2sec(.45) 1.5 1.2*earth.R 1.2*jupiter.R];

output = EJGAtraj(x0, "dv",1);