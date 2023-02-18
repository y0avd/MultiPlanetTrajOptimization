function [r,v] = extractEphem(epoch,body,vCalcEnable)
% extractEphem - returns the state vectors of a planet at a given date
% This function uses the planetary ephemeris model described by
%
%      Howard D. Curtis, Orbital Mechanics for Engineering Students 
%      (Third Edition), Butterworth-Heinemann, 2014. (Page 388)
%
% Given an epoch date and a body of interest, this function returns
% the position and velocity vectors of the body in inertial coordinates.
%
% Syntax:  [r,v] = extractEphem(epoch,body)
%
% Inputs:
%    epoch - Julian Date of interest [1x1]
%    body - Integer representing the body of interest [1x1]
%           1 - Mercury
%           2 - Venus
%           3 - Earth
%           4 - Mars
%           5 - Jupiter
%           6 - Saturn
%           7 - Uranus
%           8 - Neptune
%           9 - Pluto
%
% Outputs:
%    r - position vector of planet in inertial coordinates (km) [1x3]
%    v - velocity vector of planet in inertial coordinates (km/s) [1x3]
%
% Example: 
%    [r,v] = extractEphem(2459945.5,3)
%       Returns the position and velocity vectors of Earth at 00:00:00 on
%       Jan 01, 2023.
%
% Other m-files required: kepler_E.m
% Subfunctions: detVelocity, calcVecs, rotMat
% MAT-files required: planetaryEphem.mat
%
% See also: twoGA_trajectory.m
% Author: Andrew Brandt
% Purdue University College of Engineering, West Lafayette, IN
% Email address: brandt32@purdue.edu
% February 2023; Last revision: 13-February-2023

%% Main Function

% INITIALIZATIONS
mu = 1.32712440018E11; % Gravitational Parameter of the sun (km^3/s^2)
AU = 149597870.7; % 1 AU in km
planetaryEphem = load('planetaryEphem.mat').planetaryEphem;
ephem = zeros(1,6);
elem = zeros(1,5);

% Extract Ephemeris of body of Orbit at epoch
Tnot = (epoch - 2451545) / 36525; % Centuries between epoch and J2000
ephemTemp = planetaryEphem((body*2-1):(body*2),:); % Ephemeris of body

% Solve for ephemeris at epoch
%                             semi-major axis given in AU, converts to km
ephem(1) = (ephemTemp(1,1) + ephemTemp(2,1)*Tnot)*AU; % semi-major axis (km)
ephem(2) = (ephemTemp(1,2) + ephemTemp(2,2)*Tnot); % eccentricity

for k = 3:6
    % inclination (deg)
    % RAAN (deg)
    % longitude of perihelion (deg)
    % mean longitude (deg)
%                           Converts arc-seconds to degrees
    ephem(k) = ephemTemp(1,k) + (ephemTemp(2,k)/3600)*Tnot;
    % Keeps angle values between 0 and 360 degrees
    if ephem(k) < 0 
        ephem(k) = mod(ephem(k),360);
    elseif ephem(1,k) > 360
        ephem(k) = mod(ephem(k),360);
    end
end

% Solve for other useful orbital elements
elem(1) = (mu*ephem(1)*(1-ephem(2)^2))^(1/2); %specific angular momentum (km^2/s)
elem(2) = deg2rad(ephem(5) - ephem(4)); % argument of perihelion (rad)
elem(3) = deg2rad(ephem(6) - ephem(5)); % mean anomaly (rad)
elem(4) = kepler_E(ephem(2),elem(3)); % Eccentric Anomaly (rad)
elem(5) = 2*atan(sqrt((1+ephem(2))/(1-ephem(2)))*tan(elem(4)/2)); % True anomaly (rad)

% Determine position vector
r = detPosition(ephem(1),ephem(2),elem(5),elem(2),deg2rad(ephem(3)),deg2rad(ephem(4)));

% Determine velocity vector
if vCalcEnable % If toggle is enabled calculate velocity
    v = detVelocity(r,ephem(1),ephem(2),elem(5),elem(2),deg2rad(ephem(3)),deg2rad(ephem(4)),mu);
else % if toggle is off computation is saved by not calculating velocity
    v = [0,0,0];
end

r = r';
v = v';

end

%% Subfunctions

function v = detVelocity(R,a,e,ta,w,i,RAAN,mu)
% detVelocity - Determines the velocity vector of a body in orbit.
% Given basic Keplerian orbital elements of a body in orbit, the function
% returns the velocity vector of the body in inertial coordinates.
%
% Syntax:  v = detVelocity(R,a,e,ta,w,i,RAAN,mu)
%
% Inputs:
%    R - position vector of the body in inertial coordinates (km) [3x1]
%    a - semi-major axis of the orbit (km) [1x1]
%    e - eccentricity of the orbit (dimensionless) [1x1]
%    ta - true anomaly of the body in orbit (rad) [1x1]
%    w - argument of periapsis of the orbit (rad) [1x1]
%    i - inclination of the orbit (rad) [1x1]
%    RAAN - right ascenion of the ascending node of the orbit (rad) [1x1]
%    mu - gravitational parameter of the central body (km^3/s^2) [1x1]
%
% Outputs:
%    v - velocity vector of the body in inertial coordiantes (km/s) [3x1]
%
% theta value for 3D rotation matrix
theta = w + ta;
% Magnitudes of r and v
r = norm(R);
v = sqrt(mu*((2/r)-(1/a)));
% Orbital elements semi-latus rectum and sp. angular mom.
p = a*(1-e^2);
h = sqrt(p*mu);
% Flight path angle, solve quadrant ambiguity based on position in orbit
if ta >= 0 && ta <= pi
    fpa = abs(acos(h/(v*r)));
else
    fpa = -abs(acos(h/(v*r)));
end
% Velocity vector in orbit fixed frame
v_orbit = [v*sin(fpa),v*cos(fpa),0];
% 3-1-3 body sequence from orbit fixed frame to inertial
rotSeq = rotMat(RAAN,3)*rotMat(i,1)*rotMat(theta,3);
% Vector result
v = rotSeq * v_orbit';
end

function r = detPosition(a,e,ta,w,i,RAAN)
% detPosition - Determines the position of a body in orbit.
% Given basic Keplerian orbital elements of a body in orbit, the function
% returns the position vector of the body in inertial coordinates.
%
% Syntax:  r = detPosition(a,e,ta,w,i,RAAN)
%
% Inputs:
%    a - semi-major axis of the orbit (km) [1x1]
%    e - eccentricity of the orbit (dimensionless) [1x1]
%    ta - true anomaly of the body in orbit (rad) [1x1]
%    w - argument of periapsis of the orbit (rad) [1x1]
%    i - inclination of the orbit (rad) [1x1]
%    RAAN - right ascenion of the ascending node of the orbit (rad) [1x1]
%
% Outputs:
%    r - position vector of the body in inertial coordiantes (km/s) [3x1]
%
% theta value for 3D rotation matrix
theta = w + ta;
% Semi-latus rectum
p = a*(1-e^2);
% Polar equation for conic sections, gives magnitude of r
r = p/(1 + (e*cos(ta)));
% Position vector in the orbit fixed frame
posOrbitFrame = [r;0;0];
% 3-1-3 body sequence from orbit fixed frame to inertial
rotSeq = rotMat(RAAN,3)*rotMat(i,1)*rotMat(theta,3);
% Vector result
r = rotSeq * posOrbitFrame;
end

function R = rotMat(angle,axis)
% rotMat - Outputs a 3x3 rotation matrix given a body axis and angle.
% This function outputs a 3x3 rotation matrix given an angle and an axis of
% rotation. Intended to be used to generate a DCM for an Euler Sequence
% about body axes.
%
% Syntax:  R = rotMat(angle,axis)
%
% Inputs:
%    angle - angle of rotation (rad) [1x1]
%    axis - integer of the body axis of rotation (1, 2, or 3) [1x1]
%
% Outputs:
%    R - resulting rotation matrix of the given angle and axis [3x3]
%
switch axis
    case 1
        R = [1,0,0;
             0,cos(angle),-sin(angle);
             0,sin(angle),cos(angle)];
    case 2
        R = [cos(angle),0,sin(angle);
             0,1,0;
             -sin(angle),0,cos(angle)];
    case 3
        R = [cos(angle),-sin(angle),0;
             sin(angle),cos(angle),0;
             0,0,1];
end
end