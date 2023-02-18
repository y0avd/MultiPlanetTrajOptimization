% K_launchDate = input(1) %0 to 1
% K_TOFleg1 = input(2); %0 to 1
% K_TOFleg2 = input(3); %0 to 1
% K_TOFleg3 = input(4); %0 to 1

% min_var = 1 (TOF) 2 (dV)

function output = twoGA_Trajectory(input,launchDateRange,sequence,min_var)
%% Initializations
mu = 1.32712440018E11;

K_launchDate = input(1); %0 to 1
K_TOFleg1 = input(2); %0 to 1

launchDate = launchDateRange(1) +...
    K_launchDate*(launchDateRange(2)-launchDateRange(1)); %days
TOFleg1 = K_TOFleg1*getHohmannTOF(sequence(1),sequence(2)); %days

%% Calculating Trajectory
[r1, vp1] = extractEphem(launchDate,sequence(1),true);
r2 = extractEphem(launchDate+TOFleg1,sequence(2),false);

[v1,~,~,~] = lambert(r1,r2,TOFleg1,0,mu);

if min_var == 1
    output = -(TOFleg1 + TOFleg2 + TOFleg3);
elseif min_var == 2
    % Calculate deltaV from LEO
    earth_vinf = norm(v1 - vp1);
    
    output = earth_vinf;
else
    disp('ERROR: Invalid min_var')
end