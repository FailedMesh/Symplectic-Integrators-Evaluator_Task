clear; clc;

%%
%INITIALIZATION
E = [0; 0; 0];  %Set Electric Field
B = [0; 0; 1];  %Set Magnetic Field

c = 3e+8;       %Speed of Light
m = 9.109e-31;  %Mass of electron
q = -1.60218e-19;   %Charge of electron

h = 1e-8;   %Time Step

%Initial position:
pos(1,1) = 0;
pos(2,1) = 0;
pos(3,1) = 0;

%Initial velocity:
vel(1,1) = 0;
vel(2,1) = 0.9*c;
vel(3,1) = 0;

%Initial acceleration:
acc(1,1) = q*(E(1)/m + (vel(2,1)*B(3) - vel(3,1)*B(2)));
acc(2,1) = q*(E(2)/m + (vel(3,1)*B(1) - vel(1,1)*B(3)));
acc(3,1) = q*(E(3)/m + (vel(1,1)*B(2) - vel(2,1)*B(1)));


%%
%LOOP SETUP

%Define the number of turns to loop the electron
max_turns = 2;

%Initialize loop parameters
num_turns = 0;
i = h;
index = 0;

%%
%LOOP IN DormandPrince 5th Order Method

while num_turns < max_turns
    index = int32((1/h)*(i) + 1);
    
    k0vel1 = h*q*(E(1) + (vel(2,index-1)*B(3) - vel(3,index-1)*B(2)))/m;
    k0vel2 = h*q*(E(2) + (vel(3,index-1)*B(1) - vel(1,index-1)*B(3)))/m;
    k0vel3 = h*q*(E(3) + (vel(1,index-1)*B(2) - vel(2,index-1)*B(1)))/m;
    
    vel1inc = (1/5)*h*k0vel1;
    vel2inc = (1/5)*h*k0vel2;
    vel3inc = (1/5)*h*k0vel3;
    
    k1vel1 = h*q*(E(1) + ((vel(2,index-1) + vel2inc)*B(3) - (vel(3,index-1) + vel3inc)*B(2)))/m;
    k1vel2 = h*q*(E(2) + ((vel(3,index-1) + vel3inc)*B(1) - (vel(1,index-1) + vel1inc)*B(3)))/m;
    k1vel3 = h*q*(E(3) + ((vel(1,index-1) + vel1inc)*B(2) - (vel(2,index-1) + vel2inc)*B(1)))/m;
    
    vel1inc = (3/10)*h*(k0vel1 + 3*k1vel1/4);
    vel2inc = (3/10)*h*(k0vel2 + 3*k1vel2/4);
    vel3inc = (3/10)*h*(k0vel3 + 3*k1vel3/4);
    
    k2vel1 = h*q*(E(1) + ((vel(2,index-1) + vel2inc)*B(3) - (vel(3,index-1) + vel3inc)*B(2)))/m;
    k2vel2 = h*q*(E(2) + ((vel(3,index-1) + vel3inc)*B(1) - (vel(1,index-1) + vel1inc)*B(3)))/m;
    k2vel3 = h*q*(E(3) + ((vel(1,index-1) + vel1inc)*B(2) - (vel(2,index-1) + vel2inc)*B(1)))/m;
    
    vel1inc = (4/5)*h*(11*k0vel1 - 42*k1vel1 + 40*k2vel1)/9;
    vel2inc = (4/5)*h*(11*k0vel2 - 42*k1vel2 + 40*k2vel2)/9;
    vel3inc = (4/5)*h*(11*k0vel3 - 42*k1vel3 + 40*k2vel3)/9;
    
    k3vel1 = h*q*(E(1) + ((vel(2,index-1) + vel2inc)*B(3) - (vel(3,index-1) + vel3inc)*B(2)))/m;
    k3vel2 = h*q*(E(2) + ((vel(3,index-1) + vel3inc)*B(1) - (vel(1,index-1) + vel1inc)*B(3)))/m;
    k3vel3 = h*q*(E(3) + ((vel(1,index-1) + vel1inc)*B(2) - (vel(2,index-1) + vel2inc)*B(1)))/m;
    
    vel1inc = (8/9)*h*(4843*k0vel1 - 19020*k1vel1 + 16112*k2vel1 - 477*k3vel1)/1458;
    vel2inc = (8/9)*h*(4843*k0vel2 - 19020*k1vel2 + 16112*k2vel2 - 477*k3vel2)/1458;
    vel3inc = (8/9)*h*(4843*k0vel3 - 19020*k1vel3 + 16112*k2vel3 - 477*k3vel3)/1458;
    
    k4vel1 = h*q*(E(1) + ((vel(2,index-1) + vel2inc)*B(3) - (vel(3,index-1) + vel3inc)*B(2)))/m;
    k4vel2 = h*q*(E(2) + ((vel(3,index-1) + vel3inc)*B(1) - (vel(1,index-1) + vel1inc)*B(3)))/m;
    k4vel3 = h*q*(E(3) + ((vel(1,index-1) + vel1inc)*B(2) - (vel(2,index-1) + vel2inc)*B(1)))/m;
    
    vel1inc = h*(477901*k0vel1 - 1806240*k1vel1 + 1495424*k2vel1 + 46746*k3vel1 - 45927*k4vel1)/167904;
    vel2inc = h*(477901*k0vel2 - 1806240*k1vel2 + 1495424*k2vel2 + 46746*k3vel2 - 45927*k4vel2)/167904;
    vel3inc = h*(477901*k0vel3 - 1806240*k1vel3 + 1495424*k2vel3 + 46746*k3vel3 - 45927*k4vel3)/167904;
    
    k5vel1 = h*q*(E(1) + ((vel(2,index-1) + vel2inc)*B(3) - (vel(3,index-1) + vel3inc)*B(2)))/m;
    k5vel2 = h*q*(E(2) + ((vel(3,index-1) + vel3inc)*B(1) - (vel(1,index-1) + vel1inc)*B(3)))/m;
    k5vel3 = h*q*(E(3) + ((vel(1,index-1) + vel1inc)*B(2) - (vel(2,index-1) + vel2inc)*B(1)))/m;
    
    vel1inc = h*(12985*k0vel1 + 64000*k2vel1 + 92750*k3vel1 - 45927*k4vel1 + 18656*k5vel1)/142464;
    vel2inc = h*(12985*k0vel2 + 64000*k2vel2 + 92750*k3vel2 - 45927*k4vel2 + 18656*k5vel2)/142464;
    vel3inc = h*(12985*k0vel3 + 64000*k2vel3 + 92750*k3vel3 - 45927*k4vel3 + 18656*k5vel3)/142464;
    
    k6vel1 = h*q*(E(1) + ((vel(2,index-1) + vel2inc)*B(3) - (vel(3,index-1) + vel3inc)*B(2)))/m;
    k6vel2 = h*q*(E(2) + ((vel(3,index-1) + vel3inc)*B(1) - (vel(1,index-1) + vel1inc)*B(3)))/m;
    k6vel3 = h*q*(E(3) + ((vel(1,index-1) + vel1inc)*B(2) - (vel(2,index-1) + vel2inc)*B(1)))/m;
    
    vel(1, index) = vel(1, index-1) + h*(1921409*k0vel1 + 9690880*k2vel1 + 13122270*k3vel1 - 5802111*k4vel1 + 1902912*k5vel1 + 534240*k6vel1)/21369600;
    vel(2, index) = vel(2, index-1) + h*(1921409*k0vel2 + 9690880*k2vel2 + 13122270*k3vel2 - 5802111*k4vel2 + 1902912*k5vel2 + 534240*k6vel2)/21369600;
    vel(3, index) = vel(3, index-1) + h*(1921409*k0vel3 + 9690880*k2vel3 + 13122270*k3vel3 - 5802111*k4vel3 + 1902912*k5vel3 + 534240*k6vel3)/21369600;
    
    
    %%
    %LOOP UPDATES
    
    %Update position
    pos(:, index) = pos(:, index-1) + vel(:, index)*h;
    
    %Update acceleration
    acc(1,index) = q*(E(1) + (vel(2,index)*B(3) - vel(3,index)*B(2)))/m;
    acc(2,index) = q*(E(2) + (vel(3,index)*B(1) - vel(1,index)*B(3)))/m;
    acc(3,index) = q*(E(3) + (vel(1,index)*B(2) - vel(2,index)*B(1)))/m;
    
    %%
    %Calculate Radius Geometrically
    dx = pos(1, index) - pos(1, index-1);
    dy = pos(2, index) - pos(2, index-1);
    dz = pos(3, index) - pos(3, index-1);
    ds = sqrt(dx^2 + dy^2 + dz^2);
    
    u = [dx dy dz];
    v1 = vel(:, index-1)';
    v2 = vel(:, index)';
    
    a1 = atan2(norm(cross(u,v1)), dot(u,v1));
    a2 = atan2(norm(cross(u,v2)), dot(u,v2));
    dtheta = a1 + a2;
    
    radius(index) = ds/dtheta;
    
    % Update the momentum magnitude at each point
    momentum_norms(index) = m*norm(vel(:,index));
    
    %Turns:
    if pos(2, index) > 0 && pos(2, index-1) < 0
        num_turns = num_turns + 1;
    end
    
    i = i + h;
    
end
    

%%
%Plotting the Data
X_min = min(pos(1, :));
X_max = max(pos(1, :));
Y_min = min(pos(2, :));
Y_max = max(pos(2, :));
Z_min = min(pos(3, :));
Z_max = max(pos(3, :));

figure;

for i = 2:50000:size(vel,2)
    plot3(pos(1, 1:i-1), pos(2, 1:i-1), pos(3, 1:i-1));
    hold on;
    plot3(pos(1,i), pos(2,i), pos(3,i), '.', 'MarkerSize', 20);
    xlim([X_min X_max]);
    ylim([Y_min Y_max]);
    %zlim([Z_min Z_max]);
    hold off;
    pause(1e-25);
end

xlabel("X-position in metres");
ylabel("Y-position in metres");
zlabel("Z-position in metres");
title("Electron Trajectory using Dormand Prince 5th Order");

pause(0.5);

figure;
plot(momentum_norms(10:end));

xlabel("Time travelled in seconds");
ylabel("Momentum of the electron in kg-m/s");
title("Momentum Scaling in Dormand Prince 5th Order");
    
    
