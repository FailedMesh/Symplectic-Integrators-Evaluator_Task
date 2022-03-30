clear; clc;

%%
%INITIALIZATION
E = [0; 0; 0];  %Set Electric Field
B = [0; 0; 1];  %Set Magnetic Field

c = 3e+8;       %Speed of Light
m = 9.109e-31;  %Mass of electron
q = -1.60218e-19;   %Charge of electron

h = 1e-16;   %Time Step

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
%LOOP IN Leapfrog Method

while num_turns < max_turns
    index = int32((1/h)*(i) + 1);
    
    vel(:, index) = vel(:, index-1) + acc(:, index-1)*h;
    
    acc(1,index) = q*(E(1) + (vel(2,index)*B(3) - vel(3,index)*B(2)))/m;
    acc(2,index) = q*(E(2) + (vel(3,index)*B(1) - vel(1,index)*B(3)))/m;
    acc(3,index) = q*(E(3) + (vel(1,index)*B(2) - vel(2,index)*B(1)))/m;
    
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
title("Electron Trajectory using Leapfrog Method");

pause(0.5);

figure;
plot(momentum_norms(10:end));
    
xlabel("Time travelled in seconds");
ylabel("Momentum of the electron in kg-m/s");
title("Momentum Scaling in Leapfrog Method"); 
