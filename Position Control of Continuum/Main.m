%% Advanced Surgical Robotics - Fall 2018
%  Homework 3 - Problem 3 - Numerical inverse kinematics

close all
clear
clc

%% Create a robot with 3 links
nLinks = 3;
OD = [1.8 1.6 1.4]; % link diameters [mm]

r = Robot(nLinks, OD);

%% Create a configuration c and plot the robot pose
kappa1 = 0.001;    % curvature of link 1 [m^-1]
phi1   = 0;     % base rotation of link 1 [rad]
l1     = 1;  % length of link 1 [m]

kappa2 = 0.05;    % curvature of link 2 [m^-1]
phi2   = 0;     % base rotation of link 2 [rad]
l2     = 1;  % length of link 2 [m]

kappa3 = 0.2;    % curvature of link 3 [m^-1]
phi3   = 0;     % base rotation of link 3 [rad]
l3     = 4;  % length of link 3 [m]

c = [kappa1 phi1 l1 kappa2 phi2 l2 kappa3 phi3 l3]';
[T, links] = r.fwkinematics(c);

figure, axis equal, grid on, hold on
colorPalette = distinguishable_colors(nLinks);
worldRef = triad('Scale', 3, 'linewidth', 5);
robotRef = triad('matrix', T, 'Scale', 3, 'linewidth', 5);
xlabel('X [mm]', 'fontsize', 16)
ylabel('Y [mm]', 'fontsize', 16)
zlabel('Z [mm]', 'fontsize', 16);
xlim([-7 7]), ylim([-7 7]), zlim([0 20]);
title('Position control of a 3-link continuum robot', 'fontsize', 16);
view(110.4, 33.2);

s = cell(nLinks);

for jj = 1 : nLinks
    [X, Y, Z] = gencyl(links{jj}, OD(jj) / 2 * ones(1,size(links{jj}, 2)));
    s{jj} = surf(X, Y, Z);
    set(s{jj}, 'FaceColor', colorPalette(jj,:), 'FaceAlpha', 1);
end

%% Define a path in task space

t = 10:pi/5:10*pi;
st = linspace(1,5,size(t,2)).*sin(t);
ct = linspace(1,5,size(t,2)).*cos(t);

trajectory = [st; ct; t/2];

% Plot the trajectory
scatter3(trajectory(1,:), trajectory(2,:), trajectory(3,:));
all_config= zeros(9,size(trajectory,2));

% Perform inverse kinematics on each of the waypoints
for ii = 1 : size(trajectory, 2)
    
    % get the next waypoint
    pTarget = trajectory(:,ii);
    
    % numeric inverse kinematics
    while true
        % calculate the current location
        [T, links] = r.fwkinematics(c);
        pCurrent = T(1:3,end);
        
        % calculate the difference between current and target location
        err = norm(pTarget - pCurrent);
        
        % if the difference < epsilon, return
        if err < 1e-1, break; end
        
        % else, update the "joint" variables using the inverse of the
        % jacobian
        J = r.jacob(c);
        Jp = J(1:3,:) - skew(pCurrent) * J(4:6,:); 
        Jpinv = pinv(Jp);
        
        K = diag([1 1 1 1 1 1 0 1 1]);
        deltaQ = K*Jpinv * (pTarget - pCurrent);
        c = c + deltaQ;  
        all_config(:,ii)=c;
    end
end
for i=1:size(trajectory,2)-1
     interpol_configs(:,:,i)=(interp1([1 2],[all_config(:,i)';all_config(:,i+1)'],linspace(1,2,10)))';
end
interpol_configs=reshape(interpol_configs,[9,(size(trajectory,2)-1)*10]);
for j=1:size(interpol_configs,2)
    % plot
    [~, links] = r.fwkinematics(interpol_configs(:,j));
    
    for jj = 1 : nLinks
        [X, Y, Z] = gencyl(links{jj}, OD(jj) / 2 * ones(1, size(links{jj}, 2)));
        set(s{jj}, 'XData', X);
        set(s{jj}, 'YData', Y);
        set(s{jj}, 'ZData', Z);
    end
    
    robotRef.Matrix = T;
    
    drawnow();
end

%% Comment on K matrix
% The K matrix is a gain matrix controlling whether a change in a
% particular dimension would be applied or not. So in our case of the 7th
% element being one, the K matrix is not letting the curvature of the
% third link change. Changing it to one, makes the last link more flexible
% and leads to a different configuration being followed. To illustrate
% further, making the 6th or the 3rd element zero would lead to no change
% in the original lengths of the 2nd and 1st elements respectively. 