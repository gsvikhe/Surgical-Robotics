%% Question 1. 

cutouts_struct = struct('x',1,'h',1,'w',1.6,'orient',[0 0]);
wrist_1 = Wrist(1.6,1.8,4,cutouts_struct);
configurations= [0.0001,0,0;
                 0.1,0,0;
                 0.5,0,0;
                 0.5,20,0;
                 0.5,90,0;
                 0.5,90,5];
for i=1:6
    T_wrist_1=FwKin(wrist_1,configurations(i,:));
    pose_wrist_1(:,:,i) = permute(T_wrist_1(1:3,4,2:6),[1 3 2]);
end

x0=10;
y0=10;
width=800;
height=400
set(gcf,'units','points','position',[x0,y0,width,height])

subplot(2,3,1)
plot3(pose_wrist_1(1,:,1),pose_wrist_1(2,:,1),pose_wrist_1(3,:,1),'-o')

title('Plot for [0 0 0] configuration')
axis([-7 7 -7 7 -2 10])
grid on
view(119,19)

subplot(2,3,2)
plot3(pose_wrist_1(1,:,2),pose_wrist_1(2,:,2),pose_wrist_1(3,:,2),'-o')
title('Plot for [0.1,0,0] configuration')
axis([-7 7 -7 7 -2 10])
grid on
view(119,19)

subplot(2,3,3)
plot3(pose_wrist_1(1,:,3),pose_wrist_1(2,:,3),pose_wrist_1(3,:,3),'-o')
title('Plot for [0.5,0,0] configuration')
axis([-7 7 -7 7 -2 10])
grid on
view(119,19)

subplot(2,3,4)
plot3(pose_wrist_1(1,:,4),pose_wrist_1(2,:,4),pose_wrist_1(3,:,4),'-o')
title('Plot for [0.5,20,0] configuration')
axis([-7 7 -7 7 -2 10])
grid on
view(119,19)

subplot(2,3,5)
plot3(pose_wrist_1(1,:,5),pose_wrist_1(2,:,5),pose_wrist_1(3,:,5),'-o')
title('Plot for [0.5,90,0] configuration')
axis([-7 7 -7 7 -2 10])
grid on
view(119,19)

subplot(2,3,6)
plot3(pose_wrist_1(1,:,6),pose_wrist_1(2,:,6),pose_wrist_1(3,:,6),'-o')
title('Plot for [0.5,90,5] configuration')
axis([-7 7 -7 7 -2 10])
grid on
view(119,19)


%% Question 2. Strain Calculations
for i=1:6
    temp(i) = calcmaxstrain(wrist_1,configurations(i,:));
end
max_strain=temp
% The strain at which the max strain exceeds the recoverable strain of
% Nitinol
max_strain_ex = calcmaxstrain(wrist_1,[0.86,0,0])



% %% Question 3. Two Tendon displacement. 
% cutouts_struct_2 = struct('x',1,'h',1,'w',1.6,'orient',[0 180]);
% wrist_2 = Wrist(1.6,1.8,4,cutouts_struct_2);
% T_wrist_2=FwKin2(wrist_2,[0.46 0.46 0 0]);
% pose_wrist_2 = permute(T_wrist_2(1:3,4,2:6),[1 3 2])
% figure(2)
% plot3(pose_wrist_2(1,:),pose_wrist_2(2,:),pose_wrist_2(3,:),'-o')
% axis([-7 7 -7 7 -10 10])
% grid on
% view(119,19)
% 
% %% Note for Question 3. Tendon Displacement
% % Tried various tendon displacements but S-shape didn't form. Thought it may be wrong transformation matrix calculation
% % However confirmed the calculation by putting [0 0] in cutouts_struct_2
% % instead of [0 180]. I am not sure where I have gone wrong exactly. I
% % changed the notch structure so that notch orientations are like 
% % [0 180 0 180] (lines 218 to 221). This seems to work partially. Maybe
% % increasing the number of notches would work. 
% 
% 
% 
%% Question 3.Extra Credit. 
cutouts_struct_3 = struct('x',1,'h',1,'w',1.6,'orient',[0 90 180 270]);
wrist_3 = Wrist(1.6,1.8,8,cutouts_struct_3);
T_wrist_3=FwKin3(wrist_3,[0.6 0.6 0.6 0.6 0 0]);
pose_wrist_3 = permute(T_wrist_3(1:3,4,2:10),[1 3 2])
figure(3)
x0=10;
y0=10;
width=800;
height=400
set(gcf,'units','points','position',[x0,y0,width,height])

subplot(1,2,1)
plot3(pose_wrist_3(1,:),pose_wrist_3(2,:),pose_wrist_3(3,:),'-o')
title('Top View')
axis([-7 7 -7 7 -2 16])
grid on
view(90,90)
subplot(1,2,2)
plot3(pose_wrist_3(1,:),pose_wrist_3(2,:),pose_wrist_3(3,:),'-o')
title('View showing the helical Structure')
axis([-7 7 -7 7 -2 16])
grid on
view(102,34)

%% Note about extra credits
% This wrist design has eight notches. Although this has a helical structure
% The design could be improved by changing the intermediate notch heights
% and distance between intermediate notches except the first and last.
% A spiral structure can be achieved by reducing the notch height and
% inter notch distance from the base to the tip.
