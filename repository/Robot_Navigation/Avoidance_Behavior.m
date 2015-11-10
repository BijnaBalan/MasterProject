%Example of using Dr. Hoffmann's avoidance behavior controller
%This simulation will simulate the robot for 60 seconds, record robot data, and display
%some graphs at the end

%INSTRUCTIONS:
%Make sure "src/pkg" is in your matlab path%!!
%Run this script

close all;
clear all;

%Show the robot position on the graph in real-time as it drives?
SHOW_LIVE_DATA = false;

%Create a simulation instance
s = Demonstration.Simulator1();
mapPath = 'F:\Dropbox\Share_Rashik\temp\simDataCentre\deadEnd_FH\1.map';
s.buildWorldFromMobilsimMap(mapPath);

% Load the default map
%---------------i m commenting
% s.buildTestWorld1();


% %add small isolated obstacles
% isolatedObstacles=  [1 1; 1 5; 3 5; 5 5; 5 4; 5 3; 9 5; 3 1; 4 1; 8 4; 7 1; 9 1];
% for ii = 1:size(isolatedObstacles, 1)
%     s.addPointObstacle(isolatedObstacles(ii, 1), isolatedObstacles(ii, 2));
% end

%-------------------------------i m commenting

 %step the simulation one step with no motion to perform raycasting
 %to generate initial ultrasound readings
s.setVelocity(0);
s.setAngularVelocity(0);
s.step()   

%create an instance of Hoffmann's avoidance behavior controller
controller = Common.UsConAvoidanceBehv1();

% possible startpoints for the robot (to be chosen randomly)
startpoints = [-16   4;
             -10   -5  ;
             -5   -5;
             -10   3;
             -4    5;
              0    0;
             -16.5 5;
              3   -2;
             -1    4]'; 

%initial robot position from one of the random startpoints, and random angle
% s.setPose([startpoints(:,randi(length(startpoints))); rand(1)*2*pi-pi]);

s.setPose([-10 0 0]);

%get a plot handle and visualize walls
fig = figure(1);
s.plotWorld(fig)


disp('Starting simulation');
stepCount = 1;
while stepCount <= 250/0.1       %simulate for 50 seconds
    %give the controller inputs and have it process them
    controller.processInputs(s.getUsDistancesRaw(), s.getPose(), s.getVelocity(), s.getAngularVelocity());
    %get output velocities from the controller
    s.setVelocity(controller.getVelocity());
    s.setAngularVelocity(controller.getAngularVelocity());

    %report if the robot is in a collision state
    if s.getCollisionState()
        disp('Collision')
    end

    s.step();   %advance the simulation one time step with given inputs

    if SHOW_LIVE_DATA
        s.plotRobotAndSonar(fig);   %plot robot position and sonar raycasts
    end

    %collect position data
    poseData{stepCount, 1} = s.getPose();

    %collect velocity and steering angle
    velocityData{stepCount, 1} = s.getVelocity();
    angVelocityData{stepCount, 1} = s.getAngularVelocity();

    %collect ultrasound data
    usData{stepCount, 1} = s.getUsDistances();          %sampled with sensor model
    usDataRaw{stepCount, 1} = s.getUsDistancesRaw();    %unsampled raw distances
    stepCount = stepCount + 1;
end


%convert position, velocity and steering angle data from cell array to matrix
poseData = cell2mat(poseData);
velocityData = cell2mat(velocityData);
angVelocityData = cell2mat(angVelocityData);

%convert ultrasound data from cell arrays to matrices
usData = cell2mat(usData);
usDataRaw = cell2mat(usDataRaw);

%visualize robot path
%---------------------------------------i m commeting-------------
hold on;
plot(poseData(:,1), poseData(:,2));

%plot ultrasound channels
figure(2);
subplot(2,2,1);
hold on;
plot(usDataRaw(:,1), 'or')      %raw raycasting distances (not using statistical model)
plot(usDataRaw(:,4), 'og')
plot(usDataRaw(:,8), 'ob')
title('Sonar distances (non-stochasitic)')
ylabel('Distance [meters]')
xlabel('Sample (0.1 second dt)')

subplot(2,2,2);
hold on;
plot(usData(:,1), 'or')
plot(usData(:,4), 'og')
plot(usData(:,8), 'ob')
title('Sonar distances (stochasitic)')
ylabel('Distance [meters]')
xlabel('Sample (0.1 second dt)')

%plot robot angle
subplot(2,2,3);
plot(poseData(:,3).*(180/pi));
title('Robot angle')
ylabel('Angle [deg]')
xlabel('Sample (0.1 second dt)')

%plot robot steering angle and velocity
subplot(2,2,4);
hold on;
plot(angVelocityData(:,1).*(180/pi));
title('Robot angular velocity')
ylabel('Anglular Velocity [deg/s]')
xlabel('Sample (0.1 second dt)')
%--------------------i m commenting----------------------------------------

usData = 1000.*usData;
velocityData = 1000.*velocityData;
angVelocityData = rad2deg(angVelocityData);
poseData(:,1:2) = 1000.*poseData(:,1:2);
poseData(:,3) = rad2deg(poseData(:,3));
SonarRawData = [usData poseData velocityData angVelocityData];
Time = clock;
sonarfilename =  ['SonarReadingsSimulated' num2str(date) '_' num2str(Time(4)) num2str(Time(5)) '.mat'];
save(sonarfilename,'SonarRawData')



