function [ SonarFeatures ] = GetSonarFeatures( SonarValues )
%GetSonarFeatuers Summary of this function goes here
%   Detailed explanation goes here
% SegmentedBinaryImage contains the binary segmented image 
%   1. Left Minimum Distance 
%   2. Left Minimum Distance - Angle
%   3. Right Minimum Distance 
%   4. Right Minimum Distance - Angle
%   5. Critical Minimum Distance 
%   6. Critical Minimum Distance - Angle
%   7. Safe Direction

S = SonarValues;
AngRight = [-50 -90];
AngLeft = [50 90];
AngCrit = [30 10 -10 -30];

% right Zone  - 0 to 60

[rightMinDist,Ind] = min(S(1,7:8));
rightMinAngle = AngRight(Ind);
   
% left Zone  -  120 to 180

[leftMinDist,Ind] = min(S(1,1:2));
leftMinAngle = AngLeft(Ind);

% crtical Zone - 60 to 120

[critMinDist,Ind] = min(S(1,3:6));
critMinAngle = AngCrit(Ind);

[Safe,Ind] = max(S(1,3:6));
safeDirection = AngCrit(Ind);


% aggregate the features
SonarFeatures = [leftMinDist leftMinAngle rightMinDist rightMinAngle critMinDist critMinAngle safeDirection];

end

