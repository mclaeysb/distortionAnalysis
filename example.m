% This script launches a distortion analysis based on two sets of ground control points
% Written by Manuel Claeys Bouuaert, 2015

clear
if(~isdeployed)
    addpath(genpath('.'))
end

% Load ground control points 
points = dlmread('data_sample/points.txt')';
gcps_d = points(4:5,:);
gcps_c = points(2:3,:);

% Set preferences
spatRes_d = 500;                % Spatial Resolution in the domain when creating mesh. 500 for sample data
spatBuffer_d = 8000;            % Spatial buffer for mesh in the codomain. 8000 for sample data
spatResType = 'absolute';       % Are spatRes_d and spatBuffer_d input 'absolute' or 'relative' numbers (see how those are treated in distortionAnalysis)
scalingReference = 'helmert';   % 'none' or 'helmert'

doDisplacementVectors = 1;              % Should warped grid be written out?
doDistortionGrid = 1;                   % Should warped grid be written out?
doDifferentialDistortionAnalysis = 1;   % Should Differential Distortion Analysis be written out?
doIndicatrices = 1;                     % Should Tissot indicators be written out?

doPlots = 0;    % Should plots be made?

% Perform analysis
% distortionAnalysis(gcps_d,gcps_c) % Simplest call to distortionAnalysis.m. Default values specified in distortionAnalysis.m will be used in stead of used defined preferences.
distortionAnalysis(gcps_d,gcps_c,spatRes_d,spatBuffer_d,spatResType,scalingReference,doDisplacementVectors,doDistortionGrid,doDifferentialDistortionAnalysis,doIndicatrices,doPlots)
