% This script launches a distortion analysis or triangle area scaling analsis based on two sets of ground control points
% Written by Manuel Claeys Bouuaert, 2015

clear
if(~isdeployed)
    addpath(genpath('.'))
end

% Load ground control points 
pathToData = '/Users/Manuel/Google Drive/PhD Shared Ferraris';
% Ferraris data copyForHPC
%{
gcps_d = dlmread('data_aquaterra_copyForHPC/real_cass_noexcl.csv')';
gcps_c = dlmread('data_aquaterra_copyForHPC/old_ferrAqLam_noexcl.csv')';
%}
% Ferraris data copyForHPC, reverse

gcps_d = dlmread('data_aquaterra_copyForHPC/old_ferrAqLam_noexcl.csv')';
gcps_c = dlmread('data_aquaterra_copyForHPC/real_lam_noexcl.csv')'; % real_lam, not real_cass!

% Ferraris data
%{
gcps_d = dlmread(fullfile(pathToData,'Data_DGARNE/data/real_cass.csv'))';
gcps_c = dlmread(fullfile(pathToData,'Data_DGARNE/data/old_ferrAqLam.csv'))';
%}
% DGARNE data
%{
gcps_d = dlmread(fullfile(pathToData,'Data_DGARNE/data/real_cass.csv'))';
gcps_c = dlmread(fullfile(pathToData,'Data_DGARNE/data/old_ferrAqLam.csv'))';
%}
% Append data
%{
gcps_d = dlmread(fullfile(pathToData,'Data_Append/data/real_cass.csv'))';
gcps_c = dlmread(fullfile(pathToData,'Data_Append/data/old_ferrAqLam.csv'))';
%}
% MapAnalyst-Basel data
%{
points = dlmread(fullfile(pathToData,'MapAnalyst_Basel/MapAnalystDataEditted/points.txt'))';
gcps_d = points(4:5,:);
gcps_c = points(2:3,:);
%}
% MapAnalyst-Basel data sample
%{
points = dlmread('data_sample/points.txt')';
gcps_d = points(4:5,:);
gcps_c = points(2:3,:);
%}
% Geraardsbergen data (in both directions!)
%{
gcps_d = dlmread(fullfile(pathToData,'Geraardsbergen/data/CdC_old_ferrAqLam.csv'))';
gcps_c = dlmread(fullfile(pathToData,'Geraardsbergen/data/CM_old_cassToise.csv'))';
%}
%{
% Meetpunten data
data = dlmread(fullfile(pathToData,'Meetpunten/Alle_meetpunten.csv'),';',2,3); % Alternative loading code: fid = fopen('../Meetpunten/Alle_meetpunten.csv'); comment = textscan(fid,'%s',1, 'delimiter',';'); headers = textscan(fid,'%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s',1, 'delimiter',';'); data = textscan(fid,'%d %s %s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f','delimiter',';'); fclose(fid)
data = data(:,[1:4, 7:8, 11:28]); % Column selection
% 1:2   = mod topo - Cstad
% 3:4   = mod topo - CdT
% 5:6   = 1775
% 7:8   = NOUV
% 9:10   = TAB
% 11:12  = CM gen
% 13:14 = CM kb - Cstad
% 15:16 = CM kb - CdT
% 17:18 = CdC - Cstad
% 19:20 = CdC - Cstad - Helmert to mod topo
% 21:22 = CdC - CdT
% 23:24 = CdC - CdT - Helmert to mod topo
% Option 1: select columns here, and select rows based on shared non-zero's
data = data(:,[15:16,3:4]); % Only if Option 1
data(any(data==0,2),:)=[]; % Always: select only rows with no zero's
% Option 2: select columns here (if option 1, then this is 1:2 and 3:4)	
gcps_d = data(:,1:2)';
gcps_c = data(:,3:4)';
%}

% Set preferences
spatRes_d = 500; % Spatial Resolution in the domain when creating mesh. Default: 500
spatBuffer_d = 8000; % Spatial buffer for mesh in the codomain. Default: 8000
spatResType = 'absolute'; % Are spatRes_d and spatBuffer_d input 'absolute' or 'relative' numbers
scalingReference = 'helmert'; % 'none' or 'helmert'

doDisplacementVectors = 1; % Should warped grid be written out?
doDistortionGrid = 1; % Should warped grid be written out?
doDifferentialDistortionAnalysis = 1; % Should Differential Distortion Analysis be written out?
doIndicatrices = 0; % Should Tissot indicators be written out?

doTriangles = 0; % Should Triangle Area Scale Analysis be written out?

doPlots = 0; % Should plots be made?

% Perform analysis
distortionAnalysis(gcps_d,gcps_c,spatRes_d,spatBuffer_d,spatResType,scalingReference,doDisplacementVectors,doDistortionGrid,doDifferentialDistortionAnalysis,doIndicatrices,doPlots)

if doTriangles
	triangleAreaScaleAnalysis(gcps_d,gcps_c,scalingReference,doPlots)
end
