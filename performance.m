% This script computes the performance of the interpolation as function of the amount of gcps
% We could do repetetions to sample CV-like, but assume the dataset is large enough for a sample to be random.
% Written by Manuel Claeys Bouuaert, 2015

clear
rng(1)

% addpath(genpath(fullfile('.','functions')))

% Load points
% Ferraris data
gcps_fd = dlmread('data_aquaterra_copyForHPC/real_cass_noexcl.csv')'; % real_cass
gcps_td = dlmread('data_aquaterra_copyForHPC/old_cass_noexcl.csv')'; % old_cass
% MapAnalyst-Basel data
%{
Basel_points = dlmread('../MapAnalyst_Basel/MapAnalystDataEditted/points.txt')'; % real_cass; % Use in order to load dummy data
gcps_fd = Basel_points(4:5,:);
gcps_td = Basel_points(2:3,:);
%}
%
nPoints = size(gcps_fd,2);

% Make train and testset
testPoints = randsample(1:nPoints,floor(nPoints/10));
gcps_fd_test = gcps_fd(:,testPoints);
gcps_td_test = gcps_td(:,testPoints);
nTestPoints = length(testPoints);
trainPoints = setdiff(1:nPoints,testPoints);
gcps_fd_train = gcps_fd(:,trainPoints);
gcps_td_train = gcps_td(:,trainPoints);
nTrainPoints = length(trainPoints);
disp(['Using ',num2str(nTestPoints),' testPoints on a total of ',num2str(nPoints),' points.'])

% Iterate and compute mean error
quantileProb = [0.05 0.2 0.5 0.8 0.95];
nSteps = 20;
meanErrors = zeros(1,nSteps);
quantiles = zeros(length(quantileProb),nSteps);
nTrainPoints = floor((1:nSteps)/nSteps*nTrainPoints);
for i=1:nSteps
    disp(['At step ',num2str(i),' of ',num2str(nSteps),' building tps from ',num2str(nTrainPoints(i)),' trainPoints.'])
    tic
    
    trainPoints_i = randsample(trainPoints,nTrainPoints(i));
    gcps_fd_train_i = gcps_fd(:,trainPoints_i);
    gcps_td_train_i = gcps_td(:,trainPoints_i);

    tps_train_i = tpaps(gcps_fd_train_i, gcps_td_train_i, 1);
    gcps_td_test_i = fnval(tps_train_i,gcps_fd_test);
    errors_i = sqrt(sum((gcps_td_test-gcps_td_test_i).^2,1));
    meanErrors(i) = mean(errors_i);
    quantiles(:,i) = quantile(errors_i,[0.025 0.25 0.50 0.75 0.975])';
    
    disp(['Mean error of ',num2str(meanErrors(i)),' computed'])
    
    toc
end
dlmwrite('output/meanErrors.txt',[nTrainPoints',meanErrors',quantiles'])

% % Plot results
% clf; hold on;
% axis([0 nSteps 0 max(max(quantiles))]);
% plot(meanErrors','k');
% plot(quantiles','b');
% hold off;
