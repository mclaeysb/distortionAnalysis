function [ps_ref_to_target,a,translation,rotation,scaling,v,s0] = helmert(ps_ref,ps_target,verbose,doPlots)

% This script calculates the Helmert transform parameters of two datasets
% ps_... inputs must both have dimensions 2*nPoints
% Written by Manuel Claeys Bouuaert, 2015

% Process data input
if nargin == 2
    verbose = 0; doPlots = 0;
end
if size(ps_ref)~=size(ps_target)
    error('Loaded reference and target points have different dimensions')
end
nPoints=size(ps_target,2);
if verbose
    disp(['Performing Helmert transform with ',num2str(nPoints),' points'])
end

% Compute transformation parameters
l = [ps_target(1,:),ps_target(2,:)]';
A = [[ones(nPoints,1),zeros(nPoints,1),ps_ref(1,:)',-ps_ref(2,:)'];[zeros(nPoints,1),ones(nPoints,1),ps_ref(2,:)',ps_ref(1,:)']];
a = pinv(A)*l;

% Compute location of transformed ref points to target domain
ps_ref_to_target = [a(1)+a(3)*ps_ref(1,:)-a(4)*ps_ref(2,:) ; a(2)+a(4)*ps_ref(1,:)+a(3)*ps_ref(2,:)];
% dlmwrite('ps_ref_to_target.csv',ps_ref_to_target') % write out

% Compute translation, rotation, scaling and average error
translation = a(1:2);
rotation = atan(a(4)/a(3));
scaling = sqrt(a(4)^2+a(3)^2);

% compute distance errors and mean error
v = ps_target-ps_ref_to_target;
s0 = sqrt(norm(v,'fro')^2/(2*nPoints-4));

if verbose
    a,translation,rotation,scaling
end

if doPlots
    disp('Plotting target points (black) and reference points transformed to target domain (blue).')
    clf; hold on;
    axis equal;
    plot(ps_target(1,:)',ps_target(2,:)','.k')
    plot(ps_ref_to_target(1,:)',ps_ref_to_target(2,:)','.b')
    % plot(ps_ref(1,:)',ps_ref(2,:)','.r')
    hold off;
end

end
