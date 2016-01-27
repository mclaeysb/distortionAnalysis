function distortionAnalysis(gcps_d,gcps_c,varargin)

% This function performs a distortion analysis on two sets of ground control points
% It takes a list of ground control points and constructs a mapping between them using a thin plate spline
% Then, it calculates distortion measures on a mesh and in the control points
% Distortion measures include: Displacement Vectors, Distortion Grid, Indicatrices of Tissot and Differential Distortion Analysis
% 
% This function requires the Octave splines package or Matlab CurveFitting Toolbox
% Install the splines package in Octave using: pkg install -forge -global -auto splines
% 
% gcps_d and gcps_c are expected to be 2-by-n matrices, where n is the number of ground control points
% 
% Nomenclature:
% '..._d' = points in the 'domain' = 'from'-domain = 'centers' = 'reference'.
% '..._c' = points in the 'codomain' = 'to'-domain = 'values' = 'target'.
% 'gcps' = ground control points
% 'mps' = mesh points
% 'aps' = all points (mps + gcps)
% 
% Computed distortion values are displayed in codomain
% If more than 729 points are loaded, Matlab's Thin Plate Spline fuction tpaps() will compute the spline approximately and no longer exactly
% To prevent this, use an edited version of tpaps.m or use the Octave splines package
% 
% Written by Manuel Claeys Bouuaert, 2015

% Detect if Octave or Matlab is used
isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;
if isOctave
    pkg load splines
    page_output_immediately(1, 'local');
end

% Process data input: overwrite default values with user specified values
numvarargs = length(varargin);
optargs = {100 20 'relative' 'helmert' 1 1 1 1 0};
if numvarargs > length(optargs)
    error ('Too many input arguments')
end
[optargs{1:numvarargs}] = varargin{:};
[spatRes_d,spatBuffer_d,spatResType,scalingReference,doDisplacementVectors,doDistortionGrid,doDifferentialDistortionAnalysis,doIndicatrices,doPlots] = optargs{:};

% Check ground control points
if any(size(gcps_d)~=size(gcps_c))
    error('Loaded ground control points have different dimensions')
elseif size(gcps_d,1)~=2
    error('Input is expected to be 2-by-nGcps matrices')
elseif size(gcps_d,2)==0
    error('Empty list of ground control points')
end
nGcps = size(gcps_d,2);
gcps_row = zeros(1,nGcps); gcps_col = zeros(1,nGcps);
gcps_type_cell=cell(nGcps,1);
gcps_type_cell(:)={'gps'};
disp(['Loaded ',num2str(nGcps),' ground control points (gcps)'])

% Get scaling parameter and displacement points from Helmert transform
[gcps_d_to_c,~,~,~,scaling_helmert,~,~]=helmert(gcps_d,gcps_c);
% Set reference scaling (for reference value of a, b, sigma)
if strcmp(scalingReference,'none')
    scaling_ref = 1;
elseif strcmp(scalingReference,'helmert')
    scaling_ref = scaling_helmert;
else
    error('Unknown value for scalingReference')
end
areascaling_ref = scaling_ref^2;
% Process resolution
if strcmp(spatResType,'relative')
    % Compute absolute spatRes_d and spatBuffer_d values based on relative input and length of largest dimension
    lengthLongest_d = (max(max(gcps_d(1,:))-min(gcps_d(1,:)),max(gcps_d(2,:))-min(gcps_d(2,:))));
    spatRes_d = lengthLongest_d/spatRes_d;
    spatBuffer_d = lengthLongest_d/spatBuffer_d;
    % From this point on, spatRes_d and spatBuffer_d are absolute
elseif ~strcmp(spatResType,'absolute')
    error('Unknown value for scalingReference')
end
spatRes_c = spatRes_d*scaling_helmert;
spatBuffer_c = spatBuffer_d*scaling_helmert;
disp(['A Helmert transform from domain to codomain yields scaling of ',num2str(scaling_helmert),''])

% Only continue if needed
if doDistortionGrid|doDifferentialDistortionAnalysis|doIndicatrices
	% Compute thin plate spline function
	disp(['A reference scaling of ',num2str(scaling_ref),' will be used for distortion measures'])
	disp(['Building Thin Plate Spline from ',num2str(nGcps),' gcps'])
	tic
	% Note: in tpaps, in both Matlab and Octave, the third parameter is the smoothing parameter p, where p=0 means smoothing is applied, and p=1 means no smoothing is applied
	if isOctave
	    % Spline stored using coefficients
	    tps_coefs = zeros(2,nGcps+2+1);
	    for i=1:2
	        tps_coefs(i,:) = tpaps(gcps_d',gcps_c(i,:)',1,[])';
	    end
	else
	    % Spline is stored in stform
	    if nGcps>=729 		% turn this 'magic number' into a constant instead
	        warning('You are using 729 control points or more. Unless you are using an adapted version of tpaps.m, the Thin Plate Spline will be solved approximately, possibly taking a lot of time and yielding inadequate results.')
	    end
	    tps = tpaps(gcps_d, gcps_c, 1);
	end
	toc

	% Create mesh in domain: replicate grid vectors xgv and ygv to produce a full grid
	disp(['Creating mesh in domain'])
	disp(['It has a spatial resolution of ',num2str(spatRes_d),' units, and a buffer of ',num2str(spatBuffer_d),' units'])
	tic
	xgv_d = (min(gcps_d(1,:))-spatBuffer_d):spatRes_d:(max(gcps_d(1,:))+spatBuffer_d);
	ygv_d = (min(gcps_d(2,:))-spatBuffer_d):spatRes_d:(max(gcps_d(2,:))+spatBuffer_d);
	[mps_xgrid_d, mps_ygrid_d] = meshgrid(xgv_d, ygv_d);
	nRows = size(mps_xgrid_d,1);
	nCols = size(mps_xgrid_d,2);
	nMps = numel(mps_xgrid_d);
	mps_d = [reshape(mps_xgrid_d,1,nMps);reshape(mps_ygrid_d,1,nMps)];
	[mps_row,mps_col] = ind2sub([nRows,nCols],1:nMps);
	mps_type_cell=cell(nMps,1);
	mps_type_cell(:)={'mps'};
	disp(['Created ',num2str(nMps),' mesh points (mps)'])
	toc

	% Transpose mesh to codomain 
	disp(['Transposing mesh to codomain'])
	disp(['It has an average spatial resolution of ',num2str(spatRes_c),' units, and a buffer of ',num2str(spatBuffer_c),' units'])
	tic
	if isOctave
	    mps_c = zeros(2,nMps);
	    for i=1:2
	        mps_c(i,:) = tps_val(gcps_d',tps_coefs(i,:)',mps_d')';
	    end
	else
	    mps_c = fnval(tps,mps_d);
	end
	mps_xgrid_c = reshape(mps_c(1,:),size(mps_xgrid_d));
	mps_ygrid_c = reshape(mps_c(2,:),size(mps_ygrid_d));
	toc
end

% Only continue if needed
if doDifferentialDistortionAnalysis|doIndicatrices
	% Concatenate gcps and mps
	disp(['Concatenating all points (mps + gcps = aps)'])
	tic
	aps_d = [mps_d,gcps_d];
	aps_c = [mps_c,gcps_c];
	aps_row = [mps_row,gcps_row];
	aps_col = [mps_col,gcps_col];
	aps_type_cell=[mps_type_cell;gcps_type_cell];
	nAps = nMps + nGcps;
	disp(['Total of ',num2str(nAps),' points (aps)'])
	toc

	% Perform Differential Distortion Analysis
	disp(['Evaluating derivatives in aps'])
	tic
	if isOctave
	    aps_d_pd10 = zeros(2,nAps);
	    aps_d_pd01 = zeros(2,nAps);
	    for i = 1:2
	        aps_d_pd = tps_val_der(gcps_d',tps_coefs(i,:)',aps_d')';
	        aps_d_pd10(i,:) = aps_d_pd(1,:);
	        aps_d_pd01(i,:) = aps_d_pd(2,:);
	    end
	else
	    tps_10 = fnder(tps,[1,0]);
	    tps_01 = fnder(tps,[0,1]);
	    aps_d_pd10 = fnval(tps_10,aps_d);
	    aps_d_pd01 = fnval(tps_01,aps_d);    
	end
	toc
	disp(['Calculating distortion in aps'])
	tic
	aps_c_E = sum(aps_d_pd10.^2,1);
	aps_c_G = sum(aps_d_pd01.^2,1);
	aps_c_F = sum(aps_d_pd10 .* aps_d_pd01,1);
	aps_c_a = sqrt(0.5*(aps_c_E+aps_c_G+sqrt((aps_c_E-aps_c_G).^2+4*aps_c_F.^2)));
	aps_c_b = sqrt(0.5*(aps_c_E+aps_c_G-sqrt((aps_c_E-aps_c_G).^2+4*aps_c_F.^2)));
	aps_c_thetaxp = atan(aps_d_pd10(2,:)./aps_d_pd10(1,:));
	% aps_c_thetayp = atan(tps_01_aps_d(2,:)./tps_01_aps_d(1,:));
	% aps_c_tantheta = sqrt(aps_c_E.*aps_c_G-aps_c_F.^2)./aps_c_F;
	% aps_c_theta = aps_c_thetayp-aps_c_thetaxp;
	% aps_c_theta2 = atan(aps_c_tantheta); % theta and theta2 are equal under mod pi
	aps_c_alphap = sign(-aps_c_F).*asin(sqrt((1-(aps_c_a.^2./aps_c_E))./(1-(aps_c_a./aps_c_b).^2)));
	aps_c_log2sigma = log(aps_c_a.*aps_c_b/areascaling_ref)/log(2); % take log_2(x) for better interpretability
	aps_c_2Omega = 2*asin((aps_c_a-aps_c_b)./(aps_c_a+aps_c_b));
	% aps_c_Airy1 = 0.5*((aps_c_a./aps_c_b-1).^2+(aps_c_a.*aps_c_b/areascaling_ref-1).^2);
	% aps_c_Airy2 = 0.5*((aps_c_a/scaling_ref-1).^2+(aps_c_b/scaling_ref-1).^2);
	aps_c_AiryKavr = 0.5*(log(aps_c_a/scaling_ref).^2+log(aps_c_b/scaling_ref).^2);
	aps_c_signDetJ = sign(aps_d_pd10(1,:).*aps_d_pd01(2,:)-aps_d_pd10(2,:).*aps_d_pd01(1,:));
	aps_c_thetaa = aps_c_thetaxp - aps_c_alphap;
	toc
end

% Create output folder if it doesn't exists
nameOutputFolder = 'output';
if ~exist(nameOutputFolder, 'dir')
  mkdir(nameOutputFolder);
  addpath(nameOutputFolder);
end

% Write out displacement vectors
if doDisplacementVectors
	disp(['Writing out displacement vectors'])
	tic
	gcps_c_displacementVectors_cell = num2cell([1:nGcps;gcps_c;gcps_d_to_c]');
    file_displacementVectors = fopen('output/displacementVectors.wkt','w');
	fprintf(file_displacementVectors,'%s\r\n','number;Line');
	for i = 1:size(gcps_c_displacementVectors_cell,1)
	    fprintf(file_displacementVectors,'%1.0f;LINESTRING(%9.3f %9.3f, %9.3f %9.3f)\r\n',gcps_c_displacementVectors_cell{i,:});
	end
	fclose(file_displacementVectors);
	toc
end

% Write out distortion grid
if doDistortionGrid
	disp(['Writing out distortion grid'])
	tic
	file_distortionGrid = fopen('output/distortionGrid.wkt','w');
	fprintf(file_distortionGrid,'%s\r\n','type;number;Line');
	for i=1:nCols
	    fprintf(file_distortionGrid,'%s;','row');
	    fprintf(file_distortionGrid,'%d;',i);
	    fprintf(file_distortionGrid,'%s','LINESTRING');
	    str = mat2strWKT(mps_c(:,((i-1)*nRows+1):(i*nRows))');
	    fprintf(file_distortionGrid,'%s\r\n',str);
	end
	for i=1:nRows
	    fprintf(file_distortionGrid,'%s;','col');
	    fprintf(file_distortionGrid,'%d;',i);
	    fprintf(file_distortionGrid,'%s','LINESTRING');
	    str = mat2strWKT(mps_c(:,i:nRows:((nCols-1)*nRows+i))');
	    fprintf(file_distortionGrid,'%s\r\n',str);
	end
	fclose(file_distortionGrid);
	toc
end

% Write out Differential Distortion Analysis
if doDifferentialDistortionAnalysis
    disp(['Writing out Differential Distortion Analysis'])
    tic
    aps_c_differentialDistortionAnalysis_cell = [num2cell([aps_c;aps_row;aps_col;aps_c_log2sigma;aps_c_2Omega;aps_c_AiryKavr;aps_c_signDetJ;aps_c_a;aps_c_b;aps_c_thetaa]') aps_type_cell];
    file_differentialDistortionAnalysis = fopen('output/differentialDistortionAnalysis.csv','w');
    fprintf(file_differentialDistortionAnalysis,'%s\r\n','x_c,y_c,row,col,log2sigma,2Omega,AiryKavr,signDetJ,a,b,thetaa,type');
    for i = 1:size(aps_c_differentialDistortionAnalysis_cell,1)
        fprintf(file_differentialDistortionAnalysis,'%8.6f,%8.6f,%d,%d,%4.3f,%4.3f,%4.3f,%d,%3.2f,%3.2f,%3.2f,"%s"\r\n',aps_c_differentialDistortionAnalysis_cell{i,:});
    end
    fclose(file_differentialDistortionAnalysis);
    toc
end
    
% Write out Tissot indicatrices
if doIndicatrices
	disp(['Writing out Tissot indicatrices'])
	tic
	file_indicatrices = fopen('output/indicatrices.wkt','w');
	fprintf(file_indicatrices,'%s\r\n','row;col;Polygon');
	for i=1:nMps % Only run over mps = first nMps points
        fprintf(file_indicatrices,'%d;',aps_row(i));
        fprintf(file_indicatrices,'%d;',aps_col(i));
        fprintf(file_indicatrices,'%s','POLYGON(');
        t = [0:0.1:2*pi,2*pi];
        loc = repmat(aps_c(:,i),1,length(t));
        sca = (spatRes_c/4);
        rot = [cos(aps_c_thetaa(i)),-sin(aps_c_thetaa(i));sin(aps_c_thetaa(i)),cos(aps_c_thetaa(i))];
        el = (1/scaling_ref)*[aps_c_a(i)*cos(t);aps_c_b(i)*sin(t)];
        indicator = loc+sca*rot*el;
        str = mat2strWKT(indicator',10);
        fprintf(file_indicatrices,'%s)\r\n',str);
    end
	fclose(file_indicatrices);
	toc
end

% Plot
if doPlots
	if not(doDifferentialDistortionAnalysis|doIndicatrices)
		warning('Plots are used to show Differential Distortion Analysis. Nothing is plotted.')
	else
	    disp(['Plotting'])
	    tic
	    % Set plot options
	    qp = 0.01; % Cumulative probability for quantiles, in [0,0.5]. Used for soft colorbar range. Set value to 0 to work with min and max.
	    min_log2sigma = quantile(aps_c_log2sigma,qp); max_log2sigma = quantile(aps_c_log2sigma,1-qp);
	    maxDiff_log2sigma = max(abs(min_log2sigma),max_log2sigma);
	    cmin_log2sigma = -maxDiff_log2sigma; cmax_log2sigma = maxDiff_log2sigma; % colorbar centered around 0
	    cmin_2Omega = 0; cmax_2Omega = quantile(aps_c_2Omega,1-qp);
	    axisRange = [min(gcps_c(1,:))-spatBuffer_c max(gcps_c(1,:))+spatBuffer_c min(gcps_c(2,:))-spatBuffer_c max(gcps_c(2,:))+spatBuffer_c];
	 
	    % Plot 2Omega
	    clf; hold on;
	    axis(axisRange); axis equal;
	    set(gca,'color','none');
	    surf(mps_xgrid_c, mps_ygrid_c, zeros(size(mps_xgrid_c)),reshape(aps_c_2Omega(1,1:nMps),size(mps_xgrid_c)),'EdgeColor','none')
	    colormap(whitegreen);
	    caxis([cmin_2Omega cmax_2Omega]);
	    % caxis([0 1]);
	    colorbar;
	    % fnplt(tps); % Other option. Only works with Matlab
	    plot(gcps_d_to_c(1,:),gcps_d_to_c(2,:),'.k'); % Not displaying in Octave?
	    plot(gcps_c(1,:),gcps_c(2,:),'.y'); % Not displaying in Octave?
	    hold off;
	    if isOctave
	        print('output/2Omega.png')
	    else
	        savefig('output/2Omega.fig')
	    end

	    % Plot log2sigma
	    clf; hold on;
	    axis(axisRange); axis equal;
	    set(gca,'color','none');
	    surf(mps_xgrid_c, mps_ygrid_c, zeros(size(mps_xgrid_c)),reshape(aps_c_log2sigma(1,1:nMps),size(mps_xgrid_c)),'EdgeColor','none')
	    colormap(redblue);
	    caxis([cmin_log2sigma cmax_log2sigma]);
	    % caxis([-1 1]);
	    colorbar;
	    % fnplt(tps); % Other option. Only works with Matlab
	    plot(gcps_d_to_c(1,:),gcps_d_to_c(2,:),'.k'); % Not displaying in Octave?
	    plot(gcps_c(1,:),gcps_c(2,:),'.y'); % Not displaying in Octave?
	    hold off;
	    if isOctave
	        print('output/log2sigma.png')
	    else
	        savefig('output/log2sigma.fig')
	    end
	    toc
	end
end