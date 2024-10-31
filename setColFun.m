function m = setColFun(m) % functions, other than task functions, for setCol

	% Define function handles for functions in this file
	m.array = @array;
	m.calRog = @calRog;
	m.cone = @cone;
	m.dog = @dog;
	m.dogS = @dogS;
	m.radHorComb = @radHorComb;
	m.radHorConv = @radHorConv;
	m.radHorLim = @radHorLim;
	m.readFile = @readFile;
	m.rogPar = @rogPar;
	m.rog = @rog;

function loc = array(sep, wid) % calculate triangular array

%	Inputs:
%		sep = distance between nearest-neighbour nodes (deg)
%		wid = array width (deg)
% Output:
%		loc = node locations, [x, y] (deg): ls x 2, where ls is the number of nodes
% Method:
%		The code calculates a parallelogram of nodes and then crops it to a square.

	w = .5 * wid; % half-width (deg)
	u = exp(1i * pi / 3); % unit vector along oblique axis
	i = ceil((1 + real(u)) * w / sep); % number of nodes to right of centre
	x = linspace(- i, i, 2 * i + 1) * sep; % hor. dist. of nodes from centre (deg)
	i = ceil(w / (imag(u) * sep)); % number of nodes above centre
	y = linspace(- i, i, 2 * i + 1) * sep; % ver. dist. of nodes from centre (deg)
	loc = x' + y * u; % generate array of nodes
	loc = [real(loc(:)), imag(loc(:))]; % locations (deg): ls x 2, ls too big
	i = abs(loc); i = i(:, 1) <= w & i(:, 2) <= w; % nodes within square
	loc = loc(i, :); % locations of nodes within square (deg): ls x 2

function r = calRog(b, x) % calculate ratio of Gaussians model

	% Set model parameters
	kC = b(1); % magnitude of centre Gaussian
	kS = b(2); % ratio of surround to centre magnitude
	rC = b(3); % radius of centre Gaussian (deg)
	rS = b(4); % radius of surround Gaussian (deg)
	tau = b(5); % time constant (s)
	disp = b(6); % stimulus displacement from receptive field middle (deg)
	del = b(7); % transport delay (s)
	kA = b(8); % magnitude of adaptation
	if length(b) == 8 % set number of processing stages
		ss = 3; % cone, bipolar cell, ganglion cell
	else
		ss = b(9); % specified by user
	end

	% Convert predictor variables
	freqS = 2 * pi * x(:, 1); % spatial frequency (radians/deg): fs x 1
	freqT = 2 * pi * x(:, 2); % temporal frequency (radians/s): fs x 1

	% Calculate basic model
	gC = exp(-.25 * (rC * freqS) .^ 2); % centre Gaussian: fs x 1
	gS = exp(-.25 * (rS * freqS) .^ 2); % surround Gaussian: fs x 1
	c = 1 ./ (1 + 1i * tau * freqT); % cone system function: fs x 1
	r = kC * (c .^ ss) .* gC ./ (1 + kS * (c .^ 2) .* gS); % response: fs x 1

	%	Add extra factors required by empirical data
	pDisp = exp(1i * disp * freqS); % phase shift due to stim. disp. (radians)
	pDel = exp(- 1i * del * freqT); % phase shift due to transport delay (radians)
	a = (1 + kA * 1i * tau * freqT) ./ (kA * (1 + 1i * tau * freqT)); % adaptation
	r = pDisp .* pDel .* a .* r; % add factors

function s = cone(b, x) % calculate cone contrast sensitivity

	f = 2 * pi * x; % temporal frequency (radians/s)
	s0 = b(1); % maximum contrast sensitivity
	tau = b(2); % time constant (s)
	s = abs(s0 ./ (1 + 1i * tau * f)); % contrast sensitivity

function r = dog(b, x, bFix) % calculate difference of Gaussians model

	% Initialise
	i = isnan(bFix); % indices of coefficients for optimisation
	bFix(i) = b; b = bFix; % insert those coefficients
	fs = .5 * length(x); % number of frequencies
	x = x(1: fs, :); % second half of frequencies is redundant: fs x 2
	freqS = 2 * pi * x(:, 1); % spatial frequency (radians/deg): fs x 1

	% Set model parameters
	kC = b(1); % magnitude of centre Gaussian
	kS = b(2); % ratio of surround to centre magnitude
	rC = b(3); % radius of centre Gaussian (deg)
	rS = b(4); % radius of surround Gaussian (deg)
	pC = b(5); % phase of centre Gaussian (deg)
	pS = b(6); % phase of surround Gaussian (deg)
	disp = b(7); % stimulus displacement from receptive field middle (deg)

	% Calculate model
	pC = pi * pC / 180; % phase of centre Gaussian (radians)
	pC = exp(1i * pC); % centre phase shift (radians)
	gC = pC * exp(-.25 * (rC * freqS) .^ 2); % normalised centre Gaussian: fs x 1
	pS = pi * pS / 180; % phase of surround Gaussian (radians)
	pS = exp(1i * pS); % centre phase shift (radians)
	gS = pS * exp(-.25 * (rS * freqS) .^ 2); % normalised centre Gaussian: fs x 1
	r = kC * (gC + kS * gS ); % response: fs x 1

	%	Add extra factors required by empirical data
	pDisp = exp(1i * disp * freqS); % phase shift due to stim. disp. (radians)
	r = pDisp .* r; % add factors

	% Store real and imaginary parts
	rR = real(r); % real part of response: fs x 1
	rI = imag(r); % imaginary part of response: fs x 1
	r = [rR; rI]; % combine: 2 * fs x 1

function r = dogS(b, x, bFix) % calculate DOG model in spatial domain

	% Initialise
	i = isnan(bFix); % indices of coefficients for optimisation
	bFix(i) = b; b = bFix; % insert those coefficients
	xs = .5 * length(x); % number of locations
	x = x(1: xs); % second half of location vector is redundant: xs x 1

	% Set model parameters
	kC = b(1); % magnitude of centre Gaussian
	kS = b(2); % ratio of surround to centre magnitude
	rC = b(3); % radius of centre Gaussian (deg)
	rS = b(4); % radius of surround Gaussian (deg)
	pC = b(5); % phase of centre Gaussian (deg)
	pS = b(6); % phase of surround Gaussian (deg)

	% Calculate model
	pC = pi * pC / 180; % phase of centre Gaussian (radians)
	pC = exp(1i * pC); % centre phase shift (radians)
	gC = pC * exp(- (x / rC) .^ 2); % complex centre Gaussian: xs x 1
	pS = pi * pS / 180; % phase of surround Gaussian (radians)
	pS = exp(1i * pS); % surround phase shift (radians)
	gS = pS * exp(- (x / rS) .^ 2); % complex surround Gaussian: xs x 1
	r = kC * (gC + kS * gS); % response: xs x 1

	% Store real and imaginary parts
	rR = real(r); % real part of response: xs x 1
	rI = imag(r); % imaginary part of response: xs x 1
	r = [rR; rI]; % combine: 2 * xs x 1

function d = radHorComb(d, m) % combine horizontal cell field radii: Packer (02)

	% Split the data into eccentricity bins
	i = d.source == 'wide' & d.eccFun >= 16 & d.eccFun < 60; % relevant rows
	dC = d(i, :); % select wide fields and eccentricities with sufficient data
	bins = 8; % number of eccentricity bins
	[~, edge, bin] = histcounts(dC.eccFun, bins); % ecc. bin edges and indices
	eccCen = edge(1: bins) + .5 * (edge(2) - edge(1)); % centres of ecc. bins

	%	Split each eccentricity bin into radius bins, calc. radius for each bin
	binRads = 4; % number of radius bins: each cone is contacted by 4 horizontal
		% cells, Wässle (89)
	r = zeros(bins, binRads); % radius: (eccentricity bin) x (radius bin)
	for i = 1: bins % loop over eccentricity bins
		dCC = dC(bin == i, :); % rows for this bin
		rC = dCC.radius; % radii in this bin (deg)
		[~, ~, binRad] = histcounts(rC, binRads); % bin number for radii
		for j = 1: binRads % loop over radius bins
			r(i, j) = mean(rC(binRad == j)); % mean of radii in current bin (deg)			
		end
	end
	
	% Calculate radii
	switch m.radHor.method % choose method
		case 'conv' % calculate radius of horizontal-cone convergence function
			dR = radHorConv(eccCen, m); % table containing radius
			rad = dR.radius; % convergence function radius (deg): es x 1
		case 'grid' % count grid points
			%	*** This code is untested ***
			%	Create grid of points,
			% randomise horizontal cell locations, find which grid points lie
			% within union, divide number of union points by number of grid
			% points, multiply grid area by fraction

			% Loop over eccentricity bins
			for i = 1: bins
		
				% Create grid
				rMax = max(rC(i, :)); % maximum radius at this eccentricity (deg)
				x = rMax * (- 2: .1: 2); % grid x locations (deg)
				[x, y] = ndgrid(x); % grid points (deg)
				locG = [x(:), y(:)]; % grid location (deg): gs x 2
				locGs = size(locG, 1); % number of grid locations
				in = false(1, locGs); % grid point is in receptive field: 1 x gs
				rUnion = zeros(1, bins); % allocate storage
				for j = 1: binRads % loop over radius bins
		
					% Calculate receptive field locations
					mag = rC(i, j); % current field is one radius from cone (deg)
					ang = 2 * pi * rand; % and at a random angle (rad)
					locF = mag * exp(ang); % field location as vector (deg)
					for k = 1: locGs % loop over grid points
						in(k) = abs(locG - locF) <= mag;
					end
		
				end
		
				%	Calculate radius
				a = sum(in) / locGs; % union area (deg ^ 2)
				rUnion(i) = sqrt(a / pi); % radius (deg)
				disp(rUnion);
		
			end

		case 'sum' % sum squared radii
			rad = sqrt(sum(r .^ 2, 2, 'omitnan')); % sum-of-squares radii (deg)
		case 'union' % find area of union of polygons

			% Loop over eccentricity bins
			pInit = nsidedpoly(100); % polygon: circular r.f. with radius 1 deg
			rPoly = nan(1, bins); % storage for calculate radii (deg)
			for i = 1: bins

				% Loop over radius bins
				p = polyshape; % storage for receptive fields
				for j = 1: binRads % scale and translate each receptive field
					mag = r(i, j); % current radius (deg)
					if ~ isnan(mag)
						pC = scale(pInit, mag); % correct field radius (deg)
						ang = - pi + .5 * pi * (j - 1); % angle spread around origin (rad)
						loc = mag * exp(1i * ang); % field location as complex vector (deg)
						pC = translate(pC, [real(loc), imag(loc)]);
							% translate so that circumference passes through origin
						p(j) = pC;
					end
				end

				% Calculate radius for this eccentricity bin
				p = union(p);
				%	figure('windowStyle', 'docked'); plot(p);
				rPoly(i) = sqrt(area(p) / pi); % calculate radius

			end
			rad = rPoly'; % radius (deg): 1 x bins

	end

	% Store
	dC = repmat(dC(1, :), [bins, 1]); % one row for each bin
	dC.eccFun = eccCen'; % ecc. bin centres (deg)
	dC.source(:) = 'field'; % method for calculating radius
	dC.radius = rad; % radius (deg)
	d = [d; dC]; % concatenate

function d = radHorConv(ecc, m) % radius of hor.-cone converg'e func.

	% Initialise
	xs = 101; % number of x values
	if isfield(m.radHor, 'iter') % number of iterations for each eccentricity
		is = m.radHor.iter; % user-specified
	else
		is = 1; % default
	end
	if isfield(m.radHor, 'seed') % seed for randomisation
		seed = m.radHor.seed; % user-specified
	else
		seed = 0; % default
	end

	%	Set lower and upper limits for horizontal cell field diameter
	ecc = ecc'; % eccentricity at which to calculate radius (deg): es x 1
	eccMm = (1 / m.p.magRet) * ecc; % eccentricity (mm): es x 1
	d = radHorLim(m); % read diameter limits
	eccL = d.lowerMm(1, :); diamL = d.lowerMm(2, :); % diam. lower limits: 1 x es
	diamL = interp1(eccL, diamL, eccMm); % diameters at required ecc. (mm): es x 1
	eccU = d.upperMm(1, :); diamU = d.upperMm(2, :); % diam. upper limits: 1 x es
	diamU = interp1(eccU, diamU, eccMm); % diameters at required ecc. (mm): es x 1
	diamMm = [diamL, diamU]; % limits of h.c. receptive field diam. (mm): es x 2

	%	Read horizontal cell diameter
	folder = [m.folder, filesep, 'Wässle (89) Hor density']; % data folder
	load([folder, filesep, 'Density'], 'd'); % create table from folder
	densMm = interp1(d.x, d.y, eccMm); % densities are required ecc. (mm); es x 1
	d = table(eccMm, diamMm, densMm); % create table: es x 3

	%	Convert to degrees
	d.ecc = ecc; % eccentricity (deg): es x 1
	d.diam = m.p.magRet * diamMm; % diameter (deg): es x 2
	d.dens = m.p.magRet ^ -2 * densMm; % density (cells/deg^2): es x 1
	d.Properties.VariableDescriptions = ...
		{'Eccentricity (mm)', 'Diameter (mm)', 'Density (cells/mm^2)', ...
		'Eccentricity (deg)', 'Diameter (deg)', 'Density (cells/deg^2)'};

	%	Loop over eccentricities
	es = length(ecc); %	number of eccentricities
	d.x = zeros(es, xs); % x values (deg): es x xs
	d.loc = cell(es, 1); % cell locations (deg): es x 1 cell
	d.field = cell(es, is); % receptive field profiles: es x is cell
	d.fieldAtt = cell(es, is); % attenuated rec. field profiles: es x is cell
	d.conv = zeros(es, xs, is); % h.c.-cone converg'e fun. profile: es x xs x is
	for j = 1: es % loop over eccentricities

		%	Create spatial grid
		diam = d.diam(j, :); % current diameter (deg): 1 x 2
		if isfield(m.radHor, 'width') % visual field width (deg)
			wid = m.radHor.width; % user-specified
		else
			wid = diam(2); % default
		end
		x = .5 * wid * linspace(-1, 1, xs); % x values (deg): 1 x xs
		[xG, yG] = ndgrid(x); % grid locations (deg): xs x ys
		d.x(j, :) = x; % store x values (deg): 1 x xs
	
		%	Create horizontal cell array
		s = 1 / sqrt(sin(pi / 3) * d.dens(j)); % cell separation (deg)
		loc = array(s, wid); % calculate the cell locations (x, y) (deg): ls x 2
		locS = shiftdim(loc, -2); % prepare for calculation of r.f.: 1 x 1 x ls x 2
		locX = locS(1, 1, :, 1); locY = locS(1, 1, :, 2); % x, y (deg): 1 x 1 x ls
		disp = sqrt((xG - locX) .^ 2 + (yG - locY) .^ 2);
			% displacement of (x, y) from centre of rec. field (deg): xs x ys x ls
		d.loc(j) = {loc}; % store locations (deg): 1 x 1 cell
	
		%	Loop over random values
		rng(seed); % fix seed for reproducibility
		for k = 1: is % loop over iterations

			%	Calculate horizontal cell receptive fields
			ls = size(loc, 1); % number of horizontal cells
			i = rand([ls, 1]); % uniformly dist. indices, one for each h.c.: ls x 1
			diamC = diam(1) + (diam(2) - diam(1)) * i; % r.f. diam. (deg): ls x 1
			diamC = shiftdim(diamC, -2); % prepare for calculation of r.f.: 1 x 1 x ls
			rad = diamC / (2 * log(10)); % radius (deg): 1 x 1 x ls
			r = exp(- disp ./ rad); % receptive fields: xs x ys x ls
			d.field(j, k) = {r}; % store receptive fields: 1 x 1 cell
		
			%	Calculate attenuated receptive fields
			i = x == 0; % index of cone location: 1 x xs
			a = r(i, i, :); % distance-based attenuation: 1 x 1 x ls
			r = a .* r; % attenuate signal from distant horizontal cells: xs x ys x ls
			d.fieldAtt(j, k) = {r}; % store attenuated receptive fields: 1 x 1 cell
	
			%	Calculate convergence function of horizontal cells onto cone
			r = sum(r, 3); % sum over r.f to make convergence fun.: xs x ys
			r = sum(r, 2); % sum over ys to make function of x: xs x 1
			r = r / max(r); % normalise function: xs x 1
			d.conv(j, :, k) = r; % store convergence function: 1 x xs
		
			%	Estimate Gaussian radius
			%{
			fun = @(rad, x) exp(- x .^ 2 / rad ^ 2); % Gaussian function to be fitted
			model = fitnlm(x', r, fun, .8); % fit Gaussian
			d.radius(j, k) = abs(model.Coefficients.Estimate); % Gaussian radius (deg)
			%}

		end % loop over iterations

		%	Estimate Gaussian radius
		xC = repmat(x', [is, 1]); rC = d.conv(j, :)';
			% concatenation over iterations: is * xs x 1
		fun = @(rad, x) exp(- x .^ 2 / rad ^ 2); % Gaussian function to be fitted
		model = fitnlm(xC, rC, fun, .8); % fit Gaussian
		d.radius(j) = abs(model.Coefficients.Estimate); % Gaussian radius (deg)

	end % loop over eccentricities

function d = radHorLim(m) % set lower and upper limits for hor. cell diameter

	%	Lower limit
	eccL = [0, 4.03, 16]; % breakpoints for lower limit of h.c. r.f. (mm): 1 x bs
	diamL = .001 * [0, 67, 67]; % r.f. diameter at breakpoints (mm): 1 x bs
	lowerMm = [eccL; diamL]; % lower limit (mm): 2 x bs
	lower = m.p.magRet * lowerMm; % lower limit (deg): 2 x bs
	radLower = lower; % convert diameter to radius
	radLower(2, :) = lower(2, :) / (2 * log(10)); % radius (deg): 1 x bs

	%	Upper limit
	eccU = [0, 4.08, 6.36, 8.02, 10.7, 16]; % breakpoints for upper limit
	diamU = .001 * [0, 70.9, 171, 472, 673, 673]; % r.f. diameter at breakpoints
	upperMm = [eccU; diamU]; % upper limit (mm): 2 x bs
	upper = m.p.magRet * upperMm; % upper limit (deg): 2 x bs
	radUpper = upper; % convert diameter to radius
	radUpper(2, :) = upper(2, :) / (2 * log(10)); % radius (deg): 1 x bs

	%	Store
	d = table(lowerMm, upperMm, lower, upper, radLower, radUpper); % create table

function d = readFile(folder) % read raw data in folder and return table

	% Read contents file
	d = readtable([folder, '/', 'Contents'], 'textType', 'string', ...
		'delimiter', 'tab', 'consecutiveDelimitersRule', 'join'); % contents table
	f = d.file; % names of files
	c = d.content; % names of input variables
	l = d.label; % labels of output data
	v = split(d.var(1), ', ', 2); % names of output variables
	desc = split(d.desc(1), ', ', 2); % descriptions of output variables
	
	% Complile a table from each mat-file
	fs = length(f); % number of files
	dC = cell(fs, 1); % allocate storage
	for i = 1: fs % loop over files
		fC = f(i); % current file
		if contains(fC, ".csv") % comma-separated variable file
			d = readtable([folder, filesep] + fC); % file contents as table
			cs = size(d, 1); lC = repmat(l(i), [cs, 1]); % contents label
			cC = [d.Var1, d.Var2]; % make it double
		else % mat-file
			fC = f(i) + ".mat"; % name of current file
			s = load([folder, filesep] + fC); % file contents: structure
			cC = c(i); % current content
			cC = s.(cC); % file contents: double
			cs = size(cC, 1); % number of rows
			lC = repmat(l(i), [cs, 1]); % contents label
		end
		dC{i} = table(lC, cC(:, 1), cC(:, 2), 'variableNames', ["label", v]);
	end

	% Combine and store
	d = vertcat(dC{:}); % combine data tables
	d.Properties.VariableDescriptions = ["label", desc]; % describe variables

function r = rog(b, x, bFix, pulse) % initialise, execute, finalise ROG model

	%	Combine fixed and variable coefficients
	i = isnan(bFix); % indices of variable coefficients
	bFix(i) = b; % insert those coefficients: 1 x bs
	b = bFix; % store

	%	Prepare predictor variables and calculate ROG
	if nargin == 3 % frequency response
		cs = 2; % number of copies of frequencies: one each for real, imag. resp.
	else % pulse response
		cs = size(pulse, 2); % one copy for each response
	end
	fs = length(x) / cs; % number of uncopied frequencies
	x = x(1: fs, :); % uncopied frequencies: fs x 2
	r = calRog(b, x); % calculate ROG: fs x 1

	%	Finalise
	if nargin == 3 % frequency response: store real and imaginary parts
		rR = real(r); % real part of response: fs x 1
		rI = imag(r); % imaginary part of response: fs x 1
		r = [rR; rI]; % combine: 2 * fs x 1
	else % pulse response
		r = r .* pulse; % transform of pulse response: fs x rs
		r = (fs ^ 2) * r; % change to fft amplitude: fs x rs
		r = ifft(r); % inverse transform: fs x rs
		r = real(r); % model must return real value: fs x rs
		%	warn about significant imaginary components
		r = r(:); % linearise: fs * rs x 1
	end

function b = rogPar(m) % calculate ROG parameters from full model parameters

	%	Initialise
	e = m.p.ecc; % eccentricity (deg)
	radOpt = polyval(m.p.radOptCoef, e); % radius of point spread function (deg)
	radGang = polyval(m.p.radGangCoef, e); % radius of gang. cell dendrites (deg)
	c =  m.p.radBackCoef; % regression coefficients for feedback convergence fun.
	radBack = 10 ^ (c(1) + c(2) * normcdf(e, c(3), c(4))); % feedback c.f. (deg)

	%	Set ROG coefficients
	b(1) = m.p.kRect * m.p.kSens; % kCen, gang. cell cont. sens'y (Hz/cont.-unit)
	b(2) = m.p.kSur; % kSur, ratio of surround to centre magnitude
	b(3) = sqrt(radOpt ^ 2 + radGang ^ 2); % centre radius (deg)
	b(4) = sqrt(radOpt ^ 2 + radBack ^ 2); % surround radius (deg)
	b(5) = m.p.tau; % time constant (s)
	b(6) = 0; % stimulus displacement from receptive field middle (deg)
	b(7) = 0; % transport delay (s)
	b(8) = 1; % magnitude of adaptation
