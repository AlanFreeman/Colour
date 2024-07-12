function m = runColFun(m) % calculate model responses

%	The responses include the time course of generator potentials.

	% Define function handles for functions in this file
	m.calConv = @calConv;
	m.calStim = @calStim;
	m.calStruct = @calStruct;
	m.calTemp = @calTemp;
	m.doCrop = @doCrop;
	m.getVal = @getVal;
	m.ifftReal = @ifftReal;
	m.listLoc = @listLoc;
	m.readFile = @readFile;
	m.rect = @rect;
	m.setType = @setType;
	m.setVal = @setVal;
	m.solveF = @solveF;
	m.solveT = @solveT;

function loc = array(wid, sep, off) % calculate triangular array

%	Input:
%		wid = visual field width (deg)
%		sep = distance between nearest-neighbour nodes (deg)
%	Optional input:
%		off = offset, [x, y] (deg)
% Output:
%		loc = node locations, [x, y] (deg): ls x 2, where ls is the number of nodes
% Method:
%		The code calculates a parallelogram of nodes and then crops it to a square.
%		The parallelogram is wide enough to fill the visual field when cropped.

	%	Generate a parallelogram of nodes
	if ~ exist('off', 'var'), off = [0, 0]; end % default offset
	w = .5 * wid; % half-width (deg)
	u = exp(1i * pi / 3); % unit vector along oblique axis
	d = (1 + real(u)) * w + off(1); % distance to right of centre (deg)
	i = ceil(d / sep); % number of nodes to right of centre
	x = linspace(- i, i, 2 * i + 1) * sep; % hor. dist. of nodes from centre (deg)
	d = w + off(2); % distance above centre (deg)
	i = ceil(d / (imag(u) * sep)); % number of nodes above centre
	y = linspace(- i, i, 2 * i + 1) * sep; % ver. dist. of nodes from centre (deg)
	loc = x' + y * u; % generate array of nodes
	
	% Turn the nodes into locations and trim to a square
	loc = [real(loc(:)), imag(loc(:))]; % locations (deg): ls x 2, ls too big
	loc = loc + off; % include offset
	i = abs(loc); i = i(:, 1) <= w & i(:, 2) <= w; % nodes within square
	loc = loc(i, :); % locations of nodes within square (deg): ls x 2

function g = calConv(locI, locO, rad) % calculate convergence matrix

%	Input:
%		locI = locations (x, y) of input neurons (deg): is x 2
%		locO = locations (x, y) of output neurons (deg): os x 2
%		rad = radius of Gaussian attenuation (deg): 1 x 1
%	Output:
%		g = Gaussian function of distance between locI and locO (deg): os x is

	locI = permute(locI, [3, 1, 2]); % prepare for subtraction: 1 x is x 2
	locO = permute(locO, [1, 3, 2]); % os x 1 x 2
	g = locO - locI; % displacement of output neuron from input: os x is x 2
	g = sum(g .^ 2, 3); % squared distance of output neuron from input: os x is
	g = exp(- g / rad ^ 2); % Gaussian function of distance: os x is

function dpot = calDeriv(t, pot, m) % calc. derivatives of model equations

	% Sort potential into cell types
	i = 0; % current row in potential vector
	for c = m.p.cell % loop over cell types
		type = c.type; % current type
		a = c.array; % location array
		l = m.p.(a).loc; ls = size(l, 1); % locations, number of locations
		p.(type) = pot(i + (1: ls)); % generator potential (mV): ls x 1
		i = i + ls; % update row index
	end

	% Calculate derivatives for cone array
	rate = 1 / m.p.tau; % reciprocal of time constant (s^-1)
	drive = calDrive(t, m)'; % cone generator potential (mV): ls x 1
	back = m.p.kSur * m.p.wHorCone * p.hor; % feedback potential (mV): ls x 1
	dp.cone = rate * (- m.p.kSens * drive - m.p.back * back - p.cone);
		% cone potential (mV/s): ls x 1
	dp.hor = rate * (m.p.wConeHor * p.cone - p.hor); % hor. cell (mV/s): ls x 1
	dp.back = 0 * dp.hor; % feedback (mV/s): ls x 1
	dp.biOff = rate * (p.cone - p.biOff); % off-bipolar (mV/s): ls x 1
	dp.biOn = rate * (- p.cone - p.biOn); % on-bipolar (mV/s): ls x 1
	
	% Calculate derivatives for ganglion cell arrays
	for i = ["Off", "On"] % loop over centre signs
		w = "wBiGang" + i; bi = "bi" + i; gang = "gang" + i;
			% construct names
		dp.(gang) = rate * (m.p.(w) * p.(bi) + m.p.potRest - p.(gang));
			% ganglion cell potential (mV): ls x 1
	end

	%	Vectorise derivatives
	dpot = zeros(size(pot)); % initialise derivatives
	i = 0; % current row in potential vector
	for c = m.p.cell % loop over cell types
		type = c.type; % current type
		dpC = dp.(type); ls = size(dpC, 1); % derivatives, number of derivatives
		dpot(i + (1: ls)) = dpC; % generator potential (mV): ls x 1
		i = i + ls; % update row index
	end

function drive = calDrive(t, m) % cross-correlate stim. and cone conv. function

	% Calculate stimulus and attenuate with the optical point spread function
	loc = m.p.cone.loc; % cone location (deg): ls x 2
	d = calStim(t, loc, m); % stimulus (contrast units): ts x ls x cs
	fS = 2 * pi * m.p.freqS; % spatial frequency (radians/deg)
	d = exp(-.25 * (m.p.radOpt * fS) ^ 2) .* d; % attenuated: ts x ls x cs

	% Choose cone contrast appropriate to each one
	[~, ls, cs] = size(d); % number of locations, cone types
	i = m.p.cone.type; % cone type: ls x 1
	i = sub2ind([ls, cs], (1: ls)', i); % convert to indices: ls x 1
	drive = d(:, i); % cone-specific stimulus: ts x ls

function s = calStim(t, loc, m) % calculate stimulus movie

% Inputs:
%		t:
%			for temporal domain: times (s): vector with length ts
%			for frequency domain: empty
%		loc:
%			for grating: visual field locations, (x, y) (deg): ls x 2
%			for image: empty
%		m = metadata
%	Ouput:
%		s(t, loc, c) (contrast-units): ts x ls x cs, ts = 1 for frequency domain

	% Initialise
	t = t(:); % make time a column vector: ts x 1
	dir = (pi / 180) * m.p.dir; % motion dir. (radians ACW from rightward)
	fS = 2 * pi * m.p.freqS; % spatial frequency (radians/deg)
	fT = 2 * pi * m.p.freqT; % temporal frequency (radians/s)
	pS = (pi / 180) * m.p.phaseS; % spatial phase (radians)

	% Transform location to distance in direction of motion
	u = loc * [cos(dir); sin(dir)]; % dist. relative to centre (deg): ls x 1
	u = u'; % make it a row vector (deg): 1 x ls
	u = u - m.p.locS; % distance relative to stimulus location (deg): 1 x ls

	% Calculate stimulus
	if isempty(t) % frequency domain *** does this work for flash? ***
		s = exp(- 1i * (fS * u - pS)); % phase shift at cone location: 1 x ls
	else % time domain
		%	s = zeros(length(t), size(loc, 1)); % zero stimulus: ts x ls
		switch m.p.stim % type of grating
			case 'drift' % drifting grating
				s = cos(fS * u - pS - fT * t); % stimulus: ts x ls
			case 'flash' % stationary grating
				s = cos(fS * u - pS) * (t <= m.p.dur); % stimulus: ts x ls
					% *** should be .* ***
		end
	end

	% Multiply by contrast
	c = m.p.cont; % contrast: 1 x cs
	c = shiftdim(c, -1); % prepare to multiply: 1 x 1 x cs
	s = c .* s; % multiply: ts x ls x cs

function m = calStruct(m) % calculate the model structure

	% Initialise
	rng(m.p.seed); % make the results reproducible
	e = m.p.ecc; % current eccentricity (deg)
	
	% Calculate the densities and convergence function radii
	m.p.radOpt = polyval(m.p.radOptCoef, e); % point spread function (deg)
	m.p.radHor = exp(polyval(m.p.radHorCoef, e)); % hor. cell rec. field (deg)
	m.p.radGang = polyval(m.p.radGangCoef, e); % ganglion cell dendritic (deg)
	m.p.radCen = sqrt(m.p.radOpt ^ 2 + m.p.radGang ^ 2); % centre radius (deg)
	m.p.radSur = sqrt(m.p.radOpt ^ 2 + 2 * m.p.radHor ^ 2); % surround rad. (deg)

	% Calculate the arrays and convergence functions
	m = cone(m); % cone array
	m = gang(m); % ganglion cell array
	m = weight(m); % convergence functions

function m = calTemp(m) % calculate the sample times and temporal frequencies

	% Calculate the number of samples and fundamental frequency
	n = m.p.ts; % number of sample times and temporal frequencies
	nH = floor(.5 * n); % half the number, converted to integer
	n = 2 * nH; % make sure that the number of samples is even
	switch m.p.solver % which domain?
		case 'solveF' % frequency domain
			fund = m.p.freqT; % fundamental frequency (Hz)
		case 'solveT' % time domain
			fund = 1 / m.p.time; % fundamental frequency (Hz)
	end

	% Calculate the sample frequencies
	ind = linspace(- nH, nH, n + 1); % sample indices
	ind(end) = []; % make the sequence open-ended
	f = fund * ind; % sample frequencies (Hz): n x 1
	m.p.f = fftshift(f); % begin with zero frequency, for compatibility with fft

	% Calculate the sample times
	ind = linspace(0, 1, n + 1); % sample indices
	ind(end) = []; % make the sequence open-ended
	t = (1 / fund) * ind; % sample times (s): n x 1
	m.p.t = t; % begin with time 0, for intuitive use

function m = cone(m) % calculate the cone array

	% Calculate the cone locations
	load([m.folder, filesep, 'Cone density'], 'd'); % load cone density vs ecc.
	i = d.quad == "temporal"; d = d(i, :); % temporal retina only
	ecc = m.p.magRet * d.ecc; % eccentricity (deg)
	dens = d.dens / m.p.magRet ^ 2; % cone density (deg ^ -2)
	dens = interp1(ecc, dens, m.p.ecc); % density at specified ecc. (deg ^ -2)
	s = 1 / sqrt(sin(pi / 3) * dens); % cone separation (deg)
	loc = array(m.p.wid, s); % calculate the cone locations (x, y) (deg): cs x 2
	m.p.cone.loc = loc; % cone locations (x, y) (deg): cs x 2
		
	% Assign L and M cone types to all locations, ignoring S cones
	cs = size(loc, 1); % number of cones
	type = ones(cs, 1); % set all cones to L type: cs x 1
	ts = cs * m.p.ratCone; % number of L, M, S cones
	ls = ts(1); ms = ts(2); ss = ts(3); % number of cones of each type
	ms = round((ms / (ls + ms)) * cs); % number of M cones
	i = randperm(cs); % randomise cone order: 1 x cs
	type(i(1: ms)) = 2; % set fraction of cones to M type: cs x 1
	
	% Make an offset triangular array for the S cones and reassign these loc. to S
	dens = (ss / cs) * dens; % S cone density (cells/deg^2)
	s = 1 / sqrt(sin(pi / 3) * dens); % S cone separation (deg)
	off = .5 * s * tan(pi / 6) * [1, 1]; % S cone array offset, [x, y] (deg)
	locS = array(m.p.wid, s, off); % S cone array
	i = knnsearch(loc, locS); % find the nearest locations in the L, M array
	type(i) = 3; % set those locations to S type: cs x 1
	m.p.cone.type = type; % cone types: cs x 1
	
function [r, loc] = doCrop(r, loc, rad, dim) % crop border locations

%	Input:
%		r = response: multi-dimensional
%		loc = response locations: ls x 2
%		rad = half-width of square border in which to retain responses
%		dim = location dimension in r
%	Output:
%		r = cropped response
%		loc = retained locations

	s = size(r); % size of r
	es = prod(s(1: dim - 1)); % number of elements before dimension dim
	r = reshape(r, es, s(dim), []); % prepare for cropping
	i = abs(loc(:, 1)) <= rad & abs(loc(:, 2)) <= rad; % indices of loc. to keep
	r = r(:, i, :); % crop
	s(dim) = sum(i); % update size
	r = reshape(r, s); % restore original shape
	loc = loc(i, :); % update locations

function m = gang(m) % calculate the ganglion cell array

	% Load the data
	load([m.folder, filesep, 'Ganglion cell density'], 'd'); % midget density
	ecc = d.ecc; % functional eccentricity (deg): es x 1
	dens = d.densDeg; % midget ganglion cell density (deg ^ -2): es x 1
	rat = m.p.ratGang; % ratio of off-, on-midget to all ganglion cells
	rat = rat / sum(rat); % ratio of off-, on-centre cells to all midget g. cells
	r.off = rat(1); r.on = rat(2); % name the ratios

	% Calculate the ganglion cell locations at the current eccentricity
	for s = ["Off", "On"] % off- then on-centre
		sC = lower(s); % current centre type
		densC = r.(sC) * dens; % density of current centre type: (deg ^ -2): es x 1
		densC = interp1(ecc, densC, m.p.ecc); % density at specified ecc. (deg ^ -2)
		sep = sqrt(1 / (sin(pi / 3) * densC)); % cell separation (deg)
		loc = array(m.p.wid, sep); % calculate cell locations (x, y) (deg): ls x 2
		switch m.p.align % adjust locations
			case 1 % align ganglion cells with bipolar cells
				locCone = m.p.cone.loc; % cone locations (deg): cones x 2
				j = knnsearch(locCone, loc); % indices of nearest cones: cs x 1
				loc = locCone(j, :); % match location to that of nearest cone: cs x 2
			otherwise % randomise ganglion cell locations
				loc = loc + normrnd(0, m.p.kGangDev * sep, size(loc)); % perturb (deg)
		end
		n = "gang" + s; % name of current array
		m.p.(n).loc = loc; % ganglion cell location (x, y) (deg): ls x 2
		m.p.(n).type = sC; % ganglion cell type: 1 x 1
	end

function [name, val, vals] = getVal(m, task) % get the values of the model par.

% Inputs:
%		m.task.par, where par is a model parameter: matrix of parameter values,
%			with one row per value, e.g. m.task.cont = [1, 0, 0; 0, 1, 0]
%	Outputs:
%		name = names of specified model parameters
%		val = values for the specified model parameters: cell array with
%			one member for each model parameter
%		vals = number of values for each specified model parameter

	% Find the model parameters to be varied
	name = fieldnames(m.(task)); % names in m.task
	nameAll = fieldnames(m.p); % names of model parameters
	i = ismember(name, nameAll); % m.task fields that are model parameters
	name = name(i); % names of specified parameters
	names = length(name); % number of parameters
	
	% Find the values of the response parameters, which may be multi-column
	val = cell(1, names); % values of parameters
	vals = zeros(1, names); % number of values
	for i = 1: names % loop over response parameters
		nameC = name{i}; % name of current parameter
		valC = m.(task).(nameC); % values of parameter
		val{i} = valC; % store values
		vals(i) = size(valC, 1); % update number of values
	end

function r = ifftReal(r, dim) % calculate real part of inverse transform

	if nargin == 1 % the frequency dimension is unspecified
		dim = 1; % assume that the frequency dimension is 1
	end
	r = ifft(r, [], dim); % take inverse transform
	mReal = max(abs(rmoutliers(real(r(:))))); % maximum of real components
	mImag = max(abs(rmoutliers(imag(r(:))))); % maximum of imaginary components
	ratio = mImag / mReal; % ratio of imaginary to real maxima
	r = real(r); % make it real
	if ratio > 1e-5 % issue warning
		i = 'Col:LargeImaginary'; % warning ID
		w = 'Imaginary components of inverse transform are relatively large: %g';
		warning(i, w, ratio); % issue warning
		warning('off', i); % turn off further warnings with this ID
	end

function loc = listLoc(xs, m) % create a list of visual field locations

	x = .5 * m.p.wid * linspace(-1, 1, xs); % x values (deg): 1 x xs
	[x, y] = ndgrid(x); % grid locations (deg): xs x ys
	loc = [x(:), y(:)]; % locations (x, y) (deg): ls x 2

function [d, m] = readFile(m) % set data folder, read or create data file

	% Set names of data folder and data file
	m.folder = []; % default
	file = m.p.file; % data file name
	f = userpath; % top level of search path
	if contains(f, 'Alan', 'ignoreCase', 1) % restrict to Alan's computer
		m.folder = [f, filesep, 'Data', filesep, m.p.project]; % folder
		if ~ isempty(file) % there is a data file name
			file = [m.folder, filesep, file]; % data file with full path
		end
	else % anywhere else
		m.folder = cd; % current folder
	end
	
	% Read or create file
	if exist(file, 'file') % data file exists
		load(file, 'd'); % load data file
	else
		d = table; % default is empty data file
	end

function r = rect(r, dim) % rectify a frequency-domain response

	if nargin == 1 % the frequency dimension is unspecified
		dim = 1; % assume that the frequency dimension is 1
	end
	rC = ifftReal(r, dim); % transform to time domain and make it real
	if min(rC, [], 'all') < 0 % time domain response has negative elements
		rC = max(0, rC); % rectify
		r = fft(rC, [], dim); % return to frequency domain
	end

function m = setType(m) % create row vector of all cell types

	a = m.p.array; % cell array of neuronal array names
	c = cell(size(a)); % cell array to hold structure for each neuronal array
	for i = 1: length(a) % loop over neuronal arrays
		aC = a{i}; % name of current array
		c{i} = struct('type', m.p.(aC).stage, 'array', aC); % stuct. for this array
	end
	m.p.cell = horzcat(c{:}); % concatenate across structures

function m = setVal(m, name, val, iVal) % set the model parameter values

%	Inputs:
%		m = metadata
%		name = names of specified model parameters
%		val = values for the specified model parameters: cell array with
%			one member for each model parameter
%		iVal = linear index of current set of values

	names = length(name); % number of parameters
	s = arrayfun(@(x) size(x{1}, 1), val); % number of values for each parameter
	sub = cell(1, names); % value subscripts
	[sub{:}] = ind2sub(s, iVal); % subscripts for each parameter
	for i = 1: names % loop over parameters
		nameC = name{i}; % parameter's name
		valC = val{i}; % parameter's values
		subC = sub{i}; % subscript of the current value
		m.p.(nameC) = valC(subC, :); % set value for this model parameter
	end

function p = solveF(m) % solve the model equations in the frequency domain

	% Calculate system functions
	fT = 2 * pi * m.p.freqT; % temporal frequency (radians/s): 1 x 1
	aT = 1 / (1 + 1i * m.p.tau * fT); % temporal attenuation: 1 x 1
	c = aT; % cone system function: 1 x 1
	h = m.p.wConeHor * aT; % horizontal cell system function: ls x ls
	back = m.p.kSur * m.p.wHorCone; % feedback system function: ls x ls
	
	% Calculate cone signal
	d = calDrive([], m); % cross-corr. of stim., PSF (contrast units): 1 x ls
	d = - m.p.kSens * d .'; % drive (mV): ls x 1
	ls = length(d); % number of cones
	pC = (eye(ls) + m.p.back * c * back * h) \ (c * d); % cone resp. (mV): ls x 1
	
	%	Convert cone signal to fft format
	f = 2 * pi * m.p.f; % fft frequencies (radians/s): 1 x fs
	fs = length(f); % number of frequencies
	p.cone = zeros(ls, fs); % cone signal, all frequencies: ls x fs
	p.cone(:, 2) = pC; % fundamental component (mV): ls x fs
	p.cone(:, end) = conj(pC); % negative fundamental (mV): ls x fs
	p.cone = .5 * fs * p.cone; % fundamental fft component (mV): ls x fs
	
	% Calculate remaining cone array signals
	aT = 1 ./ (1 + 1i * m.p.tau * f); % temporal attenuation: 1 x fs
	p.hor = m.p.wConeHor * (aT .* p.cone); % horizontal cell signal (mV): ls x fs
	p.back = back * p.hor; % feedback signal (mV): ls x fs
	p.biOff = aT .* p.cone; % off-bipolar signal (mV): ls x fs
	p.biOn = - aT .* p.cone; % on-bipolar signal (mV): ls x fs
	
	% Calculate ganglion cell array signals
	for i = ["Off", "On"] % loop over centre signs
		w = "wBiGang" + i; bi = "bi" + i; gang = "gang" + i;
			% construct weight name
		p.(gang) = m.p.(w) * (aT .* p.(bi)); % ganglion cell fund. (mV): ls x fs
		p.(gang)(:, 1) = fs * m.p.potRest; % resting potential (mV): ls x fs
	end

	% Return stages, not arrays
	s = m.p.cell; % structure with one element for each stage
	for i = 1: length(s) % loop over stages
		sC = s(i).type; % current stage
		s(i).resp = p.(sC); % add response
	end
	p = s; % 1 x zs

function p = solveT(m) % solve the model equations in the time domain

	% Set initial values for potentials
	for c = m.p.cell % loop over stages
		s = c.type; % current stage
		a = c.array; % current location array
		l = m.p.(a).loc; ls = size(l, 1); % locations, number of locations
		p.(s) = zeros(ls, 1); % default resting potential (mV): ls x 1
	end
	pRest = m.p.potRest; % ganglion cell resting potential (mV)
	p.gangOff(:) = pRest; p.gangOn(:) = pRest; % ganglion cells
	pI = struct2cell(p); % make it a cell array
	pI = vertcat(pI{:}); % concatenate across stages: ls x 1, all locations

	% Numerically integrate the equations
	fun = @(t, p)calDeriv(t, p, m); % function to calculate derivatives
	t = m.p.t; % sampling times (s): 1 x ts
	[~, pot] = ode45(fun, t, pI); % potential (mV): ts x ls
	
	% Unpack potential
	i = 0; % current column in potential
	for s = m.p.cell % loop over stages
		sC = s.type; % current stage
		ls = length(p.(sC)); % number of locations
		pC = pot(:, i + (1: ls)); % store potentials in named array: ts x ls
		p.(sC) = pC'; % transpose to standard dimension order: ls x ts
		i = ls + i; % update column
	end

	% Add unintegrated signals
	p.back = m.p.kSur * m.p.wHorCone * p.hor; % f'back (mV):
		%	(ls x ls) x (ls x ts) = ls x ts

	% Rearrange potential to match old format
	s = m.p.cell; % structure with one element for each stage
	for i = 1: length(s) % loop over stages
		sC = s(i).type; % current stage
		s(i).resp = p.(sC); % add response: ls x ts
	end
	p = s; % structure, including responses: 1 x zs
	
function m = weight(m) % calculate convergence function synaptic weights

	% Calculate cone to horizontal attenuation
	locI = m.p.cone.loc; % cone locations: cs x 2
	locO = locI; % horizontal cells are colocated with cones: cs x 2
	r = m.p.radHor; % horizontal cell field radius: 1 x 1
	a = calConv(locI, locO, r); % distance-based attenuation: cs x cs
	
	% Set the weight to zero for absent inputs, and store
	i = m.p.cone.type == 3; % cone is S-type: cs x 1
	w = a; w(:, i) = 0; % H1 hor. cells receive little input from S-cones: cs x cs
	w = w ./ sum(w, 2); % sum of input weights is 1: cs x cs
	m.p.wConeHor = w; % cone to horizontal cell synaptic weights: cs x cs
	w = a'; w = w ./ sum(w, 2); % sum of input weights is 1: cs x cs
	m.p.wHorCone = w; % horizontal cell to cone synaptic weights: cs x cs

	% Calculate bipolar to off-ganglion cell attentuation
	locI = m.p.cone.loc; % bipolar cell locations: cs x 2
	locO = m.p.gangOff.loc; % ganglion cell locations: gs x 2
	r = m.p.radGang; % ganglion cell dendritic radius: 1 x 1
	w = calConv(locI, locO, r); % attenuation for off-cells: gFs x cs
	w = w ./ sum(w, 2); % sum of input weights is 1: gFs x bs
	m.p.wBiGangOff = w; % bipolar cell to off-gang. synaptic weights: gFs x cs
	
	% Calculate bipolar to on-ganglion cell attentuation	
	locO = m.p.gangOn.loc; % ganglion cell locations: gs x 2
	w = calConv(locI, locO, r); % attenuation for on-cells: gNs x cs
	w(:, i) = 0; % on-cen. g.c. don't receive input from S-bipolars: gNs x bs
	w = w ./ sum(w, 2); % sum of input weights is 1: gNs x bs
	m.p.wBiGangOn = w; % bipolar cell to on-gang. synaptic weights: gNs x bs
