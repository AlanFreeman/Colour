function m = runColTask(m) % execute an analysis task

	% Define function handles for functions in this file
	h = localfunctions; % handles of functions in this file
	hs = length(h); % number of handles
	for i = 1: hs % loop over handles
		hC = h{i}; % current handle
		task = func2str(hC); % name of task
		m.(task).fun = hC; % store handle in metadata
	end

function [d, m] = centre(d, m) % correct resp. phase for distance from origin

	% Find the model parameters to be varied
	[name, val, vals] = m.getVal(m, 'resp'); % model parameter names and values
	vs = prod(vals); % total number of values
	r = d.resp; % response: 1 x fs x ls x 1 x vs, vs can be multi-dimensional
	s = size(r); sC = [s(1: 4), vs]; % vs is 1-dimenional: 1 x fs x ls x 1 x vs
	r = reshape(r, sC); % reshape the response

	% Determine response phase shift for stimulus centred on cell
	for v = 1: vs % loop over variable values

		% Calculate distance of cell from origin
		m = m.setVal(m, name, val, v); % set values of model parameters
		%	m = m.calTemp(m); % recalculate temporal par. if temp. freq.	changes
		l = d.loc; % cell location (deg): 1 x ls x 2
		l = shiftdim(l, 1); % make it a matrix: ls x 2
		dir = pi * m.p.dir / 180; % grating direction (radians)
		freqS = 2 * pi * m.p.freqS; % grating spatial frequency (radians/deg)
		u = [cos(dir); sin(dir)]; % unit vector in grating direction: 2 x 1
		p = l * u; % distance of cell from origin in grating direction (deg): ls x 1

		% Correct phase
		p = freqS .* p; % phase shift from origin (radians): ls x 1
		p = exp(1i * p); % phase correction (complex): ls x 1
		p = shiftdim(p, -2); % prepare to correct phase: 1 x 1 x ls
		r(:, :, :, :, v) = p .* r(:, :, :, :, v); % correct phase: 1 x fs x ls x vs

	end

	% Store
	r = reshape(r, s); % restore response shape: vs can be multi-dimensional
	d.resp = r; % response with phase corrected: 1 x fs x ls x vs

function [d, m] = colour(d, m) % add colour map for specified neuronal arrays

	for a = string(m.colour.array) % loop over neuronal arrays
		aC = m.p.(a); % structure for neuronal array
		l = aC.loc; % cell locations (deg): ls x 2
		ls = size(l, 1); % number of cells
		switch a % set colours
			case 'cone' % cones
				c = eye(3); % RGB colour map for L, M, and S cones: 3 x 3
				j = aC.type; % cell type: ls x 1
				c = c(j, :); % marker colours: ls x 3
			case 'gangOff' % off-centre ganglion cells
				c = [0, 0, 1]; % blue for off-centre: 1 x 3
				c = repmat(c, [ls, 1]); % one row for each cell: ls x 3
			case 'gangOn' % on-centre ganglion cells
				c = [1, 0, 0]; % red for on-centre: 1 x 3
				c = repmat(c, [ls, 1]); % one row for each cell: ls x 3
		end
		cC.(a) = c; % store
	end
	c = struct2cell(cC); % convert to cell for concatenation
	d.colour = vertcat(c{:}); % combine arrays

function [d, m] = comp(d, m) % calc. cone components in centre and sur. resp.

	r = d.resp; % fundamental for all cone types (mV): 1 x 1 x ls x zs x cs
	r = abs(r); % response magnitude (mV): 1 x 1 x ls x zs x cs
	r = r ./ sum(r, 5); % normalise by sum over cone types: 1 x 1 x ls x zs x cs
	d.resp = r; % normalised cone components: 1 x 1 x ls x zs x cs
	d.Properties.VariableDescriptions{'resp'} = 'Cone ratio'; % describe

function [d, m] = crop(d, m) % crop locations to smaller field

	% Prepare response
	r = d.(m.z); % response: ps x ls x qs; ps and qs can be multi-dimensional
	loc = shiftdim(d.loc, 1); % locations (deg): ls x 2
	if isfield(m.crop, 'radius') % radius is specified
		rad = m.crop.radius; % radius (deg)
	else % default radius
		rad = m.p.radSur; % surround radius (deg)
	end
	rad = .5 * m.p.wid - rad; % half-width of cropped field (deg)
	if isfield(m.crop, 'dim') % set location dimension
		dim = m.crop.dim; % user specified
	else % default
		dim = d.Properties.CustomProperties.RespDim; % response dimensions
		[~, dim] = ismember('loc', dim); % index of location in response
	end

	%	Crop and store
	[r, loc] = m.doCrop(r, loc, rad, dim);
	d.(m.z) = r; % store response: ps x ls x qs, ps and qs can be multi-dim.
	d.loc = shiftdim(loc, -1); % locations (deg): ls x 2

function [d, m] = desc(d, m) % describe the variables in the data table

	desc = d.Properties.VariableDescriptions; % variable descriptions
	undesc = find(strcmp('', desc)); % indices of undescribed variables
	for i = undesc % loop over undescribed variables
		name = d.Properties.VariableNames{i}; % name of current variable
		switch name % describe this variable
			case 'cont', desc{i} = 'Contrast';	
			case 'dir', desc{i} = 'Stimulus direction (deg)';
			case 'freqS', desc{i} = 'Spatial frequency (cycles/deg)';
			case 'freqSPref', desc{i} = 'Preferred spatial frequency (cycles/deg)';
			case 'freqT', desc{i} = 'Temporal frequency (Hz)';
			case 'loc', desc{i} = 'Visual field location (deg)';
			case 'resp', desc{i} = 'Generator potential (mV)';
			case 'stage', desc{i} = 'Processing stage';
			case 'time', desc{i} = 'Time (s)';
			case 'x', desc{i} = 'Horizontal location (deg)';
			case 'xs', desc{i} = 'Number of locations';
			case 'y', desc{i} = 'Vertical location (deg)';
		end
	end
	d.Properties.VariableDescriptions = desc; % store

function [d, m] = doMax(d, m) % find max. resp. and correponding variable values

	% Maximise over time if it's a time course
	z = d.(m.z); % fund. resp.: 1 x ts x ls x 1 x ds x fs
	dim = d.Properties.CustomProperties.RespDim; % response dimensions
	if ismember('time', dim) % it's a time course
		z = max(z, [], 2); % maximise over time: 1 x 1 x ls x 1 x ds x fs
	end

	%	Maximise over remaining variables
	[~, ~, ls, ~, ds, fs] = size(z); % number of locations, directions, freq.
	v = {'dir', 'freqS'}; % names of response variables across which to maximise
	i = contains(dim, v); % resp. dimensions of maximising variables: truth
	i = find(i); % response dimensions of maximising variables: indices
	[zMax, i] = max(z, [], i, 'linear'); % maxima and their resp. dim. indices

	%	Calculate tuning curves of maximising variables
	zSub = reshape(z, [], ds, fs); % reshape for ind2sub: ls x ds x fs
	[~, iDir, iFreq] = ind2sub(size(zSub), i);
		% preferred values of variables (indices): 1 x 1 x ls
	zDir = zeros(ls, ds); zFreq = zeros(ls, fs); % allocate storage
	for l = 1: ls % loop over locations
		zDir(l, :) = z(:, :, l, :, :, iFreq(l)); % direction tuning: ls x ds
		zFreq(l, :) = z(:, :, l, :, iDir(l), :); % spatial frequency resp.: ls x fs
	end

	%	Store
	if isfield(m.doMax, 'out') % choose output variable
		z = m.doMax.out; % user-specified name
	else, z = m.z; % default: same as input name
	end
	switch z
		case m.z % default
			d.(m.z) = zMax; % store maximised response: 1 x 1 x ls
			dim = dim(1: 3); % update output dimensions
		case 'dirPref' % preferred direction
			iDir = shiftdim(iDir, 1); % preferred direction (index): 1 x ls
			d.dirPref = d.dir(iDir); % preferred direction (deg): 1 x ls
			dim = {'', 'loc'}; % update output dimensions
		case 'freqSPref' % preferred spatial frequency
			iFreq = shiftdim(iFreq, 1); % preferred spatial frequency (index): 1 x ls
			d.freqSPref = d.freqS(iFreq); % preferred spatial freq. (cycles/deg): 1 x ls
			dim = {'', 'loc'}; % update output dimensions
			d.Properties.VariableDescriptions{'freqSPref'} = ...
				'Preferred spatial frequency (cycles/deg)';
		case 'dirTun' % direction tuning
			d.dirTun = shiftdim(zDir, -1); % direction tuning: 1 x ls x ds
			dim = {'', 'loc', 'dir'}; % update output dimensions
		case 'freqSTun' % spational frequency tuning
			d.freqSTun = shiftdim(zFreq, -1); % store spatial frequency resp.: 1 x ls x fs
			dim = {'', 'loc', 'freqS'}; % update output dimensions			
	end
	d.Properties.CustomProperties.RespDim = dim; % update output dimensions

function [d, m] = field(d, m) % calculate receptive field by cross correlation

% Method: deliver a range of stimuli, weight each by response, and add

	% Obtain the stimulus characteristics and responses
	[name, val, vals] = m.getVal(m, 'resp'); % model parameter names and values
	ss = prod(vals); % total number of stimuli
	r = d.resp; % 4CBeta excitatory cell impulse rate (Hz):
		%	1 x ts x ls x 1 x ss, ls = number of neurons
	r = max(r, [], 2); % maximum response (Hz): 1 x 1 x ls x 1 x ss
	r = reshape(r, [], ss); % prepare for cross correlation: ls x ss

	% Calculate the stimulus element locations
	xs = 101; % number of x values
	x = .5 * m.p.wid * linspace(-1, 1, xs); % x values: 1 x xs
	[xG, yG] = ndgrid(x, x); % grid locations (deg): xs x ys
	loc = [xG(:), yG(:)]; % locations (x, y) (deg): %	ps x 2
	ps = size(loc, 1); % number of stimulus pixels
	
	%	Cross correlate
	s = zeros(ss, ps); % stimuli: ss x ps
	for i = 1: ss % loop over stimuli
		m = m.setVal(m, name, val, i); % set values of model parameters
		stim = m.calStim(0, loc, m); % stimulus (contrast-units): 1 x ps x cs
		s(i, :) = stim(1, :, 1); % assume that the stimulus is achromatic: ss x ps
	end
	f = r * s; % receptive field (Hz x c.-u.): (ls x ss) * (ss x ps) = ls x ps
	f = f / (ss * m.p.contMag); % normalise (Hz): ls x ps

	%	Store
	d.locS = shiftdim(loc, -1); % pixel locations (deg): 1 x ps x 2
	d.resp = shiftdim(f, -1); % r.f. (Hz): 1 x ls x ps
	d.Properties.CustomProperties.RespDim = {'', 'loc', 'locS'}; % resp. dim.

function [d, m] = fig(d, m) % calculate functions for figure plot

	n = fieldnames(m.fig)'; % cell array of names including function types: 1 x ts
	for fun = string(n) % loop over function types
		switch fun % choose function type
			case 'conv' % convergence functions
	
				% Calculate convergence functions
				xs = 101;
				if isfield(m.fig, 'xLim') % x limit (deg)
					xLim = m.fig.xLim; % specified by user
				else, xLim = .1; % default
				end
				x = xLim * linspace(-1, 1, xs); % visual field locations (deg): 1 x xs
				rad = [m.p.radOpt; m.p.radHor; m.p.radGang];
					% radii (deg): rs x 1
				f = exp(- x .^ 2 ./ rad .^ 2); % convergence function: rs x xs
			
				% Create a table for each radius
				source = ["opt", "hor", "gang"]; ss = length(source);
				dC = cell(1, ss); % tables
				for i = 1: ss % loop over sources
					dC{i} = table(source(i), x, rad(i), f(i, :), ...
						'variableNames', {'source', 'x', 'rad' 'f'}); % make table
				end
			
				% Concatenate tables
				d = vertcat(dC{:}); % concatenate
				d.source = categorical(d.source); % for listing
				d.Properties.VariableDescriptions{'source'} = 'Function source';
				d.Properties.VariableDescriptions{'x'} = 'Visual field location (deg)';
				d.Properties.VariableDescriptions{'rad'} = 'Function radius (deg)';
				d.Properties.VariableDescriptions{'f'} = 'Convergence function';
	
			case 'course' % stimulus time course
	
				s = d.stim; % stimulus: 1 x ts x xs x ys
				t = d.t; % times (s): 1 x ts
				i = t <= 1 / m.p.freqT; % reduce to first cycle (s): 1 x ts, new ts
				x = d.x; y = d.y; % locations (deg): 1 x xs
				xs = length(x); % number of locations
				[x, y] = ndgrid(x, y); % location grids (deg): xs x ys
				loc = [x(:), y(:)]; % locations (deg): xs * ys x 2
				j = knnsearch(loc, [0, 0]); % find location nearest to middle of visual field
				[j, k] = ind2sub([xs, xs], j); % middle location (deg): 1 x 2
				s = s(1, i, j, k); % time course of stimulus at middle: 1 x ts
				d.t = t(i); % first cycle only (s): 1 x ts
				d.course = s; % store: 1 x ts
	
			case 'mech' %  % circles with radii equal to centre, surround

				r = [m.p.radCen, m.p.radSur]; % centre, surround radii (deg): 1 x 2
				n = 101; % number of points around circle
				theta = linspace(0, 2 * pi, n)'; % angles around circle (rad): n x 1
				x = exp(1i * theta); % vectors around circle: n x 1
				y = r .* x; % complex values: n x 2
				%	y = [y(:, 1), y(2, :)]; % row vector: 1 x 2 * n
				loc = m.fig.mech; % circle location
				if ~ isempty(loc), y = complex(loc(1), loc(2)) + y; end
				d = table(theta, y); % create table
				d.ecc = repmat(m.p.ecc, [n, 1]); % add eccentricity
				d.array = repmat("cone", [n, 1]); % add array
	
		end
	end

function [d, m] = fitGauss(d, m) % fit Gaussian to spatial frequency response

	b(1) = d.resp(1); % peak response
	rCen = sqrt(m.p.radOpt ^ 2 + m.p.radGang ^ 2); % estimated centre rad. (deg)
	rSur = sqrt(2 * m.p.radHor ^ 2); % estimated surround radius (deg)
	switch m.fitGauss.comp % centre or surround?
		case 'cen', r = rCen; % centre component
		case 'sur', r = rSur; % surround component
	end
	switch m.fitGauss.coef % number of regression coefficients
		case 1 % peak response only
			fun = @(b, f) b * exp(-.25 * (r * 2 * pi * f) .^ 2); % Gaussian
		case 2 % peak response and radius
			b(2) = r; % radius (deg)
			fun = @(b, f) b(1) * exp(-.25 * (b(2) * 2 * pi * f) .^ 2); % Gaussian
		case 4 % both centre and surround components
			b(2) = rCen; b(3) = .5 * b(1); b(4) = rSur; % estimate coefficients
			fun = @(b, f) b(1) * exp(-.25 * (b(2) * 2 * pi * f) .^ 2) - ...
				b(3) * exp(-.25 * (b(4) * 2 * pi * f) .^ 2); % difference of Gaussians
	end
	model = fitnlm(d, fun, b, 'predictorVars', {'freqS'}, 'responseVar', 'resp');
		% fit data
	m.model{m.group} = model;
	
function [d, m] = fund(d, m) % calculate fundamental Fourier response

% Assumptions:
%		* d has 1 row
%		* frequency is the second dimension of d.resp

	% Calculate fundamental
	r = d.resp; % response: 1 x fs x vs where vs can be multi-dimensional
	s = size(r); % response size
	fs = s(2); % number of frequencies
	r = r(1, 2, :); % keep only fundamental component: 1 x 1 x vs
	r = reshape(r, [1, 1, s(3: end)]); % restore shape
	r = (2 / fs) * r; % convert from fft units to temporal units: 1 x 1 x vs
	
	% Store
	if isfield(m.fund, 'prop') % response property: amplitude, complex or phase
		p = m.fund.prop; % user-specified
	else
		p = 'complex'; % default: leave response as complex
	end
	switch p % amplitude, complex or phase?
		case 'amp' % amplitude
			d.resp = abs(r); % store amplitude: 1 x fs x vs
			%	desc = 'Fundamental amplitude'; % response description
		case 'complex'
			d.resp = r; % store: 1 x fs x vs
			%	desc = 'Fundamental component';
		case 'phase'
			d.resp = (180 / pi) * angle(r); % store phase (deg): 1 x fs x vs
			d.Properties.VariableDescriptions{'resp'} = 'Fundamental phase (deg)';
	end
	d.freq = d.freq(2); % update frequency

function [d, m] = hist(d, m) % calculate histogram

	% Initialise
	if isfield(m.hist, 'edges') % set edges or number of bins
		arg = m.hist.edges; % user-specified
	else
		arg = 10; % number of bins
	end
	
	% Calculate histogram
	x = d.(m.x); % data values
	[n, e] = histcounts(x, arg); % obtain histogram counts and edges
	w = e(2) - e(1); % bin width
	c = .5 * w + e(1: end - 1); % bin centres
	d = d(1, :); % keep only the first row
	d.(m.x) = c; % bin centres
	d.edge = e; % bin edges
	d.(m.z) = n; % counts
	d.Properties.VariableDescriptions{'edge'} = 'Edge';
	d.Properties.VariableDescriptions{m.z} = 'Count';

function [d, m] = interp(d, m) % fit tuning curve with an interpolated model

	% Calculate x values for interpolation
	x = d.(m.x); % tuning variable (deg or cycles/deg): 1 x vs
	z = d.(m.z); % response: 1 x ls x vs
	z = shiftdim(z, 1)'; % transpose for interp1: vs x ls
	if m.x == "dir" % direction is cyclic: add starting point at end, for interp.
		x(end + 1) = x(1) + 360; % close open-ended interval (deg): 1 x (vs + 1)
		z(end + 1, :) = z(1, :); % cyclic z: (vs + 1) x ls
	end
	if isfield(m.interp, 'xs') % set number of interpolation points
		xs = m.interp.xs; % user-defined
		xs = 2 * floor(.5 * xs) + 1; % make it odd, for finding bandwidth later
	else, xs = 101; % default
	end
	xI = linspace(x(1), x(end), xs); % interpolation points: 1 x xs

	% Interpolate, find maximum and preferred direction
	if isfield(m.interp, 'method'), method = m.interp.method; % set method
	else, method = 'spline'; end % default
	z = interp1(x, z, xI', method); % interpolate: xs x ls
	[~, i] = max(z); % maximum response, and its index: 1 x ls
	xPref = xI(i); % preferred value of tuning curve: 1 x ls

	%	Store
	d.(m.x) = xI; % tuning variable: 1 x xs
	zC = permute(z, [3, 2, 1]); % response: 1 x ls x xs
	if isfield(m.interp, 'out') % choose output variable
		z = m.interp.out; % user-specified name
	else, z = m.z; % default: same as input name
	end
	switch z
		case 'dirPref' % preferred direction
			d.dirPref = xPref; % preferred direction (deg): 1 x ls
			dim = {'', 'loc'}; % update output dimensions
		case 'freqSPref' % preferred spatial frequency
			d.freqSPref = xPref; % preferred spatial freq. (cycles/deg): 1 x ls
			dim = {'', 'loc'}; % update output dimensions
		case 'dirTun' % direction tuning
			d.dirTun = zC; % direction tuning: 1 x ls x ds
			dim = {'', 'loc', 'dir'}; % update output dimensions
		case 'freqSTun' % spational frequency tuning
			d.freqSTun = zC; % store spatial frequency resp.: 1 x ls x fs
			dim = {'', 'loc', 'freqS'}; % update output dimensions			
	end
	d.Properties.CustomProperties.RespDim = dim; % update output dimensions

function [d, m] = interpWid(d, m) % calculate tuning bandwidth *** fix ***

	% Read data, close direction interval, calculate x values for interpolation
	dir = d.(m.x); % direction (deg): 1 x ds
	dir(end + 1) = dir(1) + 360; % close open-ended interval: 1 x (ds + 1)
	z = d.(m.z); % response: 1 x 1 x ... x ds
	z = reshape(z, 1 , []); % remove extraneous dimensions: 1 x ds
	z(end + 1) = z(1); % assume cyclic z: 1 x (ds + 1)
	if isfield(m.interp, 'xs') % set number of interpolation points
		xs = m.interp.xs; % user-defined
		xs = 2 * floor(.5 * xs) + 1; % make it odd
	else, xs = 101; % default
	end
	x = linspace(dir(1), dir(end), xs + 1); % interp. points: 1 x (xs + 1)
	x = x(1: xs); % open-ended interval: 1 x xs

	% Interpolate, find maximum and preferred direction
	z = interp1(dir, z, x, m.interp.method); % interpolate: 1 x xs
	[zMax, i] = max(z); % maximum response, and its index
	dirPref = x(i); % preferred direction

	%	Calculate bandwidth
	iC = ceil(.5 * xs); % central index
	zCent = circshift(z, iC - i); % centre max. *** remove repeated point at end
	zCrit = zMax / sqrt(2); % bandwidth criterion
	w = zeros(2, 1); % allocate storage
	for j = 1: 2 % left then right half-bandwidths
		zC = zCent; % response with centred maximum
		if j == 1, zC = flip(zC); end % flip left for right
		zC = zC(iC: end); % retain reponse to right of maximum
		i = zC > zCrit; % indices of values above bandwidth criterion
		i = find(~ i, 1); % index of first value below criterion
		i = i - 1: i; % indices containing criterion
		i = interp1(zC(i), i, zCrit); % index of criterion
		w(j) = (i - 1) * (x(2) - x(1)); % bandwidth for directions > maximum (deg)
	end
	w = mean(w); % average left and right half-bandwidths

	%	Store
	d.(m.x) = x; % interpolated directions: 1 x xs
	d.(m.z) = z; % model values: 1 x xs
	d.respMax = zMax; % maxiumum response: 1 x 1
	d.dirPref = dirPref; % direction at which response is maximum: 1 x 1
	d.band = w; % bandwidth (deg)

function [d, m] = map(~, m) % prepare for map of neuron locations

	for i = string(m.p.array) % loop over arrays
		aC = m.p.(i); % current array
		l = aC.loc; % cell locations (deg): ls x 2
		ls = size(l, 1); % number of cells
		switch i % set colours
			case 'cone' % cones
				c = eye(3); % RGB colour map for L, M, and S cones: 3 x 3
				j = aC.type; % cell type: ls x 1
				c = c(j, :); % marker colours: ls x 3
			case 'gangOff' % off-centre ganglion cells
				c = [0, 0, 1]; % blue for off-centre: 1 x 3
				c = repmat(c, [ls, 1]); % one row for each cell: ls x 3
			case 'gangOn' % on-centre ganglion cells
				c = [1, 0, 0]; % red for on-centre: 1 x 3
				c = repmat(c, [ls, 1]); % one row for each cell: ls x 3
		end
		ecc = repmat(m.p.ecc, [ls, 1]); % eccentricity (deg)
		array = repmat(i, [ls, 1]); % array name
		array = categorical(array); % make it categorical
		x = l(:, 1); % x locations (deg): ls x 1
		y = l(:, 2); % y locations (deg): ls x 1
		dC.(i) = table(ecc, array, x, y, c); % store
	end
	dC = struct2cell(dC); % convert to cell for listing
	d = vertcat(dC{:}); % combine arrays
	d.Properties.VariableDescriptions{1} = ''; % set all descriptions to empty

function [d, m] = prep(d, m) % prepare for plot

	%	Prepare
	n = fieldnames(m.prep)'; % cell array of names including preparations: 1 x ps
	for p = string(n) % loop over preparations
		switch p % choose preparation
			case 'image' % prepare for plotting with inbuilt functions image/sc

				%	Inputs:
				%		m.x = name of location variable
				%		m.prep.image = name of image variable: ps x ls x qs, where l is
				%		a list of locations, and ps and qs can be multidimensional

				%	Initialise
				loc = m.x; % name of location variable
				[lim(1), lim(2)] = bounds(d.(loc), 'all'); % lower and upper loc. limits
				d.x = lim; d.y = lim; % store: 1 x 2
				%	d.x = .5 * d.width * [-1, 1]; d.y = d.x; % visual field limits: 1 x 2
				dim = d.Properties.CustomProperties.RespDim; % response dimensions
				i = contains(dim, loc); % true for location dimension
				j = find(i); % index of location dimension
				iPre = 1: j - 1; iPost = j + 1: length(dim); % dim. before & after loc.
		
				% Reshape response
				im = m.prep.image; % name of image variable
				z = d.(im); % image variable
				s = size(z); % image size
				zs = s(i); % number of locations
				xs = round(sqrt(zs)); % number of loc. on a side, assuming square array
				z = reshape(z, [s(iPre), xs, xs, s(iPost)]); % map: ps x xs x ys x qs
				z = permute(z, [iPre, j + 1, j, iPost + 1]); % image: ps x ys x xs x qs
				d.(im) = z; % store
				dim = [dim(iPre), 'x', 'y', dim(iPost)]; % update dimensions
				d.Properties.CustomProperties.RespDim = dim; % store

			case 'norm' % normalise by maximum
				v = d.(m.z); % values of variable: 1 x vs
				d.total = sum(v); % total number of values: 1 x 1
				maxC = max(v, [], 2); % 1 x 1
				d.(m.z) = v ./ maxC; % normalised values: 1 x vs
				d.Properties.VariableDescriptions{'total'} = 'Total count';
				d.Properties.VariableDescriptions{m.z} = 'Normalised count';
			case 'orient' % convert direction to orientation
				z = d.(m.z); % direction (deg), range = [-180, 180)
				z = z + 90; % range = [-90, 270)
				z = mod(z, 180); % range = [0, 180)
				d.(m.z) = z - 90; % orientation (deg), range = [-90, 90)
			case 'real' % make resp real *** check ***
				d.resp = real(d.resp); % make it real
			case 'resize' % resize image to, for example, smooth it *** fix ***

				% Obtain inputs
				if isempty(m.prep.resize) % scale is specified by user
					s = 4; % default: image width and height will be multiplied by s
				else, s = m.prep.resize; % user-specified
				end
			
				% Scale
				z = d.(m.z); % image: 1 x xs x ys x ps, ps can be multi-dimensional
				z = shiftdim(z, 1); % shift image to first two dimensions: xs x ys x ps
				xs = size(z, 1); % number of x values
				xs = s * (xs - 1) + 1; % increased number of x values
				z = imresize(z, [xs, xs]); % interpolate: xs x ys x ps
				%	x = d.xStim; x = linspace(x(1), x(end), xs); % recalculate x values
				d.(m.z) = shiftdim(z, -1); % store
				%	d.xStim = x; % store

			case 'unwrap' % unwrap phase
				r = d.resp; % response: 1 x vs (vs can be multidimensional)
				dimU = m.prep.unwrap; % response dimension to unwrap
				dim = d.Properties.CustomProperties.RespDim; % response dimensions
				[~, i] = ismember(dimU, dim); % index of dimension to unwrap
				r = pi * r / 180; % convert to radians
				r = unwrap(r, [], i); % unwrap
				d.resp = 180 * r / pi; % convert back to degrees
		end
	end

function [d, m] = radius(d, m) % calculate mechanism radius

	% Calculate and store radius
	f = d.freqS; % spatial frequency (cycles/deg): 1 x fs
	r = shiftdim(d.resp, 3); % spatial frequency response (mV): 1 x fs
	fC = interp1(r, f, r(1) * exp(-1)); % interpolate for 1/e point (cycles/deg)
	rad = 1 / (pi * fC); % radius (deg)
	d.ecc = m.p.ecc; d.radius = rad; % store
	d.Properties.VariableDescriptions{'ecc'} = 'Eccentricity (deg)';
	d.Properties.VariableDescriptions{'radius'} = 'Radius (deg)';

	% Append file
	dC = d; % current data
	file = m.save.name; % name of file to which current data will be appended
	if exist(file, 'file') % there are existing data
		load(file, 'd'); % load the existing data
		d = [d; dC]; % append current data
	end

function [d, m] = resp(~, m) % calculate the model time course

% Inputs:
%		m.resp.par, where par is a model parameter: matrix of parameter values,
%			with one row per value, e.g. m.resp.cont = [1, 1, 0; 0, 0, 1]
%		m.resp.stage: names of stages to return: char array or cell array of chars
%	Outputs:
%		Single neuronal array: a single-row table with all responses in a
%			single table variable
%		Multiple neuronal arrays: a multi-row table, with one row for each
%			location in each array

	% Find the model parameters to be varied
	[name, val, vals] = m.getVal(m, 'resp'); % model parameter names and values
	vs = prod(vals); % total number of values

	%	Obtain the temporal values
	switch m.p.domain % set temporal variable
		case 'freq', t = m.p.f; % frequency (Hz): 1 x ts
		case 'time', t = m.p.t; % time (s): 1 x ts
	end
	ts = length(t); % number of frequencies or times
	
	% Initialise the response array
	if ~ iscell(m.resp.stage) % m.resp.stage isn't a cell array
		m.resp.stage = {m.resp.stage}; % make it a cell array
	end
	i = contains({m.p.cell.type}, m.resp.stage); % indices of required stages
	stage = {m.p.cell(i).type}; % stage names
	ss = length(stage); % number of stages
	array = {m.p.cell(i).array}; % names of corresponding neuronal arrays
	array = unique(array); % names of unique neuronal arrays
	arrays = length(array); % number of neuronal arrays
	if arrays == 1 % initialise for one neuronal array
		array = m.p.(array{1}); % required array, struct
		loc = shiftdim(array.loc, -1); % channel locations (deg): 1 x ls x 2
		ls = size(loc, 2); % number of locations in array
		p = zeros(ts, ls, ss, vs); % allocate storage
		dim = 1; % temporal dimension
	else % initialise for multiple neuronal arrays
		s = cell(1, ss); % allocate storage for stage names
		l = cell(1, ss); % allocate storage for locations in each stage
		for i = 1: ss % loop over stages
			aC = string(array{i}); % name of this stage's array
			l{i} = m.p.(aC).loc; % store locations for this stage (deg)
			sC = stage{i}; % name of current stage
			s{i} = repmat(string(sC), [size(l{i}, 1), 1]); % store stage name for each loc.
		end
		s = vertcat(s{:}); % concatenate array names: ls x 1
		s = categorical(s); % save space
		loc = vertcat(l{:}); % concatenate locations (deg): ls x 2
		ls = size(loc, 1); % number of locations across all required stages
		p = zeros(ls, ts, vs); % allocate storage
		dim = 2; % temporal dimension
	end

	% Calculate and store the responses
	for v = 1: vs % loop over variable values

		% Set the model parameter values, and the temporal parameters
		m = m.setVal(m, name, val, v); % set values of model parameters
		m = m.calTemp(m); % recalculate temporal parameters

		% Calculate the response
		switch m.p.solver % choose solver
			case 'solveF' % solve in frequency domain
				pC = m.solveF(m); % potential (mV), struct: 1 x ss, all stages in model
				dom = "freq"; % frequency domain response
			case 'solveT' % solve in time domain
				pC = m.solveT(m); % potential (mV), struct: 1 x ss, all stages in model
				dom = "time"; % temporal domain response
		end

		%	Store the response
		j = contains({pC.type}, stage); % indices of stages to keep
		pC = pC(j); % keep required stages, struct: 1 x ss, reduced ss
		if arrays == 1 % store one neuronal array
			pC = cat(3, pC.resp); % concatenate: ls x ts x ss
			pC = permute(pC, [2, 1, 3]); % data file format: ts x ls x ss
			p(:, :, :, v) = pC; % potential (mV): ts x ls x ss x vs
		else % store multiple neuronal arrays
			pC = cat(1, pC.resp); % concatenate: ls x ts, total ls across stages
			p(:, :, v) = pC; % potential (mV): ls x ts x vs
		end

	end

	% Convert to impulse rate if required
	if ~ m.p.pot % convert to impulse rate
		if dom == "freq" % potential is in the frequency domain
			p = m.ifftReal(p, dim); % convert to temporal domain
			dom = "time"; % temporal domain response
		end
		p = m.p.kRect * max(0, p); % impulse rate (Hz)
	end

	% Set output response domain
	switch m.p.domain % output frequency or temporal response?
		case 'freq' % frequency response
			if dom == "time", p = fft(p, [], dim); end % convert temporal response
		case 'time' % temporal response
			if dom == "freq", p = ifft(p, [], dim); end % convert frequency response
	end
	
	% Store variables and response in table
	% *** check multi-row output ***
	ecc = m.p.ecc; width = m.p.wid; % eccentricity (deg), visual field width (deg)
	d = table(ecc, width, t); % store
	if arrays == 1 % store one neuronal array
		d.loc = loc; d.stage = stage; % add location, stage
		for i = 1: length(name) % loop over variables
			nameC = name{i}; % current variable
			d.(nameC) = shiftdim(val{i}, -1); % store its values in table
		end
		d.resp = reshape(p, [1, ts, ls, ss, vals]); % store potential:
			% 1 x ts x ls x ss x vs, where vs can be multi-dimensional
		dim = [{'', m.p.domain, 'loc', 'stage'}, name']; % response dimensions
	else % store multiple neuronal arrays
		d = repmat(d, [ls, 1]); % one row for each location
		d.stage = s; d.loc = loc; % add stage, location
		for i = 1: length(name) % loop over variables
			nameC = name{i}; % current variable
			valC = shiftdim(val{i}, -1); % values of this variable
			d.(nameC) = repmat(valC, [ls, 1]); % store values in table
		end
		d.resp = p; % store response: ls x ts x vs, vs can be multi-dimensional
		dim = ['loc', m.p.domain, name']; % response dimensions
	end
	
	% Set table's properties
	d.Properties.VariableNames{'t'} = m.p.domain; % make the name more meaningful
	d.Properties.Description = m.p.project; % research project
	if m.p.pot % response quantity is generator potential
		d.Properties.VariableDescriptions{'resp'} = 'Generator potential (mV)';
	else
		d.Properties.VariableDescriptions{'resp'} = 'Impulse rate (Hz)';
	end
	d = addprop(d, 'RespDim', 'table'); % add property: response dimensions
	d.Properties.CustomProperties.RespDim = dim; % response dimensions
	
function [d, m] = showMovie(d, m)
		
	% Show stimulus movie
	
	% Initialise
	s = d.stim; % stimulus: 1 x ts x xs x ys x cs
	s = shiftdim(s, 1); % remove first dimension: ts x xs x ys x cs
	frames = size(s, 1); % number of movie frames
	x = d.x; y = x; % x and y values
	a = gca; % axes handle
	a.NextPlot = 'replaceChildren'; % retain axis settings during movie

	% Compile and show movie
	switch 'loop'
		case 'loop'
			figure(gcf); % make the figure visible
			for i = 1: frames % loop over frames
				sC = shiftdim(s(i, :, :, :), 1); % current stimulus: xs x ys x cs
				image(x, y, sC); % draw image
				pause(.05); % pause between frames (s)
			end
		case 'movie'
			f = figure(gcf); % figure handle
			f.Visible = 'off'; % hide the figure
			mov(frames) = struct('cdata',[],'colormap',[]); % allocate storage
			for i = 1: frames % loop over frames
				sC = shiftdim(s(i, :, :), 1)'; % current stimulus, x is horizontal
				imagesc(x, y, sC); % draw image
				drawnow; % draw the frame
				mov(i) = getframe; % store the frame
			end
			f.Visible = 'on'; % show the figure
			movie(mov); % play the movie
	end
	
function [d, m] = stim(~, m) % calculate the stimulus

	% Read or calculate the stimulus and animate it
	t = m.p.t; % simulation times: 1 x ts
	ts = length(t); % number of times
	if isfield(m.stim, 'xs') % set the number of x values
		xs = m.stim.xs; % user-specified
	else
		xs = 101; % default
	end
	x = m.p.widH * linspace(-1, 1, xs); % x values: 1 x xs
	[xG, yG] = ndgrid(x, x); % grid locations (deg): xs x ys
	loc = [xG(:), yG(:)]; % locations (x, y) (deg): ls x 2
	s = m.calStim(t, loc, m); % stimulus: ts x ls x cs
	s = reshape(s, 1, ts, xs, xs, []); % image: 1 x ts x xs x ys x cs
	s = permute(s, [1, 2, 4, 3, 5]); % permute from matrix view to user view
	s = .5 * (1 + s); % change range [-1, 1] to [0, 1]: 1 x ts x xs x ys x cs
	d = table(t, x, x, s, 'variableNames', {'t', 'x', 'y', 'stim'}); % store
