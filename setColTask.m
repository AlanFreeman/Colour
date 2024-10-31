function m = setColTask(m) % define function handles for setCol tasks

	h = localfunctions; % handles of functions in this file
	hs = length(h); % number of handles
	for i = 1: hs % loop over handles
		hC = h{i}; % current handle
		task = func2str(hC); % name of task
		m.(task).fun = hC; % store handle in metadata
	end

function [d, m] = count(d, m) % cones per midget ganglion cell dendritic field

	% Get cone density, Wässle (89), and g.c. dend. field radius, Watanabe (89)
	eF = d.eccDeg; % functional eccentricity (deg): es x 1
	dens = d.densDeg; % cone density (deg^-2): es x 1
	rat = m.p.ratCone; % ratio of cone types, [L, M, S]
	dens = (1 - .5 * rat(3)) * dens; % remove on-S-bipolar connections
	r = polyval(m.p.radGangCoef, eF); % radius (deg): es x 1
	
	% Count
	switch 'area' % select calculation method
		case 'area' % multiply cone density by centre area
			area = pi * r .^ 2; % centre area (deg^2): es x 1
			count = dens .* area; % cones per ganglion cell field: es x 1
		case 'gaussian' % weight cones by Gaussian profile, and add
			s = 1 ./ sqrt(sin(pi / 3) * d.dens); % cone separation (deg): es x 1
			w = 6 * r; % array width (deg): Gaussian is calc. to 3 times rad.: es x 1
			es = length(s); rat = zeros(es, 1); % allocate storage
			for i = 1: es % loop over eccentricities
				l = array(s(i), w(i)); % locations of cones (deg): ls x 2
				l = sum(l .^ 2, 2); % squared dist. of cone from middle (deg^2): ls x 1
				l = exp(- l ./ r(i) ^ 2); % Gaussian weighting: ls x 1
				rat(i) = sum(l(:)); % cones per centre mechanism
			end
	end
	d.area = area; % dendritic field area (deg^2): es x 1
	d.count = count; % cones per ganglion cell dendritic field (deg): es x 1
	d.Properties.VariableDescriptions{'area'} = ...
		'Ganglion cell dendritic field area (deg^2)';
	d.Properties.VariableDescriptions{'count'} = ...
		'Cones per ganglion cell dendritic field';

function [d, m] = countCen(d, m) % cones per ganglion cell centre

	% Obtain radius data
	[dC, m] = rad(d, m); % calculate radii
	i = dC.source == 'cenConv'; % indices for centre radii
	dC = dC(i, :); % select centre data: es x 1
	r = dC.radius; % centre radius (deg): es x 1
	es = size(dC, 1); % number of eccentricities

	% Obtain density data from input table
	ecc = d.eccDeg; % eccentricity (deg): eDs x 1, eDs = number of densities
	dens = d.densDeg; % cone density (cones/deg^2): eDs x 1
	dens = interp1(ecc, dens, dC.ecc); % interp. density data (deg^-2): es x 1
	c = pi * r .^ 2 .* dens; % area * density = cone count: es x 1
	d = repmat(d(1, :), [es, 1]); % change number of rows to es
	d.eccDeg = dC.ecc; d.count = c; % store
	d.Properties.VariableDescriptions{'count'} = 'Cones per ganglion cell centre';

function [d, m] = countCone(d, m) % cone count vs eccenticity: Packer (89)

	i = d.quad == 'temporal'; d = d(i, :); % keep only temporal quadrant
	e = d.eccDeg; % eccentricity (deg): es x 1
	dens = d.densDeg; % density (deg^-2): es x 1
	dn = 2 * pi * e .* dens; % derivative of cell count (deg^-1): es x 1
	n = cumtrapz(e, dn); % integrate to find cell count: es x 1
	d.count = n; % store cell count
	d.Properties.VariableDescriptions{'count'} = 'Cone count';

function [d, m] = countGang(d, m) % g.c. count vs fun. ecc.: Wässle (89)

	% Initialise
	i = d.type == 'func'; % functional eccentricity
	d = d(~ i, :); % keep only anatomical eccentricity
	d = sortrows(d, 'eccDeg'); % sort on eccentricity
	eA = d.eccDeg; % anatomical eccentricity (deg): es x 1
	if isfield(m.p, 'densGangCoef') % regression coefficients are available
		densA = exp(polyval(m.p.densGangCoef, log10(eA))); % dens. (deg^-2): es x 1
	else
		densA = d.densDeg; % density (deg^-2): es x 1
	end
	
	% Integrate density to find cell count
	switch 'cont' % method
		case 'cont' % continuous
			dn = 2 * pi * eA .* densA; % derivative of cell count (deg^-1): es x 1
			n = cumtrapz(eA, dn); % integrate to find cell count: es x 1
		case 'disc' % discrete
			dE = diff(eA); % eccentricity increments (deg); (es - 1) x 1
			dE = [dE(1); dE]; % replicate first point (deg): es x 1
			a = pi * (eA - .5 * dE) .^ 2; % circle areas (deg^2): es x 1
			a = [a; pi * (eA(end) + .5 * dE(end)) ^ 2]; % add final pt.: (es + 1) x 1
			dA = diff(a); % annulus area (deg^2): es x 1
			n = dA .* densA; % annulus count: es x 1
			n = cumsum(n); % cumulative count: es x 1
	end

	% Convert anatomical to functional eccentricity, and store
	s = load([m.folder, filesep, 'Functional eccentricity'], 'd'); dE = s.d;
	e = interp1([0; dE.ecc], [0; dE.eccFun], eA);
		% functional eccentricity (deg), extended to origin: es x 1
	d.ecc = e; % store functional eccentricity
	d.count = n; % store ganglion cell count
	i = d.eccDeg == 0 | d.eccDeg >= dE.ecc(1); % rows in conversion range
	d = d(i, :); % keep them

	%	Describe
	d = renamevars(d, 'eccDeg', 'eccAnat'); % rename for clarity
	d.Properties.VariableDescriptions{'ecc'} = 'Functional eccentricity (deg)';
	d.Properties.VariableDescriptions{'eccAnat'} = ...
		'Anatomical eccentricity (deg)';
	d.Properties.VariableDescriptions{'count'} = 'Ganglion cell count';

function [d, m] = cover(d, m) % ganglion cell coverage: multiple authors

% Input: ganglion cell density versus functional eccentricity

	dens = d.densDeg; % density (cells/deg^2)
	switch 'cen' % source of radius
		case 'cen' % ganglion cell centre radius
			dC = rad(d, m); % calculate centre radius
			i = dC.source == 'cenConv'; dC = dC(i, :); % select centre radius
			r = dC.radius; % centre radius (deg)
			dens = interp1(d.ecc, dens, dC.ecc); % interpolate dens. at radius ecc.
			d = dC; % use radius table
		case 'dend' % ganglion cell dendrite radius
			e = d.ecc; % functional eccentricity (deg)
			r = polyval(m.p.radGangCoef, e); % g.c. dend. radius (deg)
	end
	d.cover = pi * r .^ 2 .* dens; % coverage = area * density
	d.Properties.VariableDescriptions{'cover'} = 'Ganglion cell coverage';

function [d, m] = densCone(~, m) % cone density: Packer (89)

	% Load and adjust data table
	folder = [m.folder, filesep, 'Packer (89) dens']; % data folder
	d = m.readFile(folder); % create table from folder
	s = split(d.label, ", "); % split the label into quad, near
	d.quad = s(:, 1); d.near = s(:, 2); % store
	d.quad = categorical(d.quad); d.near = double(d.near); % fix the classes
	d = removevars(d, 'label'); % no longer required
	d = movevars(d, {'quad', 'near'},'before', 'ecc'); % match old order
	d = sortrows(d, {'quad', 'ecc'}); % sort rows
	
	% Convert and store
	s = 1 - .022; % linear shrinkage
	d.ecc = (1 / s) * d.ecc; % eccentricity, corrected for shrinkage (mm)
	d.eccDeg = m.p.magRet * d.ecc; % eccentricity (deg)
	s = 1 - .05; % areal shrinkage
	dens = s * d.dens; % cone density, corrected for shrinkage (.001 x mm^-2)
	d.dens = 1000 * dens; % cone density (mm^-2)
	d.densDeg = (m.p.magRet ^ -2) * d.dens; % cone density (deg^-2)
	d.Properties.VariableDescriptions = {'Retinal quadrant', 'Near fovea', ...
		'Eccentricity (mm)', 'Cone density (mm^-^2)', 'Eccentricity (deg)', ...
		'Cone density (deg^-^2)'};

function [d, m] = densGang(~, m) % ganglion cell density: Wässle (89)

	% Load the empirical data, and convert mm to degrees
	folder = [m.folder, filesep, 'Wässle (89)']; % data folder
	d = m.readFile(folder); % create table from folder
	d = renamevars(d, 'label', 'type'); % rename label
	d.type = categorical(d.type); % for better listing
	d = sortrows(d, {'type', 'ecc'}); % sort by eccentricity
	d.dens = 10 .^ d.densLog; % density (mm^-2)
	eccDeg = m.p.magRet * d.ecc; % eccentricity (deg)
	densDeg = (m.p.magRet) ^ -2 * d.dens; % density (deg^-2)

	% Correct for shrinkage, and store
	i = d.type == 'anatNarr'; % anatWide and func data are already corrected
	s = .9; % linear shrinkage
	eccDeg(i) = (1 / s) * eccDeg(i); % shift away from fovea
	densDeg(i) = s ^ 2 * densDeg(i); % reduce density
	d.eccDeg = eccDeg; d.densDeg = densDeg; % store
	d.Properties.VariableDescriptions{'eccDeg'} = 'Eccentricity (deg)';
	d.Properties.VariableDescriptions{'dens'} = 'Ganglion cell density (mm^-^2)';
	d.Properties.VariableDescriptions{'densDeg'} = ...
		'Ganglion cell density (deg^-^2)';

	%	Prepare for fitting
	d = sortrows(d, 'ecc'); % sort by eccentricity
	d.eccDegLog = log10(d.eccDeg); d.densDegLog = log10(d.densDeg); % for log fit 
	d.Properties.VariableDescriptions{'eccDegLog'} = 'Log eccentricity (deg)';
	d.Properties.VariableDescriptions{'densDegLog'} = ...
		'Log ganglion cell density (deg^-^2)';

function [d, m] = densGangFun(d, m) % g.c. dens., Wässle (89), vs. fun. ecc.

	% Initialise
	e = d.ecc; % functional eccentricity (deg): es x 1
	n = d.count; % ganglion cell count: es x 1
	
	% Differentiate cell count to find density versus functional eccentricity
	switch 'deriv' % method
		case 'deriv' % derivative
			switch 'all' % retinal region
				case 'near' % near fovea
					i = e < 10; % eccentricities close to fovea
					e = e(i); n = n(i); % select close eccentricities
			end

			% Differentiate
			es = 100; % number of points at which to resample e, for derivative
			eR = linspace(e(1), e(end), es); % resample eA at regular spacing
			n = interp1(e, n, eR); % interpolate to obtain density: es x 1
			dens = gradient(n, eR(2) - eR(1)) ./ (2 * pi * eR); % derivative

			% Set points near fovea
			eR(1) = 0; % fovea
			%	dens(1) = interp1(eR(2: end), dens(2: end), 0, 'pchip'); % extrapolate
			dens(1) = 2 * 8930.1 / sum(m.p.ratGang); % assume twice the MGC density
			eR(2) = .5; % .5 deg
			dens(2) = 34898 / pi; % cell count in 1 deg circle divided by area

		case 'disc' % discrete
			dA = diff(pi * e .^ 2); % annulus areas (deg^2): (es - 1) x 1
			dn = diff(n); % count per annulus: (es - 1) x 1
			dens = (dn ./ dA); % density (deg^-2): (es - 1) x 1
			i = dn >= 2000; % annuli with at least 2000 cells
			e = e(i); es = length(e); % eccentricity (deg): es x 1
			dens = dens(i); % remove low-cell points
			dA = repmat(d(1, :), [es, 1]); % store
			dA.eccDeg = e; dA.densDeg = dens; % store new eccentricities, densities
	end
	
	% Store
	d = repmat(d(1, :), [es, 1]); % one row for each new eccentricity
	d.ecc = eR'; d.densDeg = dens'; % store new eccentricies, densities

function [d, m] = densGangSub(d, m) % split ganglion cells into subpopulations

	if isfield(m.densGangSub, 'z') % set name of variable to be split
		z = m.densGangSub.z; % user specified
	else, z = 'densDeg'; % default
	end
	switch 'fixed' % type of ratio
		case 'fixed' % ratio is constant across eccentricity
			d.type(:) = 'all'; % all ganglion cells
			dMid = d; % replicate
			dMid.type(:) = 'mid'; % midget ganglion cells
			r = sum(m.p.ratGang); % ratio of midget to all ganglion cells
			dMid.(z) = r * dMid.(z); % reduce variable for midget ganglion cells
			d = [d; dMid]; % concatenate
		case 'var' % ratio varies with eccentricity
			dC = ratMidget(d, m); % ratio of midget to all ganglion cells
			r = interp1(dC.eccDeg, dC.ratio, d.eccDeg); % interpolate at ecc. in d
			
			% Calculate the densities of subpopulations
			dAll = d; dAll.type(:) = 'all'; % density of all ganglion cells
			dMid = d; dMid.type(:) = 'mid'; % midget ganglion cells
			dMid.densDeg = r .* dAll.densDeg; % density of midget g.c. (deg^-2)
			dOff = d; dOff.type(:) = 'off'; % off-centre midgets
			dOff.densDeg = m.p.ratOff * dMid.densDeg; % density of off-midgets (deg^-2)
			dOn = d; dOn.type(:) = 'on'; % on-centre midgets
			dOn.densDeg = (1 - m.p.ratOff) * dMid.densDeg; % den. of on-midgets (deg^-2)
			d = [dAll; dMid; dOff; dOn]; % concatenate
	end
	d.Properties.VariableDescriptions{'type'} = 'Cell type'; % describe

function [d, m] = eccFun(d, m) % calc. func. ecc.: McGregor (18), Schein (88)
		
	% Prepare McGregor (18) data
	[d1, m] = eccFunMcG(d, m); % load McGregor (18)
	if isfield(m.eccFunFit, 'mcGregor') % limit to reliable data
		i = d1.eccFun >= m.eccFunFit.mcGregor; % reliable functional eccentricities
		d1 = d1(i, :); % remove outlier data
	end

	% Prepare Schein (88) data
	[d2, m] = offset(d, m); % load Schein (88)
	i = d2.source == 'total'; d2 = d2(i, :); % use total offset (deg): es x 1
	eF = d2.ecc; % functional eccentricity (deg): es x 1
	o = d2.offset; % offset (deg): es x 1
	e = eF + o; % anatomical eccentricity (deg): es x 1
	d2.ecc = e; % anatomical eccentricity (deg): es x 1
	d2.eccFun = eF; % functional eccentricity (deg): es x 1
	if isfield(m.eccFunFit, 'schein') % limit Schein to high eccentricities
		i = d2.ecc >= m.eccFunFit.schein; % high eccentricities
		d2 = d2(i, :); % limit
	end

	% Combine
	d = outerjoin(d1, d2, 'mergeKeys', 1); % concatenate
	d = sortrows(d, {'source', 'ecc'}); % combine and sort

function [d, m] = eccFunFit(d, m) % fit func. ecc.: McGregor (18), Schein (88)

	% Fit functional eccentricity
	switch 'reg' % select method
		case 'hist' % histogram then spline

			% McGregor (18): bin eccentricity and average function ecc. in each bin
			i = d.source == 'mcGregor'; d1 = d(i, :); % select data
			bins = 10; % number of bins in which to average func. ecc.
			[~, edge, i] = histcounts(d1.ecc, bins); % find bin indices of func. ecc.
			binWid = edge(2) - edge(1); % bin width (deg)
			ecc = zeros(1, bins); eccFun = ecc; % allocate storage
			for j = 1: bins % loop over bins
				ecc(j) = edge(j) + .5 * binWid; % ecc. at middle of bin
				eccFun(j) = mean(d1.eccFun(i == j)); % average func. ecc. for this bin
			end
			d1 = d1(1: bins, :); % keep first few lines
			d1.ecc = ecc'; d1.eccFun = eccFun'; % store
		
			%	Interpolate all values with spline
			i = d.source == 'total'; % Schein values at high eccentricities
			d2 = d(i, :); % select data
			ecc = [d1.ecc; d2.ecc]; eccFun = [d1.eccFun; d2.eccFun]; % source data
			n = 50; % number of eccentricities at which to predict fun. ecc.
			eccInt = linspace(ecc(1), ecc(end), n)'; % interpolation ecc.
			eccFunInt = spline(ecc, eccFun, eccInt); % interpolate
			
		case 'reg' % regression

			% Fit polynomial to all data
			eccLog = log10(d.ecc); % for log versus log fit
			d.eccLog = eccLog; % store
			model = fitglm(d, 'poly3', 'predictorVars', 'eccLog', ...
				'responseVar', 'eccFun', 'link', 'log');

			% Rewrite data file
			n = 20; % number of rows
			source = repmat("model", [n, 1]); % data source
			eccLog = linspace(min(eccLog), max(eccLog), n)'; % log ecc. (deg)
			ecc = 10 .^ eccLog; % log anatomical eccentricity (deg)
			d = table(source, eccLog, ecc); % make table
			d.eccFun = predict(model, d); % functional eccentricity (deg)
			d.source = categorical(d.source); % for listing
			d.Properties.VariableDescriptions = {'Data source', ...
				'Log anatomical eccentricity (log deg)', ...
				'Anatomical eccentricity (deg)', 'Functional eccentricity (deg)'};
				% describe variables
			
			% Prepare to add points at high eccentricity
			eccInt = d.ecc; eccFunInt = d.eccFun; % match name with other methods
			dC = d(1: 2, :); d = [d; dC]; % add two rows to data file

		case 'spline' % regression then spline

			% Fit polynomial to McGregor (18)
			i = d.source == 'mcGregor'; d1 = d(i, :); % select data
			model = fitlm(d1, 'poly3', 'predictorVars', 'ecc', ...
				'responseVar', 'eccFun');
			eccEnd = d1.ecc(end); % save final eccentricity (deg)
			n = 20; % number of eccentricities at which to predict fun. ecc.
			d1 = d1(1: n, :); % keep first few lines
			eccBeg = 1; % assume x-intercept (deg)
			ecc = linspace(eccBeg, eccEnd, n)'; d1.ecc = ecc; % subsample ecc.
			eccFun = predict(model, d1); % calculate model predictions
			ecc(1) = eccBeg; eccFun(1) = 0; % intercept of prediction with ecc. axis
			d1 = d1(1: length(ecc), :); % trim to non-negative fun. ecc.
			d1.ecc = ecc; d1.eccFun = eccFun; % store
		
			%	Interpolate all values with spline
			i = d.source == 'total'; % Schein values at high eccentricities
			d2 = d(i, :); % select data
			ecc = [d1.ecc; d2.ecc]; eccFun = [d1.eccFun; d2.eccFun]; % source data
			n = 50; % number of eccentricities at which to predict fun. ecc.
			eccInt = linspace(eccBeg, ecc(end), n)'; % interpolation ecc.
			eccFunInt = spline(ecc, eccFun, eccInt); % interpolate

	end
	
	% Add points at high eccentricity
	u = m.p.magRet * 2.88; % eccentricity at which offset = 0 (deg):
		% from Schein (88) Figure 14 Nasal and Temporal
	eccInt = [eccInt; u; 90]; eccFunInt = [eccFunInt; u; 90];
		% add end-points (deg): es x 1
	d = d(1: n + 2, :); d.ecc = eccInt; d.eccFun = eccFunInt; % store
	d.eccLog = log10(d.ecc); % for log versus log fit

function [d, m] = eccFunMcG(~, m) % functional vs. anat. ecc.: McGregor (18)

	folder = [m.folder, filesep, 'McGregor (18)']; % data folder
	d = m.readFile(folder); % create table from folder
	d.eccMm = d.eccUm / 1000; d.eccFunMm = d.eccFunUm / 1000; % convert to mm
	d = renamevars(d, 'label', 'source'); % rename label
	d = sortrows(d, {'source', 'eccMm'}); % sort by source, eccentricity
	d.source = categorical(d.source); % make it categorical
	d.ecc = m.p.magRet * d.eccMm; d.eccFun = m.p.magRet * d.eccFunMm; % mm to deg
	d.Properties.VariableDescriptions{'source'} = 'Source';
	d.Properties.VariableDescriptions{'ecc'} = 'Eccentricity (deg)';
	d.Properties.VariableDescriptions{'eccFun'} = 'Functional eccentricity (deg)';

function [d, m] = fitCone(d, m) % fit cone response: Baudin (19)

	model = fitnlm(d.freq, d.sens, m.cone, m.fitnlm.beta0, ...
		'coefficientNames', {'sens0', 'tau'});
	%	'predictorVars', {'freq'}, 'responseVar', 'sens');
	m.model{m.group} = model; % store model

function [d, m] = fitDogS(d, m) % fit inverse transform of DOG to ROG inv. tran.

	% Prepare data for fitting
	x = d.x'; % predictor variable is location (deg): rs x 1
	x = repmat(x, [2, 1]); % duplicate for real then imag. resp.: (2 * rs) x 1
	r = d.respS .'; % ROG response in spatial domain, complex: rs x 1
	rs = length(r); % number of responses
	y = [real(r); imag(r)]; % real then imag. responses: (2 * rs) x 1

	%	Fit
	b = m.fitDogS.coef; bFix = b; % model coefficients
	i = ~ contains(m.fitDogS.name, m.fitDogS.fix); % indices of variable coef.
	bFix(i) = nan; % replace variable coefficients with nan
	fun = @ (b, x) m.dogS(b, x, bFix); % DOG
	bOpt = b(i); % initial estimates for coefficients
	names = m.fitDogS.name(i); % names of variable coefficients
	model = fitnlm(x, y, fun, bOpt, 'coefficientNames', names);
	m.model{m.group} = model; % store model

	%	Predict response from fitted model, and calculate polar components
	r = predict(model, x); % predict responses: (2 * rs) x 1
	rR = r(1: rs); rI = r(rs + 1: end); % real and imaginary parts
	r = complex(rR, rI); % combine to obtain complex response
	d.respSReal = rR'; % real part: 1 x rs
	d.respSAmp = abs(r)'; % amplitude (Hz): 1 x rs
	p = angle(r); % phase (radians)
	d.respSPhase = (180 / pi) * p'; % phase (deg)
	d.resp = r .'; % response, complex

function [d, m] = fitPulse(d, m) % fit pulse response: Lee (94), Sinha (17)

	% Initialise
	switch m.resp.source % select data source
		case 'Lee'
			dur = unique(d.dur)'; % pulse duration (s): 1 x durs
			durs = length(dur); % number of durations
			cs = 1; % number of cells
			ts = size(d, 1) / durs; % number of times per duration
			t = d.time(1: ts); % time (s): ts x 1
			durSamp = .002; % sample duration (s)
			t = t - .5 * durSamp; % starting time is 0: ts x 1
			per = t(end) + durSamp; % simulation period (s)
		case 'Sinha'
			dur = .01; durs = 1; % pulse duration (s), number of durations
			c = unique(d.type); % cell: 1 x cs
			cs = length(c); % number of cells
			ts = size(d, 1) / cs; % number of samples per response
			per = .2; % simulation period (s)
	end

	%	Calculate frequencies and response
	nH = floor(.5 * ts); % half the number, converted to integer
	n = 2 * nH; % make sure that the number of samples is even
	fund = 1 / per; % fundamental frequency (Hz)
	ind = linspace(- nH, nH, n + 1); % sample indices
	ind(end) = []; % make the sequence open-ended
	freqT = fund * ind'; % temporal frequencies (Hz): ts x 1
	freqT = fftshift(freqT); % change to fft format
	freqS = 0 * freqT; % spatial frequencies (cycles/deg): ts x 1
	x = [freqS, freqT]; % frequencies: ts x 2
	x = repmat(x, [durs * cs, 1]); % one replicate for each response: ys x 1
	y = d.resp; % empirical pulse response (Hz): ys x 1

	% Prepare coefficients for fitting function
	b = m.fitRog.coef; % initial values of all coefficients
	bFix = b; % prepare for fixed coefficients
	i = ~ contains(m.fitRog.name, m.fitRog.fix); % indices of variable coef.
	bFix(i) = nan; % replace variable coefficients with nan
	bOpt = b(i); % variable coefficients
	name = m.fitRog.name(i); % names of variable coefficients

	% Calculate pulse transform
	f = 2 * pi * freqT; % temporal frequencies (rad/s): ts x 1
	a = f .* dur / 2; % argument for pulse transform (rad): ts x durs
	p = 2 * sin(a) ./ f; % pulse transform: ts x durs
	p = exp(- 1i * a) .* p; % shift so that pulse starts at 0: ts x durs
	i = f == 0; % indices of rows for which frequency is zero
	p(i, :) = dur(i, :); % replace nan
	p = repmat(p, [1, cs]); % one replicate for each response: ys x 1

	% Fit model and predict responses from it
	fun = @(b, x)m.rog(b, x, bFix, p); % fitting function
	model = fitnlm(x, y, fun, bOpt, 'coefficientNames', name); % fit model
	m.model{m.group} = model; % store model
	d.resp = predict(model, x); % predicted response (Hz): ts x 1

function [d, m] = fitRog(d, m) % fit DOG/ROG to freq. response: multiple authors

	% Prepare data for fitting
	x = [d.freqS, d.freqT]; % predictor variables: rs x 1
	x = repmat(x, [2, 1]); % duplicate for real then imag. resp.: (2 * rs) x 1
	r = d.resp; % responses, complex: rs x 1
	rs = size(r, 1); % number of responses
	y = [real(r); imag(r)]; % real then imag. responses: (2 * rs) x 1

	% Set up for fitting function
	b = m.fitRog.coef; bFix = b; % model coefficients
	i = ~ contains(m.fitRog.name, m.fitRog.fix); % indices of variable coef.
	bFix(i) = nan; % replace variable coefficients with nan
	switch m.fitRog.funFun % difference or ratio of Gaussians?
		case 'dog', fun = @(b, x)m.dog(b, x, bFix); % difference
		case 'rog', fun = @(b, x)m.rog(b, x, bFix); % ratio
	end
	bOpt = b(i); % coefficients for optimisation
	names = m.fitRog.name(i); % names of variable coefficients

	% Fit model
	switch 'double' % call fitnlm
		case 'double' % inputs are doubles
			model = fitnlm(x, y, fun, bOpt, 'coefficientNames', names);
		case 'table' % input is table; fails because x provided to fun is (1: n)'
			model = fitnlm(d, fun, bOpt, 'predictorVars', {'freqS', 'freqT'}, ...
				'responseVar', 'resp', 'coefficientNames', names);
	end
	m.model{m.group} = model; % store model

	%	Predict response from fitted model, and calculate polar components
	r = predict(model, x); % predict responses: (2 * rs) x 1
	rR = r(1: rs); rI = r(rs + 1: end); % real and imaginary parts
	r = complex(rR, rI); % combine to obtain complex response
	d.amp = abs(r); % amplitude (Hz)
	p = angle(r); % phase (radians)
	d.phase = (180 / pi) * p; % phase (deg)
	d.resp = r; % response, complex: rs x 1

function [d, m] = makeTable(~, m) % make table of temporal freq. for respS

	f = m.makeTable.freqT(:); % temporal frequencies supplied by user
	d = table(f, 'variableNames', {'freqT'}); % make table
	d.Properties.VariableDescriptions = {'Temporal frequency (Hz)'};

function [d, m] = mtf(~, m) % compile modulation transfer function

	% Parameters for modulation transfer function: Navarro (93)
	ecc = [	0			5			10		20		30		40		50		60]'; % eccentricity (deg)
	a =		[	.172	.245	.245	.328	.606	.82		.93		1.89]';
	b =		[	.037	.041	.041	.038	.064	.064	.059	.108]';
	c =		[	.22		.2		.2		.14		.12		.09		.067	.05]';
	
	% Calculate MTF and store
	es = length(ecc); % number of eccentricity
	freq = linspace(0, 60); % spatial frequency (cycles/deg): 1 x fs
	mtf = (1 - c) .* exp(- a .* freq) + c .* exp(- b .* freq); % MTF: es x fs
	freq = repmat(freq, [es, 1]); % spatial frequency (cycles/deg): es x fs
	d = table(ecc, a, b, c, freq, mtf); % store as table: es x 6
	d.Properties.VariableDescriptions{'ecc'} = 'Eccentricity (deg)';
	d.Properties.VariableDescriptions{'freq'} = 'Spatial frequency (cycles/deg)';
	d.Properties.VariableDescriptions{'mtf'} = 'Modulation transfer function';

function [d, m] = offset(~, m) % ganglion cell lateral offset: Schein (88)

	folder = [m.folder, filesep, 'Schein (88)']; % data folder
	d = m.readFile(folder); % create table from folder
	d = renamevars(d, 'label', 'source'); % rename label
	d = sortrows(d, {'source', 'eccMm'}); % sort by source, eccentricity
	d.source = categorical(d.source); % make it categorical
	d.ecc = m.p.magRet * d.eccMm; d.offset = m.p.magRet * d.offsetMm; % mm to deg
	d.Properties.VariableDescriptions{'source'} = 'Source';
	d.Properties.VariableDescriptions{'ecc'} = 'Eccentricity (deg)';
	d.Properties.VariableDescriptions{'offset'} = 'Displacement (deg)';

function [d, m] = predRog(d, m) % predict responses with ROG model

	%	Predict responses from fitted model, and calculate polar components
	model = m.model{m.group}; % fitted model
	x = [d.freqS, d.freqT]; % frequencies
	rs = size(x, 1); % number of rows with unique predictors
	x = repmat(x, [2, 1]); % duplicate for real then imag. resp.: (2 * rs) x 1
	r = predict(model, x); % predict responses: (2 * rs) x 1
	rR = r(1: rs); rI = r(rs + 1: end); % real and imaginary parts
	r = complex(rR, rI); % combine to obtain complex response
	d.amp = abs(r); % amplitude (Hz)
	p = angle(r); % phase (radians)
	d.phase = (180 / pi) * p; % phase (deg)
	d.resp = r; % response, complex: rs x 1

function [d, m] = prep(d, m) % prepare for plot

	% Set independent and dependent variables, and list preparations
	if isfield(m.prep, 'x'), nX = m.prep.x; x = d.(nX); % independent variable
	elseif isfield(m, 'x'), nX = m.x; x = d.(nX); end

	%	Prepare
	n = fieldnames(m.prep)'; % names, including preparations (cell): 1 x ps
	for p = string(n) % loop over preparations
		switch p % choose preparation
			case 'unwrap' % avoid 360 deg jumps in phase
				n = m.prep.unwrap; % name of phase variable
				phase = (pi / 180) * d.(n); % phase (radians)
				phase = unwrap(phase); % avoid 360 deg jumps in phase
				d.(n) = (180 / pi) * phase; % store phase (deg)
			case 'zero' % replace 0 on log x-axis by positive value *** check ***
				i = x == 0; % rows for which x is 0
				x(i) = m.prep.zero; % replace
				d.(nX) = x; % store
		end
	end

function [d, m] = prepPulse(d, m) % prepare pulse resp. for fitting: Sinha (17)

	%	Make time samples evenly spaced
	t = -.05: .002: .15; % sample times (s)
	ts = length(t); % number of samples
	r = interp1(d.time, d.resp, t); % interpolate with fixed sample interval
	d = repmat(d(1, :), [ts, 1]); % replicate first line for new length
	d.time = t'; d.resp = r'; % store

	%	Remove data before time 0
	dC = d(end, :); % save last line for extrapolation
	i = d.time == 0; i = find(i); % index for time 0
	d = d(i: end - 1, :); % use an open-ended interval starting at time 0

	% Extrapolate from last response towards steady state
	tSamp = d.time(2) - d.time(1); % sample interal (s)
	t = dC.time: tSamp: .2; % extrapolation time (s)
	ts = length(t); % number of points in extrapolatin
	dC = repmat(dC, [ts, 1]); % extrapolation table
	r = dC.resp(1) * exp(- (t - dC.time(1)) / .01); % extrapolated response
	dC.time = t'; dC.resp = r'; % store
	dC(end, :) = []; % open-ended interval
	d = [d; dC]; % extend table

function [d, m] = profHor(~, m) % horizontal cell receptive field profiles

	ecc = 30; % eccentricity (deg)
	d = m.radHorConv(ecc, m); % calculate receptive fields, cell: es x 1
	f = d.fieldAtt{1}; % attenuated receptive field: xs x ys x ls
	x = d.x; xs = length(x); % x values (deg), number of x values
	loc = d.loc{1}; % cell locations (deg): ls x 2
	ls = size(loc, 1); % number of cells
	fC = zeros(xs, ls); % allocate storage
	for i = 1: ls % loop over cells
		y = loc(i, :); % location of current cell (deg): 1 x 2
		j = knnsearch(x', y(2)); % index y value for current cell
		fC(:, i) = f(:, j, i); % section through centre: xs x 1 x ls
	end
	f = fC; % centred section along x axis: xs x ls
	d.field = shiftdim(f, -1); % store

function [d, m] = psf(d, m) % calculate point spread function

	% Calculate point spread function as inverse Fourier transform of MTF
	es = size(d, 1); % number of eccentricities
	a = d.a; b = d.b; c = d.c; % coefficients: es x 1
	w = .2; xs = 100; % location range (deg), number of locations
	x = linspace(-.5, .5, xs + 1) * w; % location (deg): 1 x (xs + 1)
	x(end) = []; % make it open-ended: 1 x xs
	y = 2 * a .* (1 - c) ./ (a .^ 2 + (2 * pi * x) .^ 2) + ...
		2 * b .* c ./ (b .^ 2 + (2 * pi * x) .^ 2); % PSF: es x xs

	% Fourier transform PSF to check on solution
	f = linspace(- .5, .5, xs + 1) * (xs / w); % spatial freq. (cycles/deg)
	f(end) = []; % make it open-ended: 1 x xs
	f = fftshift(f); % shift negative frequencies to end: 1 x xs
	z = fft(y, [], 2) ./ sum(y, 2); % Fourier transform: es x xs
		% *** fix ***
	z = abs(z); % MTF: es x xs

	% Store in table
	d.loc = repmat(x, [es, 1]); % location (deg): es x xs
	d.psf = y; % PSF: es x xs
	d.Properties.VariableDescriptions{'loc'} = 'Location (deg)'; % describe
	d.Properties.VariableDescriptions{'psf'} = 'Point spread function';
	d.freq = repmat(f, [es, 1]); % spatial frequency (cycles/deg): es x xs
	d.mtf = z; % MTF: es x xs
	d.Properties.VariableDescriptions{'freq'} = 'Spatial frequency (cycles/deg)';
	d.Properties.VariableDescriptions{'mtf'} = 'Modulation transfer function';

function [d, m] = rad(d, m) % calculate radii contibuting to centre, surround

	% Calculate optical point spread function radius
	d = radOpt(d, m); % load optical data
	switch 'fun' % choose method
		case 'emp' % use eccentricities in optics empirical data
			d = d(d.ecc >= 0, :); % retain only temporal data
			eF = d.ecc; % eccentricity (deg)
			d.optics = d.radius; % PSF radius (deg)
		case 'fun' % calculate eccentricities for fitted data
			l = log10([.006, 40]); % log eccentricity limits (deg)
			es = 30; % number of eccentricities
			eF = logspace(l(1), l(2), es)'; % eccentricity (deg): es x 1
			r = polyval(m.p.radOptCoef, eF);
			d = repmat(d(1, :), [es, 1]); % repeat first row of table
			d.ecc = eF; d.optics = r; % store
	end

	% Calculate remaining radii
	d.denGang = polyval(m.p.radGangCoef, eF); % ganglion cell dend. radius (deg)
	d.cenConv = sqrt(d.optics .^ 2 + d.denGang .^ 2); % centre radius (deg)
	if m.rad.horExp % horizontal cells have exponential receptive fields
		c =  m.p.radBackCoef; % regression coefficients for f'back convergence fun.
		d.fieldHor = 10 .^ (c(1) + c(2) * normcdf(eF, c(3), c(4))); % conv. f. (deg)
		d.surConv = sqrt(d.optics .^ 2 + d.fieldHor .^ 2); % sur. radius (deg)
	else
		d.fieldHor = exp(polyval(m.p.radHorCoef, eF)); % rec. field radius (deg)
		d.surConv = sqrt(d.optics .^ 2 + 2 * d.fieldHor .^ 2); % sur. radius (deg)
	end

	% Turn lateral array into vertical
	s = {'optics', 'denGang', 'cenConv', 'fieldHor', 'surConv'};
		% sources of radius (deg)
	ss = length(s); % number of sources
	dC = cell(ss, 1); % allocate storage
	for i = 1: ss % loop over sources
		dCC = d; % replicate existing table
		sC = s{i}; % current source
		dCC.source(:) = sC; dCC.radius = d.(sC); % current source and radius
		dC{i} = dCC; % store
	end
	d = vertcat(dC{:}); % concatenate into vertical array
	d = removevars(d, ['rrfLog', 'rrf', s]); % variables not needed
	d.source = categorical(d.source); % make it categorical, for listing purposes
	d.radiusMin = 60 * d.radius; % convert to minutes, for plotting
	d.Properties.VariableDescriptions = {'Source of radius', ...
		'Eccentricity (deg)', 'Radius (deg)', 'Radius (min)'};

function [d, m] = radCen(~, m) % compile ganglion cell centre radius data

	% Croner (95)
	load([m.folder, filesep, 'Croner (95) rad/parvo_cenrad_ecc.mat'], 'b');
	ecc = b(:, 1); % eccentricity (deg): es x 1
	es = size(ecc, 1); % number of eccentricities
	radius = b(:, 2); % ganglion cell centre radius (deg): es x 1
	source = repmat("Croner", [es, 1]); % data source
	d1 = table(source, ecc, radius); % store

	% Lee (98)
	load([m.folder, filesep, 'Lee (98)/CenLeeexp.mat'], 'CenLeeexp');
	d2 = CenLeeexp; % ecc. and radius: es x 2
	ecc = d2(:, 1); % eccentricity (deg): es x 1
	es = size(ecc, 1); % number of eccentricities
	dev = 10 .^ d2(:, 2); % standard deviation (min): es x 1
	radius = sqrt(2) * dev / 60; % gang. cell centre radius (deg): es x 1
	source = repmat("Lee", [es, 1]); % data source
	d2 = table(source, ecc, radius); % store
	
	%	Godat (22)
	d3 = radGangField(d2, m); % read data
	i = d3.type == 'centre'; d3 = d3(i, :); % keep only centre data
	d3 = renamevars(d3, 'type', 'source'); % prepare to set source
	d3.source(:) = "Godat"; d3.source = string(d3.source); % set source
	d3 = removevars(d3, {'eccLog', 'radLog'}); % remove extraneous variables

	% Concatenate
	d = [d1; d2; d3]; % concatenate
	d.source = categorical(d.source); % for listing
	d.radiusLog = log10(d.radius); % add logarithm for fitting purposes
	d.radiusMin = 60 * d.radius; % radius (min) for display purposes
	d.dev(:) = nan; d.dev(end - es + 1: end) = dev; % s.d. (min)
	d.Properties.VariableDescriptions = {'Source', ...
		'Eccentricity (deg)', 'Centre radius (deg)', ...
		'Log of centre radius (deg)', 'Centre radius (min)', ...
		'Standard deviation (min)'};

function [d, m] = radCone(~, m) % cone inner segment radius: Packer (89) diam

	% Generate a table for each quadrant
	folder = [m.folder, filesep, 'Packer (89) rad']; % data folder
	f = dir([folder, '/*.mat']); % mat-file list: fs x 1
	fs = length(f); % number of files
	dC = cell(fs, 1); % allocate storage			
	for i = 1: fs % loop over files

		% Collect data
		fC = f(i).name; % full name of current file
		q = fC(1: 2); % name of current quadrant: 2 x 1
		quad = categorical(string(q)); % store
		d = load([folder, '/', fC]); % contents of file
		d = d.(q); % eccentricity and diameter (mm, um): ds x 2
		ds = size(d, 1); % number of rows

		% Store in table
		quad = repmat(quad, [ds, 1]); % match number of rows
		dC{i} = table(quad, d(:, 1), d(:, 2), 'variableNames', ...
			{'quad', 'ecc (mm)', 'diam (um)'});

	end

	% Combine and store
	d = vertcat(dC{:}); % combine data tables
	d = sortrows(d, {'quad', 'ecc (mm)'}); % sort rows
	d.Properties.VariableDescriptions = ...
		{'Retinal quadrant', 'Eccentricity (mm)', 'Cone diameter (um)'};

function [d, m] = radHor(d, m) % hor. cell radius: Wässle (89), Packer (02)

	% Dendrite radius: Wässle (89) Horizontal
	[d, m] = radHorDend(d, m); % obtain data
	d = renamevars(d, 'type', 'source'); % rename label
	i = d.source == "H1"; d = d(i, :); % keep only H1 cells
	d.source(:) = "dend"; % source of data
	d = d(d.eccFun <= 16, :); % field only at low eccentricities
	r = d.area; % dendrite area (mm^2)
	r = sqrt(r / pi); % radius (mm): es x 1
	d.radius = m.p.magRet * r; % radius (deg): es x 1
	if isfield(m.radHor, 'double') % account for double-pass through dendrites
		d.radius = sqrt(2) * d.radius; % sum of squares rule
	end
	d.Properties.VariableDescriptions{'radius'} = 'Radius (deg)';
	d1 = d; % store
	
	% Field radius: Packer (02)
	[d, m] = radHorField(d, m); % obtain data
	d = m.radHorComb(d, m); % bin the data and calc. one radius for each ecc.
	d = outerjoin(d1, d, 'mergeKeys', 1); % concatenate

function [d, m] = radHorDend(~, m) % hor. cell dendrite area: Wässle (89) Hor

	% Compile data into table
	folder = [m.folder, filesep, 'Wässle (89) Horizontal']; % data folder
	d = m.readFile(folder); % create table from folder
	d = renamevars(d, 'label', 'type'); % rename label
	d.type = categorical(d.type); % make it categorical
	d.eccDeg = m.p.magRet * d.ecc; % eccentricity (deg)
	d.area = 10 .^ d.areaLog; % store the anti-log
	d.Properties.VariableDescriptions{'type'} = 'Cell type'; % describe variable
	d.Properties.VariableDescriptions{'eccDeg'} = 'Eccentricity (deg)';
	d.Properties.VariableDescriptions{'areaLog'} = 'Log area (mm^2)';
	d.Properties.VariableDescriptions{'area'} = 'Area (mm^2)';

	% Add functional eccentricity
	s = load([m.folder, filesep, 'Functional eccentricity'], 'd'); dE = s.d;
	d.eccFun = interp1(dE.ecc, dE.eccFun, d.eccDeg); % interpolate on fun. ecc.
	d.Properties.VariableDescriptions{'eccFun'} = 'Functional eccentricity (deg)';

function [d, m] = radHorDiam(~, m) % hor. cell field diam. limits: Packer (02)

	d = m.radHorLim(m); % limits
	d.source(:) = "field"; % add source
	
function [d, m] = radHorField(~, m) % horizontal cell field radius: Packer (02)

	% Compile data into table
	folder = [m.folder, filesep, 'Packer (02)']; % data folder
	d = m.readFile(folder); % create table from folder
	d = renamevars(d, 'label', 'source'); % rename label
	d.source = categorical(d.source); % make it categorical
	d.eccFun = m.p.magRet * d.ecc; % functional eccentricity: (deg)
	r = d.diam; % diameter (um)
	%{
	r = .5 * m.p.magRet * .001 * r; % radius (deg)
	d.radius = r / sqrt(log(10)); % convert from .1 * peak to 1 / e (deg)
	%}
	r = .5 * .001 * r; % radius (mm)
	r = m.p.magRet * r; % radius (deg)
	d.radius = r / log(10); % convert from .1 * peak to (1 / e) * peak (deg)
	d.Properties.VariableDescriptions{'source'} = 'Source of radius';
	d.Properties.VariableDescriptions{'eccFun'} = 'Functional eccentricity (deg)';
	d.Properties.VariableDescriptions{'radius'} = 'Field radius (deg)';

function [d, m] = radHorFit(d, m) % fit h.c. radius: Wässle (89), Packer (02)

	%	Calculate regression of radius on eccentricity
	d.logRadius = log10(d.radius); % for fitting: es x 1
	fun = @(b, x) (b(1) + b(2) * normcdf(x, b(3), b(4))); % function to fit
	bInit = [-1.3, 1.4, 20, 7]; % initial values of regression coefficients
	model = fitnlm(d.eccFun, d.logRadius, fun, bInit); % fit
	m.model{m.group} = model; % store
	d.logRadius = predict(model, d.eccFun); % calculate model: es x 1
	d.radius = 10 .^ d.logRadius; % calculate radius (deg): es x 1

function [d, m] = radGang(~, m) % ganglion cell denritic radius: Watanabe (89)

	% Load radius data and convert eccentricity
	folder = [m.folder, filesep, 'Watanabe (89)']; % data folder
	d = m.readFile(folder); % create table from folder
	d.label = []; % remove this variable
	d.ecc = m.p.magRet * d.eccMm; % eccentricity (deg): es x 1
	s = load([m.folder, filesep, 'Functional eccentricity'], 'd'); dE = s.d;
	d.eccFun = interp1(dE.ecc, dE.eccFun, d.ecc); % convert to functional ecc.
	d = sortrows(d, 'eccFun'); % sort on eccentricity
	
	% Store
	d.diam = 10 .^ d.diamLog; % dendritic field diameter (um): es x 1
	d.radius = (.5 * d.diam / 1000) * m.p.magRet; % radius (deg): es x 1
	d.Properties.VariableDescriptions = ...
		{'Eccentricity (mm)', 'Log of dendritic field diameter (um)', ...
			'Eccentricity (deg)', 'Functional eccentricity (deg)', ...
			'Dendritic field diameter (um)', 'Dendritic field radius (deg)'};

function [d, m] = radGangField(~, m) % g.c. centre, surround radius: Godat (22)

	% Load data and create a table
	folder = [m.folder, filesep, 'Godat (22)']; % data folder
	d = m.readFile(folder); % create table from folder
	d = renamevars(d, 'label', 'type'); % remove this variable
	d.type = categorical(d.type); % make it categorical
	d.ecc = 10 .^ d.eccLog; d.radius = 10 .^ d.radLog; % take anti-logs
	d = sortrows(d, {'type', 'ecc'}); % sort by type, eccentricity
	d.Properties.VariableDescriptions{'type'} = 'Mechanism';
	d.Properties.VariableDescriptions{'ecc'} = 'Eccentricity (deg)';
	d.Properties.VariableDescriptions{'radius'} = 'Mechanism radius (deg)';

function [d, m] = radOpt(~, m) % point spread function radius: Navarro (93)

	% Generate a table from the data
	folder = [m.folder, filesep, 'Navarro (93)']; % data folder
	d = m.readFile(folder); % create table from folder
	d = sortrows(d, {'label', 'ecc'}); % sort on version, eccentricity
	d = renamevars(d, 'label', 'source'); % rename label
	d.source(:) = "optics"; d.source = categorical(d.source); % clean label

	% Clean data
	r = 10 .^ d.rrfLog; % retinal resolution function (min)
	es = .5 * size(d, 1); % number of eccentricities
	r = mean([r(1: es), r(es + (1: es))], 2); % mean over versions
	d = d(1: es, :); % keep only mean
	d.ecc = [-60, -50, -40, -30, -20, -10, -5, 0, 5, 10, 20, 30, 40, 50, 60]';
		% clean eccentricities (deg)
	d.rrf = r; % retinal resolution function, given as	full width
		%	at half-height of Gaussian (min): es x 1

	% Convert rrf into radius
	r = r / (2 * sqrt(log(2))); % radius of point spread function (min): es x 1
	r = r / 60; % radius of point spread function (deg): es x 1
	d.radius = m.p.ratOpt * r; % convert human to macaque (deg): Harwerth (85)
	d.Properties.VariableDescriptions{'source'} = 'Source of data';
	d.Properties.VariableDescriptions{'rrf'} = 'Retinal resolution function (min)';
	d.Properties.VariableDescriptions{'radius'} = 'PSF radius (deg)';

function [d, m] = radSur(~, m) % g.c. sur.: Croner (95), Lee (98), Godat (22)

	% Compile data into table: Croner (95)
	folder = [m.folder, filesep, 'Croner (95) rad']; % data folder
	d = m.readFile(folder); % create table from folder
	d = renamevars(d, 'label', 'source'); % rename label
	d.source = categorical(d.source); % make it categorical
	d.ecc = 10 .^ d.eccLog; % eccentricity (deg)
	d.radius = 10 .^ d.radiusLog; % surround radius (deg)
	d.Properties.VariableDescriptions{'source'} = 'Source of radius'; % describe
	d.Properties.VariableDescriptions{'ecc'} = 'Eccentricity (deg)';
	d.Properties.VariableDescriptions{'radius'} = 'Surround radius (deg)';

	% Lee (98): load both centre and surround data to find surround eccentricities
	folder = [m.folder, filesep, 'Lee (98)']; % data folder
	dC = m.readFile(folder); % centre s.d. data
	dC = renamevars(dC, 'label', 'source'); % rename label
	dC.source = categorical(dC.source); % make it categorical
	dC.radCen = 10 .^ dC.radLog; % centre s.d. (min)
	dC = sortrows(dC, 'radCen'); % sort to match with surround data
	load([folder, filesep, 'CenSurLee'], 'CenSurLee'); % load surround
	dS = CenSurLee; % surround data
	dS = table(dS(:, 1), dS(:, 2), 'variableNames', {'radCen', 'radSur'}); % table
	dS = sortrows(dS, 'radCen'); % sort to match with centre data
	dC.radius = sqrt(2) * dS.radSur / 60; % radius (deg)
	dC.source(:) = "Lee"; % show source
	d = outerjoin(d, dC, 'mergeKeys', 1); % combine Croner, Lee
	
	% Godat (22)
	dC = radGangField(d, m); % read data
	i = dC.type == 'surround'; dC = dC(i, :); % keep only surround data
	dC = renamevars(dC, 'type', 'source'); % prepare to set source
	dC.source(:) = "Godat"; dC.source = categorical(dC.source); % set source
	d = outerjoin(d, dC, 'mergeKeys', 1); % add Godat

	%	Combine
	d.radiusMin = 60 * d.radius; % convert to minutes, for plotting

function [d, m] = ratCount(d, m) % ratio of cone to ganglion cell count

	eInc = .5; % eccentricity increment (deg)
	e = (0: eInc: d.ecc(end))'; % eccentricities at which to interpolate
	cG = interp1(d.ecc, d.count, e); % ganglion cell count
	dC = densCone(d, m); dC = countCone(dC, m); % cone count
	cC = interp1(dC.eccDeg, dC.count, e); % ganglion cell count
	r = cC ./ cG; % ratio of cone to ganglion cell count
	d = repmat(d(1, :), [length(e), 1]); % resize table
	d.ecc = e; d.ratio = r; % store eccentricity, ratio
	d.densDeg = cG ./ (pi * e .^ 2); % cumulative ganglion cell density (deg^-2)
	d.Properties.VariableDescriptions{'ratio'} = ...
		'Ratio of cone to ganglion cell count'; % describe

function [d, m] = ratDens(d, m) % ratio of cone to ganglion cell density

	dC = densCone(d, m); % load cone density
	i = dC.quad == 'temporal'; dC = dC(i, :); % keep only temporal quadrant
	densC = interp1(dC.eccDeg, dC.densDeg, d.ecc); % cone dens. at ecc. in d
	d.ratio = densC ./ d.densDeg; % ratio of cone to ganglion cell density
	d.Properties.VariableDescriptions{'ratio'} = ...
		'Ratio of cone to ganglion cell density'; % describe
	
function [d, m] = ratMidget(~, m) % mid. / all g.c.: % Dacey (94), Grünert (93)

	% *** unused: delete ***
	folder = [m.folder, filesep, 'Dacey (94)']; % data folder
	d = m.readFile(folder); % log ratio versus eccentricity
	d.eccDeg = m.p.magRet * d.ecc; % eccentricity (deg)
	r = 10 .^ d.ratioLog; % ratio of midget ganglion cells to all g.c.
	d.ratioPer = min(100 - 6.5, r); ...
		% Grünert et al. found that 5-8% of fovel ganglion cells were parasol
	d.ratio = .01 * d.ratioPer; % midget/ all g.c.
	d.Properties.VariableDescriptions{'eccDeg'} = 'Eccentricity (deg)';
	d.Properties.VariableDescriptions{'ratioPer'} = ...
		'Midget / all ganglion cells (%)';
	d.Properties.VariableDescriptions{'ratio'} = ...
		'Midget / all ganglion cells';

function [d, m] = refine(d, m) % refine predictor variable for model plotting

	name = m.x; % predictor name
	x = d.(name); % existing values: xs x 1
	x = unique(x); x = sort(x); % sort into ascending order: xs x 1
	if x(1) == 0 && isfield(m.refine, 'zero') % replace zero?
		x(1) = m.refine.zero; % replace zero with user-specified value
	end
	xs = m.refine.n; % number of values to generate
	x = logspace(log10(x(1)), log10(x(end)), xs)'; % refine: xs x 1
	d = repmat(d(1, :), [xs, 1]); % replicate first row of table
	d.(name) = x; % store

function [d, m] = resp(~, m) % ganglion cell frequency resp.: multiple authors

	switch m.resp.source % data source
		case 'Benardete' % Benardete (97)

			% Gain data
			folder = [m.folder, filesep, 'Benardete (97) gain'];
			d = m.readFile(folder); % gain versus spatial frequency
			d = renamevars(d, 'label', 'freqT'); % temporal frequency
			d.source(:) = categorical("Benardete"); % add source
			d = sortrows(d, {'freqT', 'freqS'}); % sort rows

			% Fix spatial frequencies
			f = [0, .143, .429, .715, 1.57, 2.15, 3.29, 4.43, 6.72, 9.01, 18.16];
				% spatial frequencies (cycles/deg)
			fC = - fliplr(f); fC(end) = []; f = [fC, f]; % add negative frequencies
			d.freqS = repmat(f', [3, 1]); % replace digitised values: fs x 1

			% Phase data and response
			folder = [m.folder, filesep, 'Benardete (97) phase'];
			dC = m.readFile(folder); % phase versus spatial frequency
			dC = renamevars(dC, 'label', 'freqT'); % temporal frequency
			dC = sortrows(dC, {'freqT', 'freqS'}); % sort rows
			p = dC.phaseSemi; % response phase (semicycles)
			d.phaseSemi = p; % store
			p = 180 * p; % response phase (deg): fs x 1

			%	Correct phase for stimulus displacement, and store response
			if isfield(m.resp, 'disp') % displacement is defined, so correct phase
				p = p + 360 * m.resp.disp * d.freqS; % correct for displacement: fs x 1
				d.freqS = abs(d.freqS); % no need for negative frequencies
			end
			i =  d.freqT > 2.2; % shift phase to reflect func. delays for higher freq.
			p(i) = p(i) - 360; % delay by a cycle
			d.phase = p; % response phase (deg)
			d.resp = d.amp .* exp(1i * (pi / 180) * d.phase);
				% complex response (Hz/unit-contrast)

			%	Finalise
			d.Properties.VariableDescriptions{'amp'} = 'Gain (Hz/unit-contrast)';
			d.Properties.VariableDescriptions{'phase'} = 'Response phase (deg)';
			
		case 'Croner' % Croner (95)
			folder = [m.folder, filesep, 'Croner (95) resp'];
			d = m.readFile(folder); % impulse rate versus spatial frequency
			d = renamevars(d, 'label', 'source'); % rename for consistency
			d.source = categorical(d.source); % for pretty printing
			d.freqS = 10 .^ d.freqLog; % spatial frequency (cycles/deg)
			d.freqT(:) = 4.22; % temporal frequency (Hz)
			d.amp = 10 .^ d.respLog; % response (Hz)
			d.phase(:) = 0; d.resp = d.amp; % store
			d.Properties.VariableDescriptions{'amp'} = 'Response amplitude (Hz)';
		case 'Wool' % Wool (18)
			f = [m.folder, filesep, 'Wool (18)/']; % data folder
			load([f, 'R1.mat'], 'R1'); % log amplitude (Hz): fs x 2
			load([f, 'Rphase.mat'], 'Rphase'); % phase (deg): fs x 2
			fs = size(R1, 1); % number of frequencies
			source = repmat("Wool", [fs, 1]); source = categorical(source); % source
			d = table(source); % make output table
			R1 = sortrows(R1); Rphase = sortrows(Rphase); % sort on frequency: fs x 2
			f = mean([R1(:, 1), Rphase(:, 1)], 2); % log frequency (cyc/deg): fs x 1
			d.freqS = 10 .^ f; % spatial frequency (cycles/deg): fs x 1
			d.freqT(:) = 2; % temporal frequency (Hz)
			a = 10 .^ R1(:, 2); % response amplitude (Hz): fs x 1
			p = Rphase(:, 2); % response phase (deg): fs x 1
			r = a .* exp(1i * pi * p / 180); % response (Hz), complex: fs x 1
			d.amp = a; d.phase = p; d.resp = r; % store
			d.Properties.VariableDescriptions{'amp'} = 'Response amplitude (Hz)';
			d.Properties.VariableDescriptions{'phase'} = 'Response phase (deg)';
		case 'Yeh' % Yeh (95)

			% Sensitivity data
			folder = [m.folder, filesep, 'Yeh (95) sens'];
			d = m.readFile(folder); % contrast sensitivity versus temporal frequency
			d = renamevars(d, 'label', 'source'); % rename for consistency
			d.source = categorical(d.source); % for listing
			d.freqS(:) = 0; % spatial frequency (cycles/deg)
			d.freqT = (.6 * 2 .^ (0: 6))'; % cleaned temporal frequency (Hz)
			a = 10 .^ d.sensLog; d.amp = a; % contrast sensitivity (Hz/%)
		
			% Phase data and response
			folder = [m.folder, filesep, 'Yeh (95) phase'];
			dC = m.readFile(folder); % response phase versus temporal frequency
			p = dC.phase; d.phase = p; % response phase (deg)
			d.resp = a .* exp(1i * pi * p / 180); % response (Hz), complex: fs x 1

			% Describe
			d.Properties.VariableDescriptions{'amp'} = 'Contrast sensitivity (Hz/%)';
			d.Properties.VariableDescriptions{'phase'} = 'Response phase (deg)';

	end
	d.Properties.VariableDescriptions{'freqS'} = 'Spatial frequency (cycles/deg)';
	d.Properties.VariableDescriptions{'freqT'} = 'Temporal frequency (Hz)';

function [d, m] = respCone(~, m) % cone frequency response: Baudin (19)

	folder = [m.folder, filesep, 'Baudin (19)'];
	d = m.readFile(folder); % contrast sensitivity versus temporal frequency
	d = renamevars(d, 'label', 'source'); % cone type and illuminance
	d.source = categorical(d.source); % improve listing
	f = [1, 2, 4, 8, 16, 20, 24, 32, 40, 48]'; % frequency (Hz)
	d.freq = repmat(f, [4, 1]); % clean frequency (Hz)

function [d, m] = respPulse(~, m) % pulse response: Lee (94), Sinha (17)

	% Switch source
	switch m.resp.source
		case 'Lee'

			% Read data
			folder = [m.folder, filesep, 'Lee (94)'];
			d = m.readFile(folder); % impulse rate versus time
			d = renamevars(d, 'label', 'dur'); % pulse duration (ms)
			d.dur = .001 * d.dur; % pulse duration (s)
			d = sortrows(d, {'dur', 'time'}); % sort rows
		
			%	Convert sample indices to times
			tSamp = .002; % duration of time sample (s)
			[dur, ~, g] = unique(d.dur); % group by duration
			durs = length(dur); % number of durations
			for i = 1: durs % loop over durations
				j = g == i; % indices for current group
				js = sum(j); % group size
				t = tSamp * (.5 + (0: js - 1)); % sample mid-point (ms)
				d.time(j) = t'; % store
			end
		
			%	Store
			d.Properties.VariableDescriptions{'dur'} = 'Pulse duration (s)';	
			d.Properties.VariableDescriptions{'time'} = 'Time (s)';

		case 'Sinha'

			% Read data
			folder = [m.folder, filesep, 'Sinha (17)'];
			d = m.readFile(folder); % impulse rate versus time
			d = renamevars(d, 'label', 'type'); % cone type
			d.type = categorical(d.type); % for readable listing
			d.Properties.VariableDescriptions{'type'} = 'Cone type';

	end

function [d, m] = respS(d, m) % inverse transform ROG spatial freq. response

	%	Initialise
	xs = 128; % number of x values
	b = m.rogPar(m); % calculate ROG parameters from full model parameters
	freqT = d.freqT; % temporal frequency (Hz)
	
	%	Calculate x values
	if isfield(m.respS, 'range') % range of x values (deg)
		range = m.respS.range; % defined by user
	else % default
		radBack = b(4); % surround radius (deg)
		range = 4 * radBack; % range of x values (deg)
	end
	x = .5 * range * linspace(-1, 1, xs + 1); % x values (deg): 1 x (xs + 1)
	x = x(1: xs); % open-ended interval (deg): 1 x xs

	%	Calculate spatial frequencies
	fund = 1 / range; % fundamental spatial frequency (cycles/deg)
	fMax = .5 * xs * fund; % maximum frequency (cycles/deg)
	f = fMax * linspace(-1, 1, xs + 1); % freq., increasing: 1 x (xs + 1)
	f = f(1: xs); % open-ended interval (deg): 1 x xs
	freqS = fftshift(f); % frequencies in fft order (cycles/deg): 1 x xs

	%	Calculate ROG and its transform
	fT = repmat(freqT, xs, 1); % one temp. freq. for each spat. freq.
	f = [freqS', fT]; % frequency (cycles/deg, Hz): xs x 2
	rF = m.calRog(b, f); % ROG response (Hz/contrast-unit): xs x 1
	rF = fund * rF; % normalise for range (Hz/(contrast-unit*deg)): xs x 1
	rF = xs * rF; % convert to fft units (Hz/(contrast-unit*deg)): xs x 1
	rS = ifft(rF); % inverse transform (Hz/(contrast-unit*deg)): xs x 1
	rS = ifftshift(rS); % responses in spatial order (Hz/(cont.-unit*deg)): xs x 1

	%	Store
	d.x = x; % x values (deg): 1 x xs
	d.freqS = freqS; % spatial frequency (cycles/deg): 1 x xs
	rF = shiftdim(rF, -1); rS = shiftdim(rS, -1); % prepare for storage: 1 x xs
	d.respF = rF; % ROG response in frequency domain (Hz/contrast-unit): 1 x xs
	d.respS = rS; % ROG response in spatial domain (Hz/contrast-unit): 1 x xs
	d.respSReal = real(rS); % real part of ROG response in spatial domain
	d.respSAmp = abs(rS); % magnitude
	rS = rad2deg(unwrap(angle(rS))); % phase (deg)
	if min(rS) < -180, rS = rS + 360; end % prevent very negative phases
	d.respSPhase = rS; % store
	d.Properties.VariableDescriptions = {'Temporal frequency (Hz)', ...
		'Location (deg)', 'Spatial frequency (cycles/deg)', ...
		'Frequency-domain response (Hz/cont.-unit)', ...
		'Spatial-domain response (Hz.contrast-unit^-1.deg^-1)', ...
		'Real response (Hz.cont-unit^-^1.deg^-^1)', ...
		'Response amplitude (Hz.cont-unit^-^1.deg^-^1)', ...
		'Response phase (deg)'};
