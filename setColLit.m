function m = setColLit(m) % set Colour metadata calculated from the literature

	% Files:
	%		'Cone density.mat'
	%		'Cortical magnification.mat'
	%		'Functional eccentricity.mat'
	%		'Ganglion cell density.mat'

	% Parameters calculated from the literature, e = functional eccentricity (deg)
	m.p.densGangCoef = [1.3497, -7.8485, 7.7392, 5.6107]; % g.c. density
		%	Wässle (89)
		% ganglion cell density (deg^-2) = exp(polyval(c, log10(eAnat)))
	m.p.kGangDev = .1241; % s.d. of gang. cell nearest-neighbour dist., Dacey (93)
	m.p.kGangGen = .47; % attenuation from g.c. to l.g.n., Kaplan (87)
	m.p.kRect = 7.2; % rectification constant (Hz/mV), Carandini (00)
	m.p.kSens = 17.8; % gang. cell contrast sens'y (mV/cont.-unit), Croner (95)
	m.p.kSur = 1.54; % surround gain, Croner (95)
	m.p.magRet = 4.73; % retinal magnification factor (deg/mm), Perry (85)
	m.p.potRest = 5.7; % g.c. resting pot. (mV), Kaplan (87): 41.1 / m.p.kRect
	m.p.radGangCoef = [-1.3583e-06, .00014598, .0012385, .012601]; % g.c. radius
		%	Watanabe (89)
		% ganglion cell dendritic radius (deg) = polyval(c, e)
	m.p.radHorCoef = [-4.5667e-05, .0030709, .029977, -3.261]; % h.c. r.f. radius
		% Wässle (89), Packer (02)
		% horizontal cell receptive field radius (deg) = exp(polyval(c, e))
	m.p.radOptCoef = [1.4241e-05, 0, .011663]; % point spread function radius
		% Navarro (93), Harwerth (85)
		% point spread function radius (deg) = polyval(c, e)
	m.p.ratCone = [0.4439, 0.4310, 0.1251]; % cone ratio: [L, M, S] / (L + M + S]
		%	Munds (22)
	m.p.ratGang = [0.485, 0.35]; % ratio of off-, on-midget to all ganglion cells
		%	Peng (19)
	m.p.ratOpt = 1.5; % ratio of PSF radius, macaque/human, Harwerth (85)
	m.p.tau = .014; % time constant (s), Lee (94)
