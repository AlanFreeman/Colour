function setCol

% Use values in the literature to set parameters for Colour model

	% Initialise
	m = setColLit; % set model parameters
	m = setColFun(m); % set handles for functions other than tasks
	m = setColTask(m); % set handles for task functions

	% Determine and display parameters
	switch 'resp.x fit' % select analysis tasks
		case 'count.ecc cone' % cone count vs eccentricity: Packer (89)
			m.tasks = 'densCone countCone plot set add';
			m.x = 'eccDeg'; m.y = 'count';
			m.set.axes = {'xLim', [0, 60]};
			m.set.axes = {'xScale', 'log', 'xLim', [.1, 60], ...
				'xTick', [.1, 1, 10, 60], 'yScale', 'log'};
			m.add.funFun = @(d, m) plot(1, interp1(d.eccDeg, d.count, 1), 'ok');
		case 'count.ecc cen' % % number of cones per ganglion cell centre
			m.tasks = 'densCone select countCen plot set';
			m.select.quad = 'temporal';
			m.x = 'eccDeg'; m.y = 'count';
			m.set.axes = {'xScale', 'log', 'xLim', [.01, 40], 'xTick', ...
				[.01, .1, 1, 10], 'yLim', [0, 20], 'yTick', [0, 10, 20], ...
				'clipping', 'on'};
		case 'count.ecc dend' % number of cones per ganglion cell dendritic field
			m.tasks = 'densCone select count prep plot set';
			m.select.quad = 'temporal';
			m.prep.zero = .1;
			m.x = 'eccDeg'; m.y = 'count';
			m.set.axes = ...
				{'xScale', 'log', 'xLim', [.1, 60], 'xTick', [.1, 1, 10, 60], ...
				'yScale', 'log', 'yLim', [1, 30], 'yTick', [1, 10, 30]};
		case 'count.ecc dens' % cones per ganglion cell via functional g.c. density
			m.tasks = 'densGang countGang densGangFun densGangSub ratDens prep';
			m.tasks = [m.tasks, ' ', 'plot set add'];
			m.x = 'ecc'; m.y = 'ratio';
			m.ratDens.group = 'type';
			m.prep.zero = .3;
			m.plot.line = 'type';
			m.set.axes = {'xScale', 'log', 'xLim', [.3, 60], ...
				'xTick', [.3, 1, 10, 60], 'yScale', 'log', 'yLim', [.3, 20], ...
				'yTick', [.5, 1, 10]};
			m.add.funFun = @(d, m) plot(.5, 4484.03 / 9275.5, 'ok');
				% from ganglion cell count
		case 'count.ecc gang' % ganglion cell count vs functional ecc.: Wässle (89)
			m.tasks = 'densGang countGang plot set add';
			m.y = 'count';
			switch 'func' % type of eccentricity
				case 'anat' % anatomical
					m.tasks = 'densGang countGang plot set add';
					m.x = 'eccAnat';
					m.set.axes = {'xScale', 'log', 'xLim', [1, 60], ...
						'yScale', 'log', 'yLim', [100, 1e7], 'clipping', 'off'};
					e = m.p.magRet * .55 / .9; n = 33000; % data for 550 um circle
					m.add.funFun = @ (d, m) plot(e, n, 'o');
				case 'func' % functional
					m.x = 'ecc';
					switch 'log' % linear or log?
						case 'lin' % linear
							lim = 1;
							m.set.axes = {'xLim', lim * [0, 1], 'xTick', lim * [0, .5, 1], ...
								'yLim', [0, 4e4], 'yTick', 4e4 * [0, .5, 1]};
						case 'log' % logarithmic
							m.set.axes = {'xScale', 'log', 'xLim', [.01, 100], ...
								'yScale', 'log', 'yLim', [100, 1e7], 'clipping', 'off'};
					end
			end
			m.add.funFun = @(d, m) plot(1, interp1(d.ecc, d.count, 1), 'ok');
		case 'cover.ecc' % ganglion cell coverage
			m.tasks = ['densGang countGang densGangFun densGangSub ', ...
				'select cover plot set'];
			m.x = 'ecc'; m.y = 'cover';
			m.select.type = 'mid';
			m.set.axes = {'xLim', [0, 60], 'xTick', [0, 30, 60], ...
				'yLim', [0, 5], 'yTick', [0, 2.5, 5]};
			m.set.axes = {'xScale', 'log', 'xLim', [.01, 40], ...
				'xTick', [.01, .1, 1, 10, 40], 'yLim', [0, 15], 'yTick', [0, 7.5, 15]};
		case 'dens.ecc cone' % cone density: Packer (89)
			m.tasks = 'densCone select prep plot set add';
			%	m.tasks = 'densCone save';
			m.select.quad = 'temporal';
			m.plot.group = 'quad';
			switch 'both' % data to plot: both, far or near
				case 'both' % both far and near
					m.x = 'eccDeg'; m.y = 'densDeg'; m.plot.arg = {'-'};
					yLim = [3, 2e4];
					switch 'log' % linear or logarithmic x axes
						case 'lin' % linear
							xLim = 55; xTick = 50;
							m.set.axes = {'xLim', xLim * [0, 1], 'xTick', ...
								xTick * [0, .5, 1], 'yScale', 'log', 'yLim', yLim};
						case 'log' % log
							z = .3; m.prep.zero = z; % foveal proxy
							m.set.axes = {'xScale', 'log', 'xLim', [z, 60], ...
								'xTick', [z, 1, 10, 60], 'yScale', 'log', 'yLim', yLim, ...
								'clipping', 'on'};
					end
				case 'far' % as published
					m.select.near = 0;
					m.x = 'ecc'; m.y = 'dens'; m.plot.arg = {'-s'}; lim = 1000;
					m.set.axes = {'xLim', [0, 20], 'xTick', [0, 10, 20], ...
						'yLim', lim * [1, 25], 'yTick', lim * [5, 15, 25]};
				case 'near' % as published
					m.select.near = 1;
					m.x = 'ecc'; m.y = 'dens'; m.plot.arg = {'-s'}; lim = 1000;
					m.set.axes = {'xLim', [0, 1.25], 'xTick', [.25, .75, 1.25], ...
						'yLim', lim * [25, 225], 'yTick', lim * [25, 125, 225]};
			end
			m.add.funFun = @(d, m) plot(.5, 14087 / pi, 'ok'); % from cone count
			m.save.name = 'Cone density.mat';
		case 'dens.ecc gang' % ganglion cell density: Wässle (89)
			m.tasks = 'densGang select plot set';
			m.select.type = {'anatNarr', 'anatWide'};
			switch 'deg' % which units?
				case 'deg' % convert to deg
					m.x = 'eccDeg'; m.y = 'densDeg';
					yLim = [3, 2e4];
					switch 'log' % linear or logarithmic eccentricity?
						case 'lin' % linear
							m.set.axes = {'xLim', [0, 55], 'xTick', [0, 25, 50], ...
								'yScale', 'log', 'yLim', yLim};
						case 'log' % logarithmic
							m.set.axes = {'xLim', [.3, 60], 'xTick', [.3, 1, 10, 60], ...
								'xScale', 'log', 'yScale', 'log', 'yLim', yLim, ...
								'clipping', 'off'};
					end
				case 'mm' % plot as published
					m.x = 'ecc'; m.y = 'dens';
					xLim = 12; xTick = [0, 6, 12]; yLim = [100, 1e6];
					m.set.axes = {'xLim', xLim * [0, 1], 'xTick', xTick, ...
						'yScale', 'log', 'yLim', yLim};
			end
			m.plot.line = 'type'; m.plot.arg = {'o', 'clipping', 'off'};
		case 'dens.ecc gang fit' % fit ganglion cell density: Wässle (89)
			switch 'plot'
				case 'check'
					m.tasks = 'densGang select plot set add';
					m.add.funFun = @(d, m)plot(d.eccDegLog, ...
						exp(polyval(m.p.densGangCoef, d.eccDegLog)));
				case 'fit'
					m.tasks = 'densGang select plot set fitglm pred add';
				case 'plot'
					m.tasks = 'densGang select plot set';
				case 'show'
					m.tasks = 'densGang select fitglm show'; % set m.p.densGangCoef here
			end
			m.select.type = {'anatNarr', 'anatWide'};
			m.x = 'eccDegLog'; m.y = 'densDeg';
			m.plot.line = 'type';
			m.plot.arg = {'o', 'clipping', 'off'};
			m.set.axes = {'xLim', log10([.3, 60]), 'yScale', 'log', ...
				'yLim', [3, 2e4], 'clipping', 'off'};
			m.fitglm.arg = {'poly3', 'predictorVars', m.x, 'responseVar', m.y, ...
						'link', 'log'};
			m.add.line = {};
		case 'dens.ecc gang fitlm' % fit ganglion cell density: Wässle (89)
			m.tasks = 'densGang select plot set fitlm show pred add';
			m.select.type = {'anatNarr', 'anatWide'};
			m.x = 'eccDegLog'; m.y = 'densDegLog';
			%	m.plot.line = 'type';
			m.plot.arg = {'o', 'clipping', 'off'};
			m.set.axes = {'xScale', 'log', 'xLim', [.3, 60], ...
				'xTick', [.3, 1, 10, 60], 'yScale', 'log', 'yLim', [3, 2e4], ...
				'clipping', 'off'}; % *** fix ***
			m.fitlm.arg = {'poly3', 'predictorVars', m.x, 'responseVar', m.y};			
		case 'dens.ecc gang fun' % g.c. density with functional ecc.: Wässle (89)
			m.tasks = 'densGang countGang densGangFun prep plot set';
			m.x = 'ecc'; m.y = 'densDeg';
			switch 'log' % x axis scaling
				case 'lin' % linear
					xLim = 55; yLim = [3, 3e4];
					m.set.axes = {'xLim', xLim * [0, 1], 'xTick', [0, 25, 50], ...
						'yScale', 'log', 'yLim', yLim, 'clipping', 'off'};
				case 'log' % logarithmic
					m.prep.zero = .1;
					m.set.axes = {'xScale', 'log', 'xLim', [.1, 60], ...
						'xTick', [.1, 1, 10, 60], 'yScale', 'log', 'yLim', [4, 2e4], ...
						'clipping', 'off'};
			end
		case 'dens.ecc gang sub' % g.c. dens. subpopulations: multiple authors
			m.tasks = 'densGang countGang densGangFun densGangSub prep plot set';
			m.tasks = 'densGang countGang densGangFun densGangSub prep plot set add';
			m.tasks = 'densGang countGang densGangFun densGangSub select save';
			z = .3; m.prep.zero = z;
			m.select.type = 'mid';
			m.x = 'ecc'; m.y = 'densDeg';
			m.plot.line = 'type';
			m.set.axes = {'xScale', 'log', 'xLim', [z, 60], 'xTick', ...
				[z, 1, 10, 60], 'yScale', 'log', 'yLim', [3, 2e4]};
			m.add.funFun = @(d, m) plot(.5, (34898 / pi) * [1, sum(m.p.ratGang)], ...
				'o'); % from ganglion cell count
			m.save.name = 'Ganglion cell density';
		case 'eccFun.ecc' % func. ecc. vs anat. ecc.: McGregor (18), Schein (88)
			m.tasks = 'eccFun plot set add'; % plot
			m.x = 'ecc'; m.y = 'eccFun';
			m.plot.line = 'source'; m.plot.arg = {'o'};
			lim = 15;
			m.set.axes = {'xLim', lim * [0, 1], 'xTick', lim * [0, .5, 1], ...
				'yLim', lim * [0, 1], 'yTick', lim * [0, .5, 1], 'clipping', 'on'};
			m.add.funFun = @(d, m)plot([0, 90], [0, 90], '--');
		case 'eccFun.ecc fit' % fun. vs anat. ecc. + fit: McGregor (18), Schein (88)
			m.tasks = 'eccFun plot set eccFunFit add'; % plot
			%	m.tasks = 'eccFun eccFunFit save'; % save
			%	m.eccFunFit.mcGregor = .03; % lower limit on McGregor fun. ecc.
			m.eccFunFit.schein = 5; % lower limit on Schein eccentricities
			m.x = 'ecc'; m.y = 'eccFun';
			m.plot.line = 'source'; m.plot.arg = {'o'};
			switch 'log' % axis scaling
				case 'lin' % linear
					lim = 10;
					m.set.axes = {'xLim', lim * [0, 1], 'xTick', lim * [0, .5, 1], ...
						'yLim', lim * [0, 1], 'yTick', lim * [0, .5, 1], 'clipping', 'on'};
				case 'log' % logarithmic
					m.set.axes = ...
						{'xScale', 'log', 'xLim', [1, 10], 'xTick', [1, 3, 10], ...
						'yScale', 'log', 'yLim', [.007, 10], 'yTick', [.01, .1, 1, 10], ...
						'clipping', 'on'};
			end
			m.add.funFun = @(d, m) plot(d.ecc, d.eccFun, 'k');
			m.save.name = 'Functional eccentricity';
		case 'eccFun.ecc mcG' % func. ecc. vs anat. ecc.: McGregor (18)
			m.tasks = 'eccFunMcG plot set';
			m.x = 'eccFunUm'; m.y = 'eccUm'; m.plot.arg = {'o'};
			m.set.axes = {'xLim', [0, 140], 'xTick', [0, 70, 140], ...
				'yLim', [100, 600], 'yTick', [100, 300, 500]};
		case 'gain' % surround gain: Croner (95)
			freq = 2 * pi * 4.22; % stimulus temporal frequency (radians/s)
			a = 1 + 1i * m.p.tau * freq; % temporal attenuation
			r = .547; % surround response / centre response at 0 spatial frequency
			switch 'fzero' % method
				case 'direct' % algebraic solution
					k = (a ^ 2) * ((1 / r) - 1);
					fprintf('Surround gain magnitude, phase (deg): %3g, %3g\n', ...
						abs(k), (180 / pi) * angle(k)); % list
				case 'fzero' % numerical solution
					fun = @ (k) abs(1 + k / a ^ 2) - 1 / (1 - r); % find zero for this fun
					m.p.kSur = fzero(fun, 1); % solve for kSur
					fprintf('Surround gain: %3g\n', m.p.kSur); % list
			end
			return
		case 'mtf.freq' % modulation transfer function: Navarro (93)
			m.tasks = 'mtf plot set';
			m.x = 'freq'; m.y = 'mtf';
			m.plot.line = 'ecc';
			m.set.axes = {'yLim', [.005, 1], 'yScale', 'log'};			
		case 'offset.ecc' % ganglion cell offset vs. eccentricity: Schein (88)
			m.tasks = 'offset plot set';
			m.plot.line = 'source'; m.plot.arg = {'clipping', 'off'};
			m.x = 'eccMm'; m.y = 'offsetMm'; % m.plot.arg = {'o'};
			m.set.axes = {'xLim', [0, 2.5], 'xTick', 0: 2, ...
				'yLim', [0, .4], 'yTick', [0, .2, .4]};
		case 'prof.x' % plot horizontal cell r.f. or conv. function profiles
			m.tasks = 'profHor plot set export';
			m.x = 'x'; m.y = 'conv';
			lim = .4;
			switch m.y % conv, field
				case 'conv' % plot h.c.-cone convergence function profiles
					m.radHor.iter = 10;
					m.radHor.width = 1.6; lim = .5 * m.radHor.width;
				case 'field' % plot horizontal cell receptive field profiles
					m.radHor.seed = 2;
			end
			m.plot.funFun = @(d, m)plot(d.(m.x), shiftdim(d.(m.y), 1));
			m.set.axes = {'xLim', lim * [-1, 1], 'xTick', lim * [-1, 0, 1], ...
				'yLim', [0, 1], 'yTick', [0, .5, 1]};
		case 'psf.loc' % point spread function: Navarro (93)
			m.tasks = 'mtf psf plot set';
			switch 'psf' % function to display
				case 'mtf' % modulation transfer function
					m.x = 'freq'; m.y = 'mtf';
					m.plot.line = 'ecc';
					m.set.axes = {'xLim', [0, 60], 'yScale', 'log', 'yLim', [.005, 1]};
				case 'psf' % point spread function
					m.x = 'loc'; m.y = 'psf';
					m.plot.line = 'ecc';
			end
		case 'rad beta' % 4CBeta radius: multiple authors
			m.tasks = 'radBeta plot set'; % individual radii
			m.tasks = 'radBeta plot set list'; % m.p.radBetaMm = convergence radius
			m.x = 'loc'; m.y = 'radius'; m.plot.arg = {'o'};
			m.plot.line = 'source';
			m.set.line = 'fill';
		case 'rad.ecc' % centre, surround radius decomposition: multiple authors
			switch 'sur' % choose mechanism
				case 'cen' % centre mechanism
					m.tasks = 'rad select plot set'; % mechanisms
					m.select.source = {'optics', 'denGang', 'cenConv'};
					%	m.tasks = 'rad select plot set radCen add export'; % empirical data
					%	m.select.source = {'cenConv'};
					yLim = [.007, .3]; yTick = [.01, .1];
				case 'sur' % surround mechanism
					m.tasks = 'rad select plot set export'; % mechanisms
					m.rad.horExp = 1; % horizontal cells have exponential receptive fields
					m.select.source = {'optics', 'fieldHor', 'surConv'};
					%	m.tasks = 'rad select plot set radSur add export'; % empirical data
					%	m.select.source = {'surConv'};
					yLim = [.007, 5]; yTick = [.01, .1, 1];
			end
			m.x = 'ecc'; m.y = 'radius';
			m.plot.line = 'source';
			m.add.arg = {'o'};
			m.set.axes = {'xScale', 'log', 'xLim', [.006, 40], ...
				'xTick', [.01, .1, 1, 10], 'yScale', 'log', 'yLim', yLim, ...
				'yTick', yTick, 'clipping', 'off'};
		case 'rad.ecc beta' % geniculocortical convergence radius: multiple authors
			m.tasks = 'radBetaEcc plot set';
			m.p.file = 'Cortical magnification';
			m.x = 'ecc'; m.y = 'radius';
			m.set.axes = {'xLim', [0, 50], 'xTick', [0, 25, 50]};
		case 'rad.ecc cen' % gang. cell centre radius: Croner (95), Lee (98)
			m.tasks = 'radCen select plot set'; % no fitting
			m.select.source = 'Lee';
			m.x = 'ecc'; m.y = 'radius';
			m.plot.line = 'source';
			m.plot.arg = {'o'};
			switch m.select.source % plot in published form
				case 'Croner' % Croner (95)
					m.set.axes = {'xLim', [0, 40], 'xTick', [0, 20, 40], 'yLim', [0, .3]};
				case 'Lee' % Lee (98)
					m.y = 'dev';
					m.set.axes = {'xLim', [0, 15], 'xTick', [0, 7.5, 15], ...
						'yScale', 'log', 'yLim', [.1, 10]};
			end
		case 'rad.ecc cen fit' % fit gang. cell centre radius vs ecc.: Croner (95)
			m.tasks = 'radCen select plot fitlm pred add set'; % show prediction
			m.tasks = 'radCen select fitlm show'; % show model
			m.tasks = 'radCen select plot add set'; % check stored parameters
			m.select.source = 'Croner';
			m.x = 'ecc'; m.y = 'radiusLog';
			m.plot.line = 'source';
			m.plot.arg = {'o'};
			m.fitlm.arg = {'poly3', 'predictorVars', 'eccFun', ...
				'responseVar', 'radiusLog'};
			if ~ contains(m.tasks, 'fitlm') % check stored parameters
				p = m.p.radCen; m.add.funFun = @(d, m)plot(d.ecc, p(1) + ...
					p(2) * d.ecc + p(3) * d.ecc .^ 2 + p(4) * d.ecc .^ 3);
			end
			m.set.axes = {'xLim', [0, 40], 'xTick', [0, 20, 40]};
		case 'rad.ecc cone' % cone inner segment diameter: Packer (89)
			m.tasks = 'radCone';
		case 'rad.ecc gang' % ganglion cell dendritic diameter: Watanabe (89)
			m.tasks = 'radGang plot set'; % plot as published
			m.x = 'eccMm'; m.y = 'diam';
			m.plot.arg = {'o', 'clipping', 'off'};
			m.set.axes = {'xLim', [0, 14], 'xTick', [0, 7, 14], ...
				'yScale', 'log', 'yLim', [1, 500]};		
		case 'rad.ecc gang field' % gang. cell cen., sur. radius: Godat (22)
			m.tasks = 'radGangField plot set';
			m.x = 'ecc'; m.y = 'radius';
			m.plot.arg = {'o'}; m.plot.line = 'type';
			m.set.axes = {'xScale', 'log', 'xLim', [.006, .1], 'xTick', [.01, .1], ...
				'yScale', 'log', 'yLim', [.01, .1]};
		case 'rad.ecc gang fit' % ganglion cell dendritic diameter: fitted curve
			switch 'show' % select analysis tasks
				case 'check' % check stored parameters
					m.tasks = 'radGang plot set add'; % check stored parameters
					m.add.funFun = @(d, m)plot(d.eccFun, ...
						polyval(m.p.radGangCoef, d.eccFun));
				case 'fit' % fit regression model
					m.tasks = 'radGang plot set fitlm pred add'; % add fitted curve
				case 'show' % show regression model
					m.tasks = 'radGang fitlm show'; % model: set m.p.radGangCoef from this
			end
			m.x = 'eccFun'; m.y = 'radius'; m.plot.arg = {'o'};
			m.fitlm.arg = {'poly3', 'predictorVars', m.x, 'responseVar', m.y};
			m.set.axes = {'xLim', [0, 60], 'xTick', [0, 30, 60], ...
				'yScale', 'log', 'yLim', [.007, .4]};
			m.set.axes = {'xScale', 'log', 'xLim', [.1, 60], ...
				'xTick', [.1, 1, 10, 60], 'yScale', 'log', ...
				'yLim', [.007, .4], 'clipping', 'off'};
		case 'rad.ecc hor dend' % hor. cell dendrite area: Wässle (89) Horizontal
			m.tasks = 'radHorDend plot set'; % plot as published
			m.x = 'ecc'; m.y = 'area';
			m.plot.group = 'type'; m.plot.arg = {'o', 'clipping', 'off'};
			m.set.axes = {'xLim', [0, 14], 'xTick', [0, 7, 14], ...
				'yScale', 'log', 'yLim', [1e-4, 2e-2]};			
		case 'rad.ecc hor dend fit' % horizontal cell dendrite area: fitted curve
			m.tasks = 'radHorDend plot fitlm pred add set'; % add fitted curve
			%	m.tasks = 'radHorDend fitlm show'; % show regression model
			%	m.tasks = 'radHorDend plot add set'; % check stored parameters
			m.x = 'eccFun'; m.y = 'areaLog';
			m.plot.group = 'type'; m.plot.arg = {'o', 'clipping', 'off'};
			m.fitlm.group = 'type'; m.pred.group = 'type';
			m.fitlm.arg = {'poly3', 'predictorVars', m.x, 'responseVar', m.y};
			if ~ contains(m.tasks, 'fitlm') % check stored parameters
				p = m.p.area; m.add.funFun = @(d, m)plot(d.eccFun, p(1) + ...
					p(2) * d.eccFun + p(3) * d.eccFun .^ 2 + p(4) * d.eccFun .^ 3);
			end
			m.set.axes = {'xLim', [0, 70], 'xTick', [0, 35, 70], ...
				'yLim', [-4, log10(.02)]};			
		case 'rad.ecc hor field' % horizontal cell receptive field: Packer (02)
			m.tasks = 'radHorField select plot set'; % no fitting
			m.select.source = {'narrow', 'wide'};
			m.x = 'ecc'; m.y = 'diam';
			m.plot.arg = {'o'}; m.plot.line = 'source';
			m.set.axes = {'xLim', [0, 16], 'xTick', [0, 8, 16], ...
				'yLim', [0, 1200], 'yTick', [0, 400, 800]};
		case 'rad.ecc hor field fit' % hor. cell radius: Wässle (89), Packer (02)
			m.tasks = 'radHor select plot set radHorFit add';
			m.x = 'eccFun'; m.y = 'radius';
			m.radHor.iter = 10;
			m.radHor.method = 'conv';
			m.select.source = {'dend', 'field'};
			ana = 'check';
			switch ana % select analysis tasks
				case 'check' % check stored parameters
					m.tasks = 'radHor select plot set add';
					m.radHor.double = 1;
					m.add.funFun = @(d, m)plot(d.eccFun, ...
						10 .^ (m.p.radBackCoef(1) + m.p.radBackCoef(2) * ...
						normcdf(d.eccFun, m.p.radBackCoef(3), m.p.radBackCoef(4))));
				case 'fit' % fit regression model
					m.radHor.double = 1;
				case 'plot' % no fit
					m.tasks = 'radHor select plot set radHorDiam add export';
					m.select.source = {'dend', 'wide'};
					m.add.funFun = @(d, m) plot(d.radLower(1, :), d.radLower(2, :), ...
						d.radUpper(1, :), d.radUpper(2, :), '-');
				case 'show' % show regression model
					m.tasks = 'radHor select radHorFit show'; % set m.p.radBackCoef here
			end
			m.plot.arg = {'o'}; m.plot.line = 'source';
			m.set.axes = {'xLim', [0, 60], 'xTick', [0, 30, 60], ...
				'yScale', 'log', 'yLim', [.035, 1.5]};
		case 'rad.ecc hor field fit Gauss' % h.c. radius: Wässle (89), Packer (02)
			m.x = 'eccFun'; m.y = 'radius';
			switch 'plot' % select analysis tasks
				case 'check' % check stored parameters
					m.tasks = 'radHor select plot set add';
					m.select.source = {'dend', 'conv'};
					m.add.funFun = @(d, m)plot(d.eccFun, ...
						exp(polyval(m.p.radHorCoef, d.eccFun)));
				case 'fit' % fit regression model
					m.tasks = 'radHor select plot set fitglm pred add';
					m.select.source = {'dend', 'field'};
				case 'plot' % no fit
					m.tasks = 'radHor select plot set';
					m.select.source = {'dend', 'wide'};
				case 'show' % show regression model
					m.tasks = 'radHor select fitglm show'; % set m.p.radHorCoef from list
					m.select.source = {'dend', 'conv'}; % *** out of date?
			end
			m.radHor.method = 'union'; % sum, union
			m.fitglm.arg = {'poly3', 'predictorVars', m.x, 'responseVar', m.y, ...
				'link', 'log'};
			m.plot.arg = {'o'}; m.plot.line = 'source';
			m.set.axes = {'xLim', [0, 60], 'xTick', [0, 30, 60], ...
				'yScale', 'log', 'yLim', [.035, 1.3], 'clipping', 'off'};
		case 'rad.ecc opt' % retinal resolution function: Navarro (93)
			m.tasks = 'radOpt plot set';
			switch 'show' % choose tasks to perform
				case 'add' % add fitted model to empirical data
					m.tasks = 'radOpt plot set fitlm pred add';
				case 'check' % check stored coefficients for fitted model
					m.tasks = 'radOpt plot set add';
					m.add.funFun = @(d, m)plot(d.ecc, polyval(m.p.radOptCoef, d.ecc));
				case 'save' % save empirical data to file
					m.tasks = 'radOpt save'; % save data
					m.save.name = 'Point spread function.mat'; % use polynomial instead
				case 'show' % show statistics for fitted model
					m.tasks = 'radOpt plot set fitlm show'; % set m.p.radOptCoef from list
			end
			m.x = 'ecc'; m.y = 'radius'; m.plot.arg = {'-o'};
			m.fitlm.arg = {'radius ~ ecc^2 - ecc'};
			switch 'psf' % data to display
				case 'rrf' % retinal resolution function, as published
					m.y = 'rrf';
					m.set.axes = {'xLim', [-60, 60], 'xTick', [-60, 0, 60], ...
						'yScale', 'log', 'yLim', [.7, 100], 'yTick', [1, 10, 100]};
				case 'psf' % radius of point spread function
					m.set.axes = {'xLim', [-60, 60], 'xTick', [-60, 0, 60], ...
						'yLim', [0, .075], 'yTick', [0, .025, .05]};
			end
		case 'rad.ecc sur' % ganglion cell surround radius: Croner (95), Lee (98)
			m.tasks = 'radSur select plot set'; % no fitting
			m.x = 'eccFun'; m.y = 'radius';
			m.plot.line = 'source'; m.plot.arg = {'o'};
			m.set.axes = {'xScale', 'log', 'xLim', [.1, 100], ...
				'yScale', 'log', 'yLim', [.01, 10], 'clipping', 'off'};
		case 'rad.ecc sur fit' % ganglion cell surround radius: fitted curve
			% *** untested code ***
			m.tasks = 'radSur select plot fitlm pred add set'; % show prediction
			m.tasks = 'radSur select fitlm show'; % show model
			m.tasks = 'radSur select plot add set'; % check stored parameters
			m.select.source = 'sur';
			m.x = 'ecc'; m.y = 'radiusLog';
			m.plot.line = 'source';
			m.plot.arg = {'o'};
			m.fitlm.arg = {'poly3', 'predictorVars', 'ecc', ...
				'responseVar', 'radiusLog'};
			if ~ contains(m.tasks, 'fitlm') % check stored parameters
				p = m.p.radSur; m.add.funFun = @(d, m)plot(d.ecc, p(1) + ...
					p(2) * d.ecc + p(3) * d.ecc .^ 2 + p(4) * d.ecc .^ 3);
			end
			m.set.axes = {'xLim', [0, 40], 'xTick', [0, 20, 40]};
		case 'rat cone' % cone ratio: munds (22)
			switch 'munds' % choose source
				case 'mult' % multiple authors
					l = .6; % L / (L + M): Dacey et al (00) "Physiology ...", temporal q.
					lm = [l, 1 - l]; % [L, M] / (L + M)
					s(1) = .09934; % S / (L + M): Martin et al. (99); ecc. = 4-67 deg temp
					s(1) = s(1) / (1 + s(1)); % S / (L + M + S): Martin et al. (99)
					s(2) = .073; % S / (L + M + S): Roorda et al. (01); ecc. = 1-1.5 deg
					s = mean(s); % average
					r = [(1 - s) * lm(1), (1 - s) * lm(2), s]; % L: M: S ratio
				case 'munds' % Munds (22)
					lm = 1.03; % L / M
					s = .143; % S / (L + M)
					r = [lm / (1 + lm), 1 / (1 + lm), s] / (1 + s);
						% [L, M, S] / (L + M + S), m.p.ratCone = r
			end
			fprintf('Cone ratio, [L, M, S]: [%4g, %4g, %4g]\n', r); % report
			return
		case 'rat gang' % ratio of midget to all ganglion cells: multiple authors
			switch 'peng' % source
				case 'mult' % multiple authors: off-midget / all midget ganglion cells
					r = [.63, .62, .53]; % Dacey (93), Peng (19), Rhoades (19)
					r = mean(r); % average over studies
					fprintf('Off-midget / all midget ganglion cells: %4g\n', r); % report
				case 'peng' % Peng (19)
					r = [.47, .37; .5, .33]; % ratio of midget to all ganglion cells
						% row 1: fovea; row 2: periphery;
						% column 1: off-centre; column 2: on-centre
					r = mean(r); % mean of fovea, periphery: m.p.ratGang = r
					fprintf('Ratio of off-, on-midget to all ganglion cells: %g, %g\n', r);
			end
			return
		case 'rat gang dev' % nearest neighbour gang. cells: s.d./dist., Dacey (93)
			dist = 145; % mean distance between nearest neighbour on-centre g.c. (um)
			dev = 18; % standard deviation of distance (um)
			r = dev / dist; % ratio, m.p.kGangDev = r
			fprintf('Nearest neighbour gang. cells distance: s.d. / dist. = %g\n', r);
			return
		case 'rat ret' % retinal magnification ratio: Perry, Cowey (85)
			% Also see De Monasterio et al. (85): 1 / .203 = 4.93
			f(1) = 1 / .223; % M. mulatta
			f(2) = 1 / .201; % M. fascicularis
			f = mean(f); % average over species (deg/mm): set m.p.magRet = f
			fprintf('Retinal magnification factor (deg / mm): %4g\n', f); % report
			return
		case 'rat.ecc cum' % ratio of cone count/dens. to ganglion cell count/dens.
			m.tasks = 'densGang countGang densGangSub ratCount plot set export';
			m.densGangSub.z = 'count';
			m.ratCount.group = 'type';
			m.x = 'ecc'; m.plot.line = 'type';
			switch 'count'
				case 'count'
					m.y = 'ratio'; m.plot.arg = {'-o'};
					m.set.axes = {'xScale', 'log', 'xLim', [.3, 60], ...
						'xTick', [.3, 1, 10, 60], 'yScale', 'log', 'yLim', [.3, 20], ...
						'yTick', [.5, 1, 10, 20]};
				case 'dens'
					m.y = 'densDeg'; m.plot.arg = {'o-'};
					m.set.axes = {'xScale', 'log', 'xLim', [.3, 60], ...
						'xTick', [.3, 1, 10, 60], 'yScale', 'log', 'yLim', [3, 2e4]};
			end
		case 'rat.ecc mid' % midget / all ganglion cells: Grünert (93), Dacey (94)
			% *** unused: delete ***
			m.tasks = 'ratMidget stop plot set';
			m.plot.arg = {'-o'};
			switch 'deg' % unit for eccentricity
				case 'deg' % degrees
					m.x = 'eccDeg'; m.y = 'ratio';
					m.set.axes = {'xLim', [0, 75], 'yScale', 'log', 'yLim', [.1, 1]};
				case 'mm' % (mostly) as published
					m.x = 'ecc'; m.y = 'ratioPer';
					m.set.axes = {'xLim', [0, 15], 'yScale', 'log', 'yLim', [1, 100]};
			end
		case 'resp.freq' % frequency response: Benardete, Croner, Wool, Yeh
			m.tasks = 'resp plot set'; % plot empirical data
			m.resp.source = 'Yeh'; % data source: Benardete, Croner, Wool, Yeh
			m.y = 'amp'; % polar component: amp or phase
			m.plot.group = 'source'; m.plot.arg = {'o-'};
			switch m.resp.source
				case 'Benardete' % Benardete (97)
					m.plot.group = 'freqT'; m.x = 'freqS';
					switch m.y
						case 'amp'
							m.set.axes = {'xLim', [-20, 20], 'yLim', [0, 150]};
						case 'phase'
							m.set.axes = {'xLim', [-20, 20], 'yLim', [-360, 360]};
					end
				case 'Croner' % Croner (95)
					m.x = 'freqS';
					m.set.axes = {'xScale', 'log', 'xLim', [.1, 100], ...
						'yScale', 'log' 'yLim', [1, 100]};
				case 'Wool' % Wool (18)
					m.x = 'freqS';
					switch m.y
						case 'amp'
							m.set.axes = {'xScale', 'log', 'xLim', [.01, 7], ...
								'yScale', 'log' 'yLim', [1, 100]};
						case 'phase'
							m.set.axes = {'xScale', 'log', 'xLim', [.01, 7], ...
								'yLim', [-180, -90], 'yTick', [-180, -135, -90]};
					end
				case 'Yeh' % Yeh (95)
					m.x = 'freqT';
					switch m.y
						case 'amp'
							m.set.axes = {'xScale', 'log', 'xLim', [.1, 100], ...
								'yScale', 'log', 'yLim', [.1, 10]};
						case 'phase'
							m.set.axes = {'xScale', 'log', 'xLim', [.1, 100], ...
								'yLim', [-720, 90], 'yTick', [-540, -360, -180, 0]};
					end
			end
		case 'resp.freq fit' % fitted freq. resp.: Benardete, Croner, Wool, Yeh
			m.tasks = 'resp plot set fitRog show prep add'; % no refinement
			m.tasks = 'resp select plot set fitRog show refine predRog prep add';
			m.y = 'amp'; % polar component: amp, phase
			m.resp.source = 'Benardete'; % data source: Benardete, Croner, Wool, Yeh
			m.fitRog.funFun = 'rog'; % function to calculate: dog or rog
			m.plot.group = 'source'; m.plot.arg = {'o'};
			switch m.resp.source % set plotting data
				case 'Benardete'
					m.resp.disp = 0.010715; % stimulus displacement correction: Benardete
					m.plot.group = 'freqT';
					m.prep.group = 'freqT'; m.prep.unwrap = 'phase';
					m.x = 'freqS'; m.refine.group = 'freqT';
					x0 = .04; % replacement for 0
					switch m.y % set plotting data
						case 'amp'
							m.set.axes = { ...
								'xScale', 'log', 'xLim', [x0, 20], 'xTick', [.1, 1, 10], ...
								'yScale', 'log', 'yLim', [10, 150], 'yTick', [10, 100]};
							%	m.set.axes = { 'xLim', [0, 20], 'xTick', [0, 10, 20], ...
							%	'yScale', 'log', 'yLim', [10, 150], 'yTick', [10, 100]};
							% plot x = 0
						case 'phase'
							m.set.axes = { ...
								'xScale', 'log', 'xLim', [x0, 20], 'xTick', [.1, 1, 10], ...
								'yLim', [-540, -120], 'yTick', [-540, -360, -180]};
							%	m.set.axes = {'xLim', [0, 20], 'xTick', [0, 10, 20], ...
							%	'yLim', [-540, -120], 'yTick', [-540, -360, -180]}; % plot x = 0
							m.add.funFun = @(d, m) plot(d.(m.x), d.(m.y) - 360);
					end
				case 'Croner'
					m.x = 'freqS';
					m.set.axes = {'xScale', 'log', 'xLim', [.1, 100], ...
						'yScale', 'log' 'yLim', [1, 100]};
				case 'Wool'
					m.x = 'freqS';
					switch m.y % set plotting data
						case 'amp'
							m.set.axes = {'xScale', 'log', 'xLim', [.01, 7], ...
								'yScale', 'log' 'yLim', [1, 100]};
						case 'phase'
							m.set.axes = {'xScale', 'log', 'xLim', [.01, 7], ...
								'yLim', [-180, 0], 'yTick', [-180, -90, 0]};
					end
				case 'Yeh'
					m.x = 'freqT';
					switch m.y
						case 'amp'
							m.set.axes = {'xScale', 'log', 'xLim', [.5, 40], ...
								'yScale', 'log', 'yLim', [.5, 3], 'yTick', [.5, 1, 2]};
						case 'phase'
							m.set.axes = {'xScale', 'log', 'xLim', [.5, 40], ...
								'yLim', [-450, 90], 'yTick', [-360, -180, 0]};
					end
			end
			switch m.fitRog.funFun % set optimisation data
				case 'dog' % difference of Gaussians
					m.fitRog.group = 'freqT'; m.predRog.group = 'freqT';
					m.fitRog.name = {'kCen', 'kSur', 'rCen', 'rSur', 'pCen', 'pSur', ...
						'disp'}; % coefficient names
					switch m.resp.source
						case 'Benardete'
							m.select.freqT = 16.92; % 2.17, 8.46, 16.92
							m.fitRog.fix = {'disp'}; % fixed coefficients
							switch m.select.freqT % optimisation depends on frequency
								case 2.17
									m.fitRog.coef = [100, .6, .03, .2, 160, 0, 0]; % initial
									m.add.funFun = @(d, m) plot(d.(m.x), d.(m.y));
								case 8.46
									m.fitRog.coef = [150, .6, .03, .4, 0, -180, 0]; % initial
								case 16.92
									m.fitRog.coef = [100, .6, .03, .4, 160, 0, 0]; % initial
							end
						case 'Croner'
							m.fitRog.coef = [40, .6, .05, .5, 0, 180, 0]; % initial
							m.fitRog.fix = {'pCen', 'pSur', 'disp'}; % fixed coef.
						case 'Wool'
							m.fitRog.coef = [45, .15, .1, 1, -180, 0, 0]; % initial
							m.fitRog.fix = {'rSur'}; % fixed coef.
						case 'Yeh' % not applicable: DOG doesn't do temporal frequency resp.
					end
				case 'rog' % ratio of Gaussians
					m.fitRog.name = {'kCen', 'kSur', 'rCen', 'rSur', ...
						'tau', 'disp', 'del', 'kA'}; % coefficient names
					switch m.resp.source
						case 'Benardete'
							%	m.select.freqT = 2.17; % 2.17, 8.46, 16.92
							m.fitRog.coef = [-150, 1, .03, .2, .014, 0, 0, 2]; % initial
							m.fitRog.fix = {'disp'}; % fixed coefficients
						case 'Croner' % *** fix ***
							m.fitRog.coef = [-245, 1.6, .08, 1, .014, .03, 0, 7]; % initial
							m.fitRog.fix = {'rSur', 'tau', 'del'}; % fixed coefficients
						case 'Wool'
							m.fitRog.coef = [-245, 1.6, .08, 1, .014, .03, 0, 7]; % initial
							m.fitRog.fix = {'rSur', 'tau', 'del'}; % fixed coefficients
						case 'Yeh'
							m.fitRog.coef = [4, 0, 0, 0, .01, 0, .01, 1]; % initial
							m.fitRog.fix = {'kSur', 'rCen', 'rSur', 'disp'}; % fixed
					end
			end
			m.refine.zero = .1; % replace 0, if it exists, with this value
			m.refine.n = 30; % number of predictor values
		case 'resp full' % predict ganglion cell output for full model (mV)
			m.p.ecc = 10; % eccentricity (deg)
			b = m.rogPar(m); % ROG parameters
			x = [5, 2]; % 5 cycles/deg, 2 Hz
			c = .3; % contrast
			a = abs(m.calRog(b, x) * c); % response amplitude (Hz)
			fprintf('ROG response amplitude: %.3g Hz\n', a); % list
			return
		case 'resp.freqT' % cone temporal frequency response: Baudin (19)
			m.tasks = 'respCone select plot set';
			%	m.select.source = 'm50';
			m.plot.group = 'source'; m.plot.arg = {'o-'};
			m.x = 'freq'; m.y = 'sens';
			m.set.axes = {'xScale', 'log', 'xLim', [1, 50], 'yScale', 'log', ...
				'clipping', 'off'};
		case 'resp.freqT fit' % fit cone temporal frequency response: Baudin (19)
			m.tasks = 'respCone select plot set fitCone show pred add';
			m.select.source = 'm5';
			m.plot.group = 'source'; m.plot.arg = {'o-'};
			m.x = 'freq'; m.y = 'sens';
			m.set.axes = {'xScale', 'log', 'xLim', [1, 50], 'yScale', 'log'};
			m.fitnlm.funFun = @cone; m.fitnlm.beta0 = [.9, .02];
			m.fitnlm.arg = {'predictorVars', {'freq'}, 'responseVar', 'sens'};
		case 'resp.t' % plot pulse response: Lee (94), Sinha (17)
			m.tasks = 'respPulse plot set'; % plot empirical data
			m.resp.source = 'Sinha'; % data source: Lee, Sinha
			switch m.resp.source
				case 'Lee' % Lee (94)
					m.plot.group = 'dur'; m.x = 'time'; m.y = 'resp';
					m.set.axes = {'xLim', [0, .2], 'xTick', [0, .1, .2], ...
						'yLim', [-100, 100], 'yTick', [-100, -50, 0, 50, 100], ...
						'clipping', 'off'};
				case 'Sinha' % Sinha (17)
					m.tasks = 'respPulse prepPulse plot set';
					m.prepPulse.group = 'type';
					m.plot.line = 'type'; m.x = 'time'; m.y = 'resp';
					%	m.plot.colourOrder = 'LM'
					m.set.axes = {'xLim', [-.05, .15], 'xTick', [0, .05, .1], ...
						'yLim', [-1, .1], 'yTick', [-1, -.5, 0], 'clipping', 'off', ...
						'colorOrder', eye(3)};
			end
		case 'resp.t fit' % plot fitted pulse response: Lee (94)), Sinha (17)
			m.tasks = 'respPulse select plot set fitPulse show add';
			m.resp.source = 'Lee'; % Lee or Sinha
			switch m.resp.source
				case 'Lee'
					m.plot.group = 'dur'; % m.plot.arg = {'o'};
					m.x = 'time'; m.y = 'resp';
					m.set.axes = {'xLim', [0, .2], 'xTick', [0, .1, .2], ...
						'yLim', [-140, 140], 'yTick', [-100, 0, 100], ...
						'clipping', 'off'};
					m.fitRog.name = {'kCen', 'kSur', 'rCen', 'rSur', 'tau', 'disp', ...
						'del', 'kA'}; % coefficient name
					m.fitRog.coef = [100, 2, 0, 0, .014, 0, .01, 20]; % intial value of coef.
					m.fitRog.fix = {'rCen', 'rSur', 'disp', 'kA'};
						% fixed coefficient
				case 'Sinha'
					m.tasks = ...
						'respPulse prepPulse select plot set fitPulse show add';
					m.prepPulse.group = 'type';
					%	m.select.type = 'M';
					m.plot.line = 'type'; % m.plot.arg = {'o'};
					m.x = 'time'; m.y = 'resp';
					m.set.axes = {'xLim', [0, .15], 'xTick', [0, .05, .1, .15], ...
						'yLim', [-1.2, .2], 'yTick', [-1, -.5, 0], 'clipping', 'on'};
					m.fitRog.name = {'kCen', 'kSur', 'rCen', 'rSur', 'tau', 'disp', ...
						'del', 'kA', 'ss'}; % coefficient name
					m.fitRog.coef = [-.3, 0, 0, 0, .04, 0, .02, 2.5, 1];
						% initial value of coefficient
					%	m.fitRog.coef = [-.3, 0, 0, 0, .04, 0, .02, 5, 2];
						% 2 cone stages
					m.fitRog.fix = {'kSur', 'rCen', 'rSur', 'disp', 'ss'};
						% fixed coefficient
			end
		case 'resp.x' % inverse transform ROG spatial frequency response
			m.tasks = 'makeTable respS plot set';
			m.x = 'x';
			m.y = 'respSReal'; % respSAmp, respSPhase, respSReal
			m.p.ecc = 10; % eccentricity (deg)
			m.makeTable.freqT = [2, 8, 16]; % temporal frequency (Hz)
			m.respS.group = 'freqT';
			m.respS.range = .3;
			lim = .5 * m.respS.range;
			x = {'xLim', lim * [-1, 1], 'xTick', lim * [-1, 0, 1]};
			switch m.y % set plot conditions
				case 'respSAmp'
					m.plot.line = 'freqT';
					y = {'yLim', [0, 1200], 'yTick', 1000 * [0, .5, 1]};
				case 'respSPhase'
					m.plot.line = 'freqT';
					y = {'yLim', [-180, 230], 'yTick', [-180, 0, 180]};
				case 'respSReal'
					m.plot.group = 'freqT';
					y = {'yLim', [-200, 1200], 'yTick', 1000 * [0, .5, 1]};
			end
			m.set.axes = [x, y];
		case 'resp.x fit' % calc. inverse transform of DOG by fitting inverse of ROG
			m.tasks = 'makeTable respS fitDogS show plot set'; % show model
			m.tasks = 'makeTable respS fitDogS plot set export';
			m.x = 'x';
			m.y = 'respSPhase'; % respSAmp, respSPhase, respSReal
			m.p.ecc = 10; % eccentricity (deg)
			m.makeTable.freqT = 2; % temporal frequency (Hz)
			m.fitDogS.name = {'kCen', 'kSur', 'rCen', 'rSur', 'pCen', 'pSur'};
				% coefficient names
			m.fitDogS.coef = [1000, .7, .0404, .0902, 0, 180]; % initial estimates
			m.fitDogS.fix = {'rCen', 'rSur'}; % fixed coefficients
			lim = .15;
			x = {'xLim', lim * [-1, 1], 'xTick', lim * [-1, 0, 1]};
			switch m.y % set plot conditions
				case 'respSAmp'
					m.plot.line = 'freqT';
					y = {'yLim', [0, 1200], 'yTick', 1000 * [0, .5, 1]};
				case 'respSPhase'
					m.plot.line = 'freqT';
					y = {'yLim', [-180, 230], 'yTick', [-180, 0, 180]};
				case 'respSReal'
					m.plot.group = 'freqT';
					y = {'yLim', [-200, 1200], 'yTick', 1000 * [0, .5, 1]};
			end
			m.set.axes = [x, y];
		case 'sens' % ganglion cell contrast sensitivity: Croner (95)
			s = 100 * .963; % cont. sens. (Hz/contrast-unit): Croner (95) Figure 9b
			s = s / .75; % convert from descending arm to peak of spatial freq. resp.
			s = s / m.p.kRect; % 17.8 mV/contrast-unit
			fprintf(['Ganglion cell contrast sensitivity, ', ...
				'm.p.kSens = %g mV/contrast-unit\n'], s); % display
			return
	end

	% Set up plotting for export
	if contains(m.tasks, 'export')
		m.set.axesLength = [5, 5]; m.set.lineWidth = 1.5; % for uniform plotting
		m.export.arg = {'colorSpace', 'cmyk'}; % for file export
	end

	% Run tasks
	m.folder = [userpath, '/Data/Col']; % data folder
	d = table; % empty data table
	if isfield(m.p, 'file') % m.p.file is defined
		load([m.folder, filesep, m.p.file], 'd'); % load data
	end
	stream(d, m); % run tasks
