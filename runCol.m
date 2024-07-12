function runCol

% Analyse the behaviour of the visual signal processing model, version "Colour".
% This file sets the metadata and initiates the stream of analysis tasks.

	% Initialise
	m = setColLit; % set metadata calculated from the literature
	m = setDef(m); % set default metadata
	switch 1 % set current analysis
		case 1
			m.p.ecc = 10;
			switch m.p.ecc % set visual field width
				case 1, m.p.wid = .3; m.p.freqSRange = [5, 15]; % lower, upper
				case 3, m.p.wid = .5; m.p.freqSRange = [4.6, 10.2];
				case 10, m.p.wid = 1; m.p.freqSRange = [1.8, 7];
				case 30, m.p.wid = 4;
			end
	end

	% Specify the analysis tasks to be performed
	switch 'array.x.y'
		case 'array.x.y' % plot cell locations and types
			m.tasks = 'map select desc plot set';
			%	m.tasks = 'map select desc plot set fig add'; % add mechanisms
			m.x = 'x'; m.y = 'y';
			m.select.array = {'gangOff', 'gangOn'}; m.select.array = 'cone';
			m.fig.mech = [];
			switch m.p.ecc % for Figure 9
				case 1, m.p.wid = .14;
				case 3, m.p.wid = .14;
					%	m.fig.mech = [.0125, .0212]; % maximum response to S-stimulation
				case 10, m.p.wid = .24; m.p.wid = .5; % for Figure 1
				case 30, m.p.wid = 1.4;
					%	m.fig.mech = [0, 0]; % nearest g.c. to innermost S-cone
			end
			m.plot.group = {'ecc', 'array'}; m.plot.group = 'ecc';
			m.plot.funFun = @(d, m)scatter(d, 'x', 'y', 'filled', 'colorVar', 'c', ...
				'sizeData', 10);
			w = .5 * m.p.wid; lim = w * [-1, 1]; tick = w * [-1, 0, 1];
			m.set.axes = {'xLim', lim, 'xTick', tick, 'yLim', lim, 'yTick', tick};
			m.add.funFun = @(d, m)plot(d.y, 'k');
		case 'comp.m.a' % plot cone components in centre and surround responses
			m.tasks = 'resp reduce fund unpack plot set';
			switch 'cen' % mechanism
				case 'cen', m.resp.stage = 'gangOn'; % centre
				case 'sur', m.resp.stage = 'back'; % surround
			end
			m.p.freqS = 0; % uniform stimulus field
			m.p.back = 0; % open loop
			m.resp.cont = .3 * eye(3); % L, M, and S cone contrasts
			m.reduce.loc = [0, 0]; % m.reduce.loc = [-.038, -.066];
			m.unpack.var = 'cont'; m.unpack.z = 'resp';
			m.x = 'resp'; m.y = 'resp'; m.plot.line = 'cont';
			m.plot.funFun = @(d, m)polarplot(squeeze(d.resp), 'o', ...
				'color', d.cont / max(d.cont), 'markerFaceColor', d.cont / max(d.cont));
			m.set.axes = {'rLim', [0, 7], 'rTick', [0, 3, 6]};
		case 'count.ratio' % plot histogram of cone components in centre and sur.
			m.tasks = 'resp fund comp crop unpack hist plot'; % unnormalised
			m.tasks = 'resp fund comp crop unpack hist prep plot set export';
			switch 'sur' % mechanism
				case 'cen', m.resp.stage = 'gangOn'; % centre
				case 'sur', m.resp.stage = 'back'; % surround
			end
			switch m.p.ecc
				case .3, m.p.wid = .8; % wid = 1 better (n = 3500) but takes long time
				case 1, m.p.wid = 1; m.p.align = 1; m.p.radGangCoef = 1e-5; % parafovea
				case 3, m.p.wid = 1.5;
				case 10, m.p.wid = 2.5;
				case 30, m.p.wid = 5;
			end
			m.p.freqS = 0; % uniform stimulus field
			m.p.back = 0; % open loop
			m.resp.cont = .3 * eye(3); % L, M, and S cone contrasts
			m.crop.z = 'resp';
			m.unpack.var = 'cont'; m.unpack.z = 'resp';
			m.hist.group = 'cont';
			m.hist.edges = linspace(0, 1, 11);
			m.prep.group = 'cont'; m.prep.norm = 'count';
			m.x = 'resp'; m.y = 'count'; m.z = 'count';
			m.plot.group = {'cont', 'ecc'};
			m.plot.funFun = @(d, m)plot(d.(m.x), d.(m.z), ...
				'color', d.cont / max(d.cont));
			m.set.axes = {'xLim', [0, 1], 'xTick', [0, .5, 1], ...
				'yLim', [0, 1], 'yTick', [0, .5, 1]};
		case 'count.resp' % plot histogram of maximal response
			m.tasks = 'resp fund doMax hist desc plot set';
			m.p.pot = 0; m.p.cont = m.p.contMag * [1, -1, 0];
			dir = linspace(-180, 180, 20 + 1); m.resp.dir = dir(1: end - 1)';
			f = m.p.freqSRange; f = log10(f);
			freqS = logspace(f(1), f(2), 10); m.resp.freqS = freqS';
			m.resp.stage = 'gangOff';
			m.fund.prop = 'amp'; % amp, complex, phase
			m.doMax.z = 'resp';
			m.x = 'resp'; m.y = 'count'; m.z = 'count';
			m.plot.group = 'ecc';
			%	m.set.axes = {'xLim', lim, 'xTick', tick, 'yLim', lim, 'yTick', tick};
		case 'count.resp cont' % plot histogram of maximal resp., multiple contrasts
			m.tasks = 'resp fund unpack doMax crop hist prep desc plot set';
			m.p.pot = 0;
			m.resp.cont = m.p.contMag * [1, -1, 0; 1, 1, 0];
			dir = linspace(-180, 180, 20 + 1); m.resp.dir = dir(1: end - 1)';
			f = m.p.freqSRange; % f = [.1, 7];
			f = log10(f); f = logspace(f(1), f(2), 10); m.resp.freqS = f';
			m.resp.stage = 'gangOff';
			m.fund.prop = 'amp';
			m.unpack.var = 'cont';
			m.doMax.group = 'cont';
			m.crop.group = 'cont';
			m.hist.group = 'cont'; m.hist.z = 'count';
			m.prep.group = 'cont'; m.prep.norm = 1; m.prep.z = 'count';
			m.x = 'resp'; m.y = 'count'; m.z = 'resp';
			m.plot.group = {'ecc', 'total'}; m.plot.line = 'cont';
			lim = 25;
			m.set.axes = {'xLim', lim * [0, 1], 'xTick', lim * [0, .5, 1], ...
				'yLim', [0, 1], 'yTick', [0, .5, 1]};
		case 'count.resp pref' % plot histogram of preferred spatial frequency
			m.tasks = 'resp fund doMax crop hist prep desc plot';
			m.tasks = 'resp fund doMax interp crop hist prep desc plot';
			m.p.pot = 0;
			dir = linspace(-180, 180, 20 + 1); m.resp.dir = dir(1: end - 1)';
			f = m.p.freqSRange;
			%	f = log10(f); f = logspace(f(1), f(2), 10); m.resp.freqS = f';
			f = linspace(f(1), f(2), 10); m.resp.freqS = f';
			m.resp.stage = 'gangOff';
			m.fund.prop = 'amp'; % amp, complex, phase
			m.doMax.z = 'resp'; m.doMax.out = 'freqSPref'; m.doMax.out = 'freqSTun';
			m.interp.x = 'freqS'; m.interp.z = 'freqSTun';
			m.interp.out = 'freqSPref';
			m.crop.z = 'freqSPref';
			m.hist.edges = 10;
			m.prep.norm = 1;
			m.x = 'freqSPref'; m.y = 'count'; m.z = 'count';
			m.plot.group = {'ecc', 'total'};
			m.set.axes = {'xScale', 'log', 'xLim', [1, 10], 'xTick', [1, 3, 10], ...
				'yLim', [0, 1], 'yTick', [0, .5, 1]};
			m.set.axes = {'xLim', [0, 10], 'xTick', [0, 5, 10], ...
				'yLim', [0, 1], 'yTick', [0, .5, 1]};
		case 'count.wid' % plot histogram of direction bandwidth
			%	*** fix ***
			m.tasks = 'resp fund doMax unpack interpWid hist desc plot set';
			m.p.pot = 0;
			dir = linspace(-180, 180, 20 + 1); m.resp.dir = dir(1: end - 1)';
			m.resp.freqS = linspace(1, 10, 10)'; m.resp.freqS = 3.6;
			m.resp.stage = 'gangOff';
			m.fund.prop = 'amp'; % amp, complex, or phase
			m.doMax.z = 'resp'; m.doMax.out = 'dirTun';
			m.unpack.z = 'dirTun'; m.unpack.var = 'loc';
			m.interp.group = 'loc'; m.interp.z = 'dirTun'; m.interp.method = 'spline';
			m.x = 'dir'; m.y = 'resp'; m.z = 'resp';
			m.plot.group = {'loc', 'stage'};
			m.plot.funFun = @(d, m)plot(d.(m.x), squeeze(d.dirTun), 'o');
			m.set.axes = {'xLim', [-180, 180], 'xTick', [-180, 0, 180]};
		case 'fig' % plot figure
			switch 'conv'
				case 'conv' % convergence functions
					m.tasks = 'fig plot set';
					lim = .2;
					m.fig.conv = 1; m.fig.xLim = lim;
					m.x = 'x'; m.y = 'f'; m.plot.group = 'source';
					m.set.axes = {'xLim', lim * [-1, 1], 'xTick', ...
						lim * [-1, 0, 1], 'yLim', [0, 1], 'yTick', [0, 1]};
				case 'course' % stimulus time course *** fix ***
					m.tasks = 'stim fig plot set';
					m.fig.course = 1;
					m.x = 't'; m.y = 'course';
					yLim = m.p.cont(1);
					m.set.axes = {'xLim', [0, .5], 'xTick', [0, .25, .5], ...
						'yLim', yLim * [-1, 1], 'yTick', yLim * [-1, 0, 1]};
			end
		case 'rad.ecc' % calculate radius from spatial frequency response
			m.tasks = 'resp reduce fund radius save list'; % save radius data
			m.tasks = 'load plot set'; % plot it
			m.p.ecc = 3; % save one eccentricity at a time
			m.p.back = 0;
			switch m.p.ecc
				case .1, m.p.wid = .2; xLim = log10([2, 50]);
				case .3, m.p.wid = .2; xLim = log10([2, 50]);
				case 1, m.p.wid = .5; xLim = log10([2, 50]);
				case 3, m.p.wid = .5; xLim = log10([2, 50]);
				case 10, m.p.wid = 1; xLim = log10([.5, 12]);
				case 30, m.p.wid = 2.5; xLim = log10([.1, 10]);
			end
			switch 'cen' % choose mechanism
				case 'cen' % centre mechanism
					m.resp.stage = 'gangOff'; % centre
					yLim = [.007, .3]; yTick = [.01, .1];
				case 'sur' % surround mechanism
					m.resp.stage = 'back'; % surround
					yLim = [.007, 5]; yTick = [.01, .1, 1];
			end
			freqS = logspace(xLim(1), xLim(2), 30); freqS = [0, freqS];
			m.resp.freqS = freqS';
			m.reduce.loc = [0, 0];
			m.fund.prop = 'amp'; % amp, complex, phase		
			m.save.name = 'Radius.mat'; m.load.name = m.save.name;
			m.x = 'ecc'; m.y = 'radius'; m.plot.arg = {'o'};
			m.plot.group = 'stage';
			m.set.axes = {'xScale', 'log', 'xLim', [.006, 40], ...
				'xTick', [.01, .1, 1, 10], 'yScale', 'log', 'yLim', yLim, ...
				'yTick', yTick, 'clipping', 'off'};
		case 'resp.dir' % plot direction tuning with optional interpolated fit
			switch 'none' % add interpolated fit?
				case 'fit' % yes
					m.tasks = 'resp reduce fund doMax desc plot set interp add';
					m.plot.funFun = @(d, m)plot(d.(m.x), squeeze(d.dirTun), 'o');
				case 'none' % no
					m.tasks = 'resp reduce fund doMax desc plot stop set';
					m.plot.funFun = @(d, m)plot(d.dir, squeeze(d.dirTun));
			end
			m.p.pot = 0;
			dir = linspace(-180, 180, 20 + 1); m.resp.dir = dir(1: end - 1)';
			f = m.p.freqSRange; m.resp.freqS = linspace(f(1), f(2), 10)';
			m.resp.stage = 'gangOff';
			m.reduce.loc = [0, 0]; % m.reduce.loc = [-.01, -.071];
			m.fund.prop = 'amp'; % amp, complex, or phase
			m.doMax.out = 'dirTun';
			m.x = 'dir'; m.y = 'resp'; m.z = 'resp';
			m.plot.group = {'loc', 'stage'}; yLim = 5;
			m.set.axes = {'xLim', [-180, 180], 'xTick', [-180, 0, 180], ...
				'yLim', yLim * [0, 1], 'yTick', yLim * [0, .5, 1], 'clipping', 'off'};
			m.interp.z = 'dirTun';
			m.add.y = 'dirTun';
			m.add.funFun = @(d, m)plot(d.(m.x), squeeze(d.dirTun));
		case 'resp.dir cont' % plot direction tuning with variable contrast
			m.tasks = 'resp reduce fund unpack doMax desc plot stop set';
			m.p.pot = 0;
			m.resp.cont = m.p.contMag * [1, -1, 0; 1, 1, 0];
			dir = linspace(-180, 180, 20 + 1); m.resp.dir = dir(1: end - 1)';
			f = m.p.freqSRange; m.resp.freqS = linspace(f(1), f(2), 10)';
			m.resp.stage = 'gangOff';
			m.reduce.loc = [0, 0]; m.reduce.loc = [-.1, -.07];
			m.fund.prop = 'amp'; % amp, complex, or phase
			m.unpack.var = 'cont';
			m.doMax.group = 'cont'; m.doMax.out = 'dirTun';
			m.x = 'dir'; m.y = 'resp'; m.z = 'resp';
			m.plot.group = {'loc', 'stage', 'cont'};
			m.plot.funFun = @(d, m)plot(d.dir, squeeze(d.dirTun));
			m.set.axes = {'xLim', [-180, 180], 'xTick', [-180, 0, 180], ...
				'yLim', [0, 5], 'yTick', 5 * [0, .5, 1]};
		case 'resp.dir.freqS' % plot direction and spatial frequency tuning
			m.tasks = 'resp reduce fund unpack desc plot set';
			m.p.pot = 0;
			%	m.p.cont = m.p.contMag * [1, -1, 0];
			%	m.resp.cont = .3 * [-1, -1, 1; 1, -1, 0; 1, 1, 0; 1, 1, 1];
			m.resp.dir = linspace(-180, 180, 20 + 1)';
			f = [0, 10]; % f = m.p.freqSRange;
			m.resp.freqS = linspace(f(1), f(end), 11)';
			m.resp.stage = 'gangOff';
			m.reduce.loc = [0, 0]; % m.reduce.loc = [.069, -.189];
			m.fund.prop = 'amp'; % amp, complex, or phase
			m.unpack.var = {'cont', 'loc'}; m.unpack.var = 'loc';
			m.x = 'dir'; m.y = 'freqS'; m.z = 'resp';
			%	m.plot.group = {'loc', 'stage', 'cont'};
			m.plot.group = {'loc', 'stage'};
			m.plot.funFun = @(d, m)imagesc(d.dir, d.freqS, squeeze(d.resp)');
			m.set.axes = {'xLim', [-180, 180], 'xTick', [-180, 0, 180], ...
				'yLim', f, 'yTick', [0, 5, 10]};
			m.set.colorbar = {};
		case 'resp.freqS' % plot spatial frequency tuning
			switch 'no' % add interpolated fit?
				case 'yes'
					m.tasks = 'resp reduce fund doMax desc plot set interp add';
					m.plot.funFun = @(d, m)plot(d.(m.x), squeeze(d.freqSTun), 'o');
				case 'no'
					m.tasks = 'resp reduce fund doMax desc plot stop set';
					m.plot.funFun = @(d, m)plot(d.(m.x), squeeze(d.freqSTun));
			end
			m.p.pot = 0;
			dir = linspace(-180, 180, 20 + 1); m.resp.dir = dir(1: end - 1)';
			%	f = m.p.freqSRange; m.resp.freqS = linspace(f(1), f(2), 10)';
			f = [1, 13]; f = log10(f); m.resp.freqS = logspace(f(1), f(2), 10)';
			m.resp.stage = 'gangOff';
			m.reduce.loc = [0, 0]; % m.reduce.loc = [.052, -.137];
			m.fund.prop = 'amp'; % amp, complex, phase
			m.doMax.out = 'freqSTun';
			m.x = 'freqS'; m.y = 'resp'; m.z = 'resp';
			m.plot.group = {'ecc', 'loc', 'stage'};
			m.set.axes = {'xScale', 'log', 'xLim', [1, 13], 'xTick', [1, 3, 10], ...
				'yScale', 'log', 'yLim', [1, 4], 'yTick', [1, 4]};
			%	m.set.axes = {'xLim', [0, 10], 'xTick', [0, 5, 10], ...
			%	'yLim', [0, 4.5], 'yTick', [0, 2, 4]};
			m.interp.z = 'freqSTun';
			m.add.y = 'freqSTun';
			m.add.funFun = @(d, m)plot(d.(m.x), squeeze(d.freqSTun));
		case 'resp.freqS fix' % plot spatial frequency tuning with fixed directions
			m.tasks = 'resp reduce fund prep unpack desc plot set'; % uncentred
			m.tasks = 'resp reduce centre fund prep unpack desc plot set';
			%	m.p.cont = [1, 0, 0]; %	m.p.pot = 0;
			m.resp.dir = (0: 15: 90)'; m.resp.dir = 0;
			xLim = [1, 10]; % subcortex: [1, 13]; cortex: [1, 20]
			x = log10(xLim); % limits for log x
			freqS = logspace(x(1), x(2), 20); m.resp.freqS = freqS';
			m.resp.back = [0; 1]; m.resp.back = 1;
			m.resp.stage = 'back'; m.resp.stage = {'gangOff'};
			m.reduce.loc = [0, 0];
			m.fund.prop = 'amp'; % amp, complex, phase
			m.unpack.var = {'dir', 'loc', 'stage'};
			m.x = 'freqS'; m.y = 'resp'; m.z = 'resp';
			m.plot.group = 'stage';
			m.plot.line = {'dir', 'loc'};
			m.plot.funFun = @(d, m)plot(d.freqS, squeeze(d.resp)');
			a = {'xScale', 'log', 'xLim', [1, 10], ...
				'xTick', [1, 3, 10], 'clipping', 'off'};
			switch m.fund.prop % set axes
				case 'amp'
					m.set.axes = [a, 'yScale', 'log', 'yLim', [1, 13], ...
						'yTick', [1, 3, 10]];
					m.set.axes = [a, 'yScale', 'log', 'yLim', [1, 8], 'yTick', [1, 3, 8]];
				case 'phase'
					m.prep.unwrap = 'freqS';
					m.set.axes = [a, 'yLim', 200 * [-1, 1], ...
						'yTick', 180 * [-1, 0, 1], 'clipping', 'off'];
					m.set.axes = [a, 'yLim', [120, 180], ...
						'yTick', [120, 150, 180], 'clipping', 'off'];
			end
		case 'resp.freqT' % plot temporal frequency response
			m.tasks = 'resp reduce fund unpack desc plot set';
			%	m.p.pot = 0;
			m.resp.stage = 'gangOff';
			xLim = log10([1, 50]); % limits for log frequency
			freqT = logspace(xLim(1), xLim(2), 30); m.resp.freqT = freqT';
			m.reduce.loc = [0, 0];
			m.fund.prop = 'amp'; % amp, complex, phase
			m.unpack.var = 'stage'; m.unpack.z = 'resp';
			m.x = 'freqT'; m.y = 'resp';
			m.plot.group = 'stage';
			m.plot.funFun = @(d, m)plot(d.freqT, squeeze(d.resp));
			m.set.axes = {'xScale', 'log'}; % 'yLim', [0, 6]};
		case 'resp.t' % plot model time course for drifting grating
			m.tasks = 'resp reduce unpack prep desc plot stop set';
			m.p.domain = 'time'; m.p.solver = 'solveF';
			%	m.p.solver = 'solveT'; m.p.ts = 128;
			m.p.pot = 0;
			m.p.dir = -126; m.p.freqS = 4.8; m.p.cont = m.p.contMag * [1, -1, 0];
			m.resp.stage = 'gangOff';
			m.reduce.loc = [0, 0]; m.reduce.loc = [-.1, -.07];
			m.unpack.var = {'stage', 'loc'};
			m.prep.real = 1;
			m.x = 'time'; m.y = 'resp'; m.z = 'resp';
			m.plot.group = 'loc';
			m.plot.line = {'stage'};
			m.plot.funFun = @(d, m)plot(d.time, squeeze(d.resp));
			switch m.p.solver, case 'solveT', off = .5; otherwise, off = 0; end
			lim = 7; x = {'xLim', off + .5 * [0, 1], 'xTick', off + .5 * [0, .5, 1]};
			if m.p.pot % potential
				y = {'yLim', off + lim * [-1, 1], 'yTick', off + lim * [-1, 0, 1]};
			else % impulse rate
				y = {'yLim', lim * [0, 1], 'yTick', lim * [0, .5, 1]};
			end
			m.set.axes = [x, y]; % m.set.axes = x;
		case 'resp.t flash' % plot model time course for flashed grating
			m.tasks = 'resp reduce unpack prep desc plot set';
			m.p.stim = 'flash'; % m.p.pot = 0;
			m.p.time = .16; m.p.domain = 'time'; m.p.solver = 'solveT';
			m.resp.stage = 'gangOn';
			m.reduce.loc = [0, 0];
			m.unpack.var = [{'loc'}, {'stage'}, fieldnames(m.resp)'];
			m.unpack.z = 'resp';
			m.x = 'time'; m.y = 'resp';
			m.prep.real = 1;
			m.plot.group = {'loc'};
			m.plot.line = {'stage'};
			m.plot.funFun = @(d, m)plot(d.time, squeeze(d.resp));
			m.set.axes = {'xLim', .15 * [0, 1], 'xTick', .15 * [0, .5, 1]};
			%	lim = 10; y = {'yLim', lim * [-1, 1], 'yTick', lim * [-1, 0, 1]};
		case 'resp.resp' % polar plot g.c. response to S-stim and lum-stim
			m.tasks = 'resp fund crop centre unpack plot set';
			%	m.p.ecc = 3; m.p.wid = .35; rLim = [0, 23]; rTick = 10;
			%	m.p.ecc = 10; m.p.wid = .5;
			m.p.ecc = 30; m.p.wid = 2; rLim =[0, 12]; rTick = 5;
			m.p.freqS = 0;
			m.p.pot = 0;
			m.resp.stage = 'gangOff';
			m.resp.cont = m.p.contMag * [0, 0, 1; 1, 1, 0]; % S-, luminance-stimulation
			m.crop.z = 'resp';
			m.unpack.var = 'cont';
			m.x = 'resp'; m.y = 'resp'; m.z = 'resp';
			m.plot.line = {'stage', 'cont'};
			m.plot.funFun = @(d, m)polarplot(squeeze(d.resp), 'o');
			m.set.axes = {'rLim', rLim, 'rTick', rTick * [1, 2]};
		case 'resp.resp cart' % g.c. response to S-stim vs lum-stim response
			m.tasks = 'resp fund crop unpack plot set';
			m.p.ecc = 3; m.p.wid = .35;
			m.p.ecc = 10; m.p.wid = .5;
			m.p.freqS = 0;
			m.resp.stage = 'gangOff';
			m.resp.cont = .3 * [0, 0, 1; 1, 1, 0]; % S-, luminance-stimulation
			m.fund.prop = 'amp';
			m.crop.z = 'resp';
			m.unpack.var = 'cont'; m.unpack.z = 'resp';
			m.x = 'resp'; m.y = 'resp'; m.plot.group = 'stage';
			m.plot.funFun = @(d, m)plot(squeeze(d.resp(2, 1, :)), ...
				squeeze(d.resp(1, 1, :)), 'o');
			m.set.axes = {'xLim', [0, 6], 'xTick', [0, 3, 6], ...
				'yLim', [0, 6], 'yTick', [0, 3, 6]};
		case 'resp.x.y' % plot receptive field: cross correlation of stim., response
			stage = 'gangOff'; m.save.name = ['Field ', stage, '.mat'];
			switch 'all' % all, plot, or save
				case 'all'
					m.tasks = 'resp field reduce unpack prep desc plot stop set';
				case 'plot'
					m.tasks = 'reduce unpack prep desc plot set';
					m.p.file = m.save.name;
				case 'save'
					m.tasks = 'select resp field save';
			end
			m.p.pot = 0;
			m.p.stim = 'flash';
			m.p.domain = 'time'; m.p.solver = 'solveT'; m.p.time = .16;
			%	m.p.cont = m.p.contMag * [1, -1, 0];
			dir = linspace(-90, 90, 10 + 1); m.resp.dir = dir(1: end - 1)';
			f = m.p.freqSRange; m.resp.freqS = linspace(f(1), f(2), 5)';
			%	m.resp.freqS = [4, 6, 8, 10]'; m.resp.freqS = 2.1;
			m.resp.phaseS = [-180, -90, 0, 90]';
			m.resp.stage = stage;
			m.reduce.loc = [0, 0]; % m.reduce.loc = [-.1, -.07];
			m.unpack.var = 'loc';
			m.prep.group = 'loc'; m.prep.image = 'resp'; m.prep.x = 'locS';
			m.x = 'x'; m.y = 'y'; m.z = 'resp'; m.plot.group = 'loc';
			m.plot.funFun = @(d, m)imagesc(d.x, d.y, shiftdim(d.resp, 1));
			w = .35 * m.p.wid; lim = [- w, w]; tick = [- w, 0, w];
			m.set.axes = {'xLim', lim, 'xTick', tick, 'yLim', lim, 'yTick', tick};
			limC = .3; m.set.caxis = limC * [-1, 1];
			m.set.colorbar = {'xLim', limC * [-1, 1], 'xTick', limC * [-1, 0, 1]};
			m.set.colormap = 'turbo';
		case 'stim.x.y' % show movie of stimulus
			m.tasks = 'stim desc plot set'; % first frame only
			%	m.tasks = 'stim desc plot set showMovie'; % all frames
			m.p.cont = .5 * [1, -1, 0]; m.p.dir = 0;
			m.plot.funFun = @(d, m) image(d.x, d.y, ...
				shiftdim(d.stim(1, 1, :, :, :), 2));
			m.x = 'x'; m.y = 'y';
			w = .5 * m.p.wid; % half-width of visual field
			lim = [- w, w]; tick = [- w, 0, w]; % axis limits and ticks (deg)
			m.set.axes = {'xLim', lim, 'xTick', tick, 'yLim', lim, 'yTick', tick};
	end

	% Set up plotting for export
	if contains(m.tasks, 'export')
		m.set.axesLength = [5, 5]; m.set.lineWidth = 1.5; % for uniform plotting
		m.export.arg = {'colorSpace', 'cmyk'}; % for file export
	end
	
	% Store function handles
	m.setDef = @setDef; % setDef function handle
	m = runColTask(m); % task function handles
	m = runColFun(m); % other function handles
	
	% Initiate analysis
	[d, m] = m.readFile(m); % set the data folder and read or create the data file
	m = m.setType(m); % create list of all cell types
	m = m.calStruct(m); % calculate model structure
	m = m.calTemp(m); % calculate temporal parameters
	w = warning; % save the current state of warnings
	stream(d, m); % execute the task list
	warning(w); % restore the original state of warnings

function m = setDef(m) % set default metadata

	% Stimulus parameters
	m.p.contMag = .3; % contrast magnitude
	m.p.cont = m.p.contMag * [1, 1, 1]; % cone contrast, [cL, cM, cS]
	m.p.dir = 0; % stimulus direction (deg anticlockwise from rightward)
	m.p.dur = .05; % stimulus duration (s)
	m.p.freqS = 4; % spatial frequency (cycles/deg)
	m.p.freqT = 2; % temporal frequency (cycles/s)
	m.p.locS = 0; % displacement in stimulus direction (deg)
	m.p.phaseS = 0; % spatial phase (deg)
	m.p.stim = 'drift'; % stimulus type
	
	% Structural parameters
	m.p.array = {'cone', 'gangOff', 'gangOn'}; % neuron arrays
	m.p.cone.stage = {'cone', 'hor', 'back', 'biOff', 'biOn'}; % cone stages
	m.p.gangOff.stage = {'gangOff'}; % off-centre stages
	m.p.gangOn.stage = {'gangOn'}; % off-centre stages

	% Computational parameters
	m.p.align = 0; % alignment of gang. cells with cones: 0 for random, 1 align
	m.p.back = 1; % feedback loop: 0 for open, 1 for closed
	m.p.domain = 'freq'; % response domain: freq or time
	m.p.ecc = 10; % functional eccentricity (deg)
	m.p.file = ''; % data file
	m.p.pot = 1; % signal type: 0 for impulse rate, 1 for potential
	m.p.project = 'Col'; % project
	m.p.seed = 0; % seed for random number generator
	m.p.solver = 'solveF'; % solver function: solveF or solveT
	m.p.time = 1; % simulation time (s)
	m.p.ts = 64; % number of sample times and temporal frequencies
	m.p.wid = .5; % visual field width and height (deg)
	m.p.widH = .5 * m.p.wid; % visual field half-width and half-height (deg)
