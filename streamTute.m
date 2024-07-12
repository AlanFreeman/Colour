function streamTute

% streamTute: tutorial in the use of stream

	% Specify tasks to be performed
	switch 'cos.time points'
		case 'cos.time axes' % plot cosine functions on different axes
			m.tasks = 'cos plot set'; % tasks
			m.plot.group = 'amp'; % grouping variable for axes
			m.x = 'time'; m.y = 'fun'; % plotting variables
			m.set.axes = {'xLim', [0, 6], 'xTick', [0, 3, 6], ...
				'yLim', [-2, 2], 'yTick', [-2, 0, 2]};
		case 'cos.time lines' % plot cosine functions on the same axes
			m.tasks = 'cos plot set'; % tasks
			m.plot.line = 'amp'; % grouping variable for lines
			m.x = 'time'; m.y = 'fun'; % plotting variables
			m.set.axes = {'xLim', [0, 6], 'xTick', [0, 3, 6], ...
				'yLim', [-2, 2], 'yTick', [-2, 0, 2]};
		case 'cos.time points' % plot cosine functions as points
			m.tasks = 'cos plot set'; % tasks
			m.plot.line = 'amp'; % grouping variable for lines
			m.x = 'time'; m.y = 'fun'; % plotting variables
			m.plot.arg = {'or'}; % plot symbols instead of lines
			m.set.axes = {'xLim', [0, 6], 'xTick', [0, 3, 6], ...
				'yLim', [-2, 2], 'yTick', [-2, 0, 2]};
	end
	
	% Define handles for task functions in this file
	m.cos.fun = @doCos;

	% Run the tasks
	d = table; % dummy table
	stream(d, m); % run the analysis tasks

function [d, m] = doCos(d, m)

% Calculate cosine functions

	% Calculate
	t = linspace(0, 2 * pi, 65); % time (s)
	t = t(1: end - 1); % open-ended interval
	f = cos(t); % cosine function
	f2 = 2 * f; % different amplitude
	
	% Store
	d.time = t; % time
	d = repmat(d, [2, 1]); % make second row
	d.amp = [1; 2]; % amplitude
	d.fun = [f; f2]; % function
