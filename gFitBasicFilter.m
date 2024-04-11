function [lp, f, dFF] = gFitBasicFilter(obj, ts, samples_per_ms, cutoff_f, noiseKillingWindow, stdThresh)
	% 
	% 	Low pass filter version -- NB, you will need to correct the sampling ratess for your application
	% 

	if nargin < 6, stdThresh = 15;end
	if nargin < 5, noiseKillingWindow = 1000;end
	if nargin < 4, cutoff_f = 1/20000; end
	fs = samples_per_ms*1000;
	% 
	% 	Start by correcting away noise from rawF
	% 
	ts = rawF;
	% 1. Smooth the whole day timeseries (window is 10 seconds)
	disp(['killing noise with ' num2str(noiseKillingWindow) ' sample window, moving']);
	Fraw = obj.smooth(ts, noiseKillingWindow, 'moving');
	disp('noise killing complete');

	% 2. Subtract the smoothed curve from the raw trace 			(Fraw)
	Fraw = ts - Fraw;

	% 2b. Eliminate all noise > 15 STD above the whole trace 		(Frs = Fraw-singularities)
		% find points > 15 STD above/below trace and turn to average of surrounding points
	ignore_pos = find(Fraw > stdThresh*nanstd(Fraw));
	% disp(['Ignored points SNc: ', num2str(ignore_pos)]);
	for ig = ignore_pos
		ts(ig) = mean([ts(ig-1), ts(ig+1)],2);
	end

	steepness = 0.95;

	[lp, f] = lowpass(ts,cutoff_f,fs, 'ImpulseResponse','iir','Steepness',steepness);
	dF = ts - lp;
	dFF = (ts - lp)./lp;

	
	figure, hold on
	plot([1:numel(ts)]./1000,ts, 'DisplayName', 'Raw Signal')
	plot([1:numel(lp)]./1000,lp, 'DisplayName', 'Low Pass');
	% plot([1:numel(dF)]./1000,dF, 'DisplayName', 'Detrended');
	legend
	xlabel('session time (s)')
	figure
	plot([1:numel(dFF)]./1000, dFF, 'DisplayName', 'dFF');
	legend

	fvtool(f, 1000, 1000)
	freqz(f, 1000, 1000)

	if setObj
		obj.gFitLP.dFF = dFF;
		obj.gFitLP.cutoff_f = cutoff_f;
		obj.gFitLP.steepness = steepness;
	end
end