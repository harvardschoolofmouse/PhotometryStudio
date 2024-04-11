function dFF = normalizedBaselineDFF(rawF, baselineWin, trialStartPositions, setBaselineZero, noiseKillingWindow, stdThresh)
	if nargin < 4, setBaselineZero = false;end
	if nargin < 5, noiseKillingWindow=1000;end
	if nargin < 6, stdThresh=15;end
	warning('did you mean to use the multi-trial verrsion of this? this method will distort your baseline if your trial outcomes are not randomized at the level of the experimental design!')	
	% 
	% 	Start by correcting away noise from rawF
	% 
	ts = rawF;
	% 1. Smooth the whole day timeseries (window is 10 seconds)
	disp(['killing noise with ' num2str(noiseKillingWindow) ' sample window, moving']);
	Fraw = smooth(ts, noiseKillingWindow, 'moving');
	% Fraw = smooth(SNc_values, 1000, 'gauss');
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
	% 
	% 	Determine the indicies of the beginnings of each baseline Period...
	% 
	baselineStart = trialStartPositions - baselineWin + 1;
	baselineStart(end+1) = numel(rawF); % tack on the full length so that we can correct the entire signal...
	% 
	% 	Initialize the baseline dFF
	% 
	nbDFF = ts;
	% 
	% 	Find the baseline indices and Run the dF/F correction...
	% 
	for iTrial = 1:numel(trialStartPositions)
		baselineIdxs_baseline((iTrial-1)*baselineWin + 1:iTrial*baselineWin) = baselineStart(iTrial):trialStartPositions(iTrial);
		% baselineIdxs_trial = baselineStart(iTrial):baselineStart(iTrial + 1)-1;
		F0 = nanmean(obj.GLM.rawF(baselineStart(iTrial):trialStartPositions(iTrial)));
		F = nbDFF(baselineStart(iTrial):baselineStart(iTrial + 1)-1);
		nbDFF(baselineStart(iTrial):baselineStart(iTrial + 1)-1) = (F - F0)/F0;
	end
	% 
	% 	If the baselines should all be set to exactly zero after dF/F, do that here:
	% 
	if setBaselineZero
		nbDFF(baselineIdxs_baseline) = 0;
		nbDFF(1:baselineIdxs_baseline(1)) = 0;
	end
	dFF = nbDFF;
	nbDFF = dFF;
end