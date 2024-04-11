function [dFF, style, edgeTrials] = normalizedMultiBaselineDFF(baselineWin, nTrials, ts,trialStartPositions, verbose,stdThresh)
	% 
	% 	Created: Allison Hamilos, ahamilos[at]wi.mit.edu; 2024
	%		git:HarvardSchoolofMouse/PhotometryStudio
	% 
	% 	Method described in Hamilos et al., eLife, 2021
	%
	%	This method establishes a baseline period for each trial (up to you--here it's configured for the self-time movement task)
	%		(baselineWin is the number of samples BEFORE the trial-start event of your choosing)
	%	It then will use as F0 the median of the current trial's baseline + the nTrials/2 nearest trials' baselines on either side of the current trial
	%	This is useful when normalizing to each trial's baseline would distort meaningful baseline info (e.g., when the animal gets to decide the trial outcome)
	% 
	% 	ts is the photometry timeseries you want to process
	% 	nTrials is the number of trials to average the baseline over (must be even #)
	%	verbose is true/false whether you want to see progress
	%	trialStartPositions is in samples relative to the ts vector
	% 
	%	The method first corrects for shot artifacts (use stdThresh to set your threshold for eliminating these by interpolation)
	% 		we recommend a threshold of 15 in headfixed animals, but check this against your own signals.
	%		using verbose=true will help you evaluate this
	%	The method then finds the baseline windows surrounding each trial and calculates the median to get F0
	% 	Finally, it gets dF/F
	%	
	% 
	if nargin < 6
        stdThresh = 15;
    end
    if nargin < 5
		verbose = false;
	end
	% 
	% 	Start by correcting away noise from rawF
	% 
	% 1. Smooth the whole day timeseries (window is 10 seconds)
	if verbose,disp('killing noise with 1,000 ms window, moving'); end
	Fraw = phot_smooth(ts, 1000, 'moving');
	% Fraw = smooth(SNc_values, 1000, 'gauss');
	if verbose,disp('noise killing complete');end

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
	%	We get F0 = median of n trials baselines
	% 
	% 	Then we normalize all the points in the trial to the n/2-trial baseline on either side (n+1 trials total!)
	% 
	% 	For the edges, we will just take median of however many trials are in range... Can exclude these if necessary
	% 
	% 	To be able to handle chop data, we need to exclude nans from the mean...
	% 
	if nargin < 3
		nTrials = 10;
	elseif rem(nTrials, 2)
		nTrials = nTrials - 1;
		warning(['nTrials must be even, because use n/2 trials on either side. nTrials reduced to ' num2str(nTrials)]);
	end
	% 
	% 	Determine the indicies of the beginnings of each baseline Period...
	% 
	baselineStart = trialStartPositions - baselineWin + 1;
	baselineStart(end+1) = numel(ts); % tack on the full length so that we can correct the entire signal...
	% 
	% 	Initialize the baseline dFF
	% 
	dFF = nan(size(ts));
	style = ['median | baseline window(ms) = ' num2str(baselineWin) ' | nTrials (not including nth trial) = ' num2str(nTrials), ' | gFit: normalized multibaseline | symmetric'];
	edgeTrials = [1:nTrials/2, numel(trialStartPositions) - nTrials/2:numel(trialStartPositions)];
	% 
	% 	Find the baseline indices and Run the dF/F correction... as well as the number of licks in range.
	% 
	F0s = nan(size(trialStartPositions));
	for iTrial = 1:numel(trialStartPositions)
		if iTrial == 1 && baselineStart(iTrial) <= 0
			nanpad = nan(1, 1-baselineStart(iTrial));
			baselineIdxs_baseline((iTrial-1)*baselineWin + 1:iTrial*baselineWin) = [nanpad, 1:trialStartPositions(iTrial)];
			F0s(iTrial) = nanmedian(ts(1:trialStartPositions(iTrial)));
		else
			baselineIdxs_baseline((iTrial-1)*baselineWin + 1:iTrial*baselineWin) = baselineStart(iTrial):trialStartPositions(iTrial);
			F0s(iTrial) = nanmedian(ts(baselineStart(iTrial):trialStartPositions(iTrial)));
		end
    end
    if verbose
        figure 
        plot(F0s, 'DisplayName', 'Median Baseline')
        movingAve = smooth(F0s, nTrials+1, 'moving');
        hold on, plot(movingAve, 'DisplayName', 'F0')
        xlabel('Trial #')
        ylabel('F (V)')
    end
	% 
	% 	Now normalize to multibaseline:
	% 
	for iTrial = 1:numel(trialStartPositions)
		if ismember(iTrial, edgeTrials)
			if iTrial < numel(trialStartPositions)/2
				nLeftSideTrials = iTrial-1;
				nRightSideTrials = nTrials/2;
			else
				nLeftSideTrials = nTrials/2;
				nRightSideTrials = numel(trialStartPositions) - iTrial;
			end
		else
			nLeftSideTrials = nTrials/2;
			nRightSideTrials = nTrials/2;
		end
		F0 = nanmean(F0s(iTrial-nLeftSideTrials:iTrial+nRightSideTrials));
		if iTrial == 1 && baselineStart(iTrial) <= 0
			nanpad = nan(1, 1-baselineStart(iTrial));					
			F = ts(1:baselineStart(iTrial + 1)-1);
			% 
			% 	Need to handle nans... this will happen automatically - any nan points will remain nans after this step. This is ok.
			% 
			dFF(1:baselineStart(iTrial + 1)-1) = (F - F0)/F0;	
		else
			F = ts(baselineStart(iTrial):baselineStart(iTrial + 1)-1);
			% 
			% 	Need to handle nans... this will happen automatically - any nan points will remain nans after this step. This is ok.
			% 
			dFF(baselineStart(iTrial):baselineStart(iTrial + 1)-1) = (F - F0)/F0;	
		end
	end
	% 
	% 	Set the pre-session with moving boxcar...
	% 
    % 11/7/2020: I think if pre-sesh is very long, this will cause
    % issue. So handle this case...
	boxwin = nTrials*20*1000*(numel(1:baselineStart(1)-1)<nTrials*20*1000) +(numel(1:baselineStart(1)-1)+ nTrials*20*1000)*(numel(1:baselineStart(1)-1)>nTrials*20*1000);
	preSesh = boxDFF(ts(1:boxwin*2), '200000');
	dFF(1:baselineStart(1)-1) = preSesh(1:baselineStart(1)-1);
    %
    %   Kill noise
    %
    dFF(dFF>stdThresh*nanstd(dFF)) = nan;
   
end
