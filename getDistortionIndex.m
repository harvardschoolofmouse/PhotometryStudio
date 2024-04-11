function testdFFstruct = getDistortionIndex(rawF, baselineWin, simType, trialStartPositions,samples_per_ms, noiseKillingWindow,stdThresh)
	% 
	% 	You will need your raw fluorescence signal and whichever post-dF/F timeseries you want to examine
	% 
	error('this has not been implemented for general use yet. currently is configued for use with HSOM datasets. please contact ahamilos[at]wi.mit.edu if you want to try this')
	if nargin < 7, stdThresh=15;end
	if nargin < 6, noiseKillingWindow = 1000;end
	if nargin < 2, baselineWin=5000;end
	if ~isfield(obj.GLM, 'rawF'), obj.loadRawF, end
	if strcmpi(simType, 'nbDFF')
		% 
		% 	Simulate data
		% 
		simNB5000 = normalizedBaselineDFF(rawF, baselineWin, trialStartPositions, true, noiseKillingWindow, stdThresh);
		% 
		% 	Filter simulated data
		% 
		[testdFFstruct.lp, testdFFstruct.f, testdFFstruct.dFF] = obj.gFitBasicFilter(simNB5000+1, cutoff_f, samples_per_ms*1000);
		% 
		%	Filter with the box 200 
		% 
		box200simNB = FX_gfitdF_F_fx_roadmapv1_4(simNB5000+1, '200000');
		% 
		% 	Bin gfit200 data
		% 
		obj.getBinnedTimeseries(obj.GLM.gfit, 'outcome', 6);
		% 
		% 	Plot CTAs...
		% 
		figure,
		subplot(1,5,1)
		obj.plot('CTA', [3,5], true, obj.Plot.smooth_kernel, 'last-to-first', true);
		title('gFit 200s Boxcar Data')
		% 
		% 	Bin lp DFF
		% 
		obj.getBinnedTimeseries(dFF, 'outcome', 6);
		subplot(1,5,2)
		obj.plot('CTA', [3,5], true, obj.Plot.smooth_kernel, 'last-to-first', true);
		title(['gFit Basic Filter: fc: ' num2str(cutoff_f) 'Hz'])
		% 
		% 	Bin similated data
		% 
		obj.getBinnedTimeseries(obj.GLM.flush.nbDFF, 'outcome', 6);
		subplot(1,5,3)
		obj.plot('CTA', [3,5], true, obj.Plot.smooth_kernel, 'last-to-first', true);
		title('Simulated Data: 5 sec Normalized Baseline Method')
		% 
		% 	box200 dFF of simData
		% 
		% obj.getBinnedTimeseries(testdFFstruct.dFF, 'outcome', 6);
		obj.getBinnedTimeseries(box200simNB, 'outcome', 6);
		subplot(1,5,4)
		obj.plot('CTA', [3,5], true, obj.Plot.smooth_kernel, 'last-to-first', true);
		title(['gFit 200s Boxcar Simulation'])
		% 
		% 	Bin dFF of simData
		% 
		% obj.getBinnedTimeseries(testdFFstruct.dFF, 'outcome', 6);
		obj.getBinnedTimeseries(testdFFstruct.dFF, 'outcome', 6);
		subplot(1,5,5)
		obj.plot('CTA', [3,5], true, obj.Plot.smooth_kernel, 'last-to-first', true);
		title(['gFit LP Simulation fc: ' num2str(cutoff_f) 'Hz'])
		% 
		% 	Add simulated data to the gFitLP struct
		% 
		nbSimulatedDFF = testdFFstruct.dFF;
		nbSimulated_fc = cutoff_f;
		nbSimulated_baselineWidth = 5000;
	elseif strcmpi(simType, 'expBaseline')
		error('Not implemented')
		% 
		% 	Simulate data
		% 
		simData = obj.GLM.rawF;
		baselineWin = 3000;
		n = numel(obj.GLM.rawF);
		t = [0:n];
		exps = (exp(-1/(n/20)*[0:n])-1 + exp(-1/(n/4)*[0:n])-1+ exp(-1/(n/2)*[0:n])-1+0.5*exp(-1/(100*n)*[0:n])-1+6.3)/5;
		% 
		% 	Determine the indicies of the beginnings of each baseline Period...
		% 
		baselineStart = obj.GLM.pos.cue - baselineWin + 1;
		baselineStart(end+1) = numel(obj.GLM.rawF); % tack on the full length so that we can correct the entire signal...
		for iTrial = 1:numel(obj.GLM.pos.cue)
			baselineIdxs_baseline = baselineStart(iTrial):obj.GLM.pos.cue(iTrial);
			simData(baselineIdxs_baseline) = exps(baselineIdxs_baseline);	
		end
		figure, hold on
		plot(obj.GLM.rawF, 'DisplayName', 'RawF')
		plot(simData, 'DisplayName', 'simData')
		legend		
		% 
		% 	Filter simulated data
		% 
		if strcmpi(gfitType, 'lopass')
			[testdFFstruct.lp, testdFFstruct.f, testdFFstruct.dFF] = obj.gFitBasicFilter(simData, cutoff_f, fs);
		elseif strcmpi(gfitType, 'box200')
			testdFFstruct.dFF = FX_gfitdF_F_fx_roadmapv1_4(simData, '200000');
		end
		% 
		% 	Bin gfit200 data
		% 
		obj.getBinnedTimeseries(obj.GLM.gfit, 'outcome', 6);
		% 
		% 	Plot CTAs...
		% 
		figure,
		subplot(1,3,1)
		obj.plot('CTA', [3,5], true, obj.Plot.smooth_kernel, 'last-to-first', true);
		title('gFit 200 Boxcar Data')
		% 
		% 	Bin similated data
		% 
		obj.getBinnedTimeseries(simData, 'outcome', 6);
		subplot(1,3,2)
		obj.plot('CTA', [3,5], true, obj.Plot.smooth_kernel, 'last-to-first', true);
		title('Simulated Data: 3 sec Exponential Baseline Method')
		% 
		% 	Bin dFF of simData
		% 
		obj.getBinnedTimeseries(testdFFstruct.dFF, 'outcome', 6);
		subplot(1,3,3)
		obj.plot('CTA', [3,5], true, obj.Plot.smooth_kernel, 'last-to-first', true);
		if strcmpi(gfitType, 'lopass')
			title(['gFit Basic Filter: fc: ' num2str(cutoff_f) 'Hz'])
		elseif strcmpi(gfitType, 'box200')
			title(['gFit 200 Boxcar on simData'])
		end
	end

end