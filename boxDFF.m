function [dFF] = boxDFF(ts, boxWin, smoothingWin, stdLevel)
% 
% 	Created: Allison Hamilos, ahamilos[at]wi.mit.edu; 2024
%		git:HarvardSchoolofMouse/PhotometryStudio
% 
% 	Method described in Hamilos et al., eLife, 2021
% 
% 	Uses a boxWin width (1/2 width on either side of current point) moving average to get F0
%	This method is very useful for datasets without clear trial-based structure (i.e., no obvious quiescent baseline period)
%	--fits the data well while introducing minimal signal processing artifacts
%	Disadvantages: it's very slow if your sampling rate is high
%					you should experiment with a variety of boxWindows given your sampling rate 
%					and frequency information in your signal of interest
%		-- for headfixed GCaMP6f, dlight and DA2m photometry, window of 200s works very well
%		NB: About 1/2 your window is the longest time over which your signals won't be distorted by the processing
%			this distortion window is a factor for any filtering method--we are essentially low-pass filtering and using that as F0
%			box method is pretty useful because you can interpret readily the time window over which your signals are minimally affected
%
%	smoothingWin: used for killing noise in F0:
%		we recommend 1000 (10s) for a 1kHz photometry dataset
%	stdLevel: used for eliminating shot noise:
%		we recommend stdLevel=15 for headfixed photometry
% ------------------------------------------------------------------

if nargin < 3, smoothingWin=1000;end
if nargin < 4, stdLevel=15;end
% 1. Smooth the whole day timeseries (window is 10 seconds)
disp(['killing noise with ' num2str(smoothingWin) ' sample window, moving']);
Fraw = smooth(ts, smoothingWin, 'moving');
disp('noise killing complete');

% 2. Subtract the smoothed curve from the raw trace (Fraw)
Fraw = ts - Fraw;

% 2b. Eliminate all noise > 15 STD above the whole trace 		(Frs = Fraw-singularities)
	% find points > 15 STD above/below trace and turn to average of surrounding points
ignore_pos = find(Fraw > stdLevel*std(Fraw));
Frs = ts;
for ig = ignore_pos
	Frs(ig) = mean([Frs(ig-1), Frs(ig+1)],2);
end


% 3. Repeat step 1 without the noise points (this way smoothing not contaminated with artifacts)
%																(Fsmooth)
disp(['Fitting boxDFF with ', boxWin, ' window, moving ', datestr(now,'HH:MM AM')]);
% Fsmooth = smooth(Frs, 500000, 'gauss');
Fsmooth = smooth(Frs, str2num(boxWin), 'moving');
disp('gfitting complete');

% 4. Now at each point, do the dF/F calculation: [Frs(point) - Fsmooth(point)]/Fsmooth(point)
dFF = (Frs - Fsmooth)./Fsmooth;