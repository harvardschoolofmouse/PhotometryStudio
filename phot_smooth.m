function sts = phot_smooth(obj, ts, OnOff, method)
	% 
	% OnOff = 0: no smoothing
	% OnOff > 0: the kernel used for smoothing (eg # of samples. For 1khz sampled photometry data we recommend 70-100 samples for data presentation)
	% 
	if nargin < 4
		method = 'gausssmooth';
	end
	if nargin < 3 || OnOff < 0
		OnOff = 100;
    end
    
    if isempty(ts), sts = []; return, end

	if strcmp(method, 'gausssmooth')
		if OnOff
			sts = gausssmooth(ts, round(OnOff), 'gauss');
		else
			sts = ts;
		end
	else
		if OnOff
			sts = smooth(ts, round(OnOff), 'moving');
		else
			sts = ts;
		end
	end
end