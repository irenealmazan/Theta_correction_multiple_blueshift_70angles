% test script for lsq.m (testpeak.m)% 5-3-94 SMB (Bren@SLAC.stanford.edu)x_array=[-2:.01:2];% elements of fit_par: center, y(center), fwhm and fraction% lorentzian for each peak, followed by intercept and slope% for backgroundfit_par= [-.5 10 .1 .5 ...		% First peak		 0 10 .2 .3 ...		% Second peak		.5 10 .4  0 ...		% Third peak		 1 .1];					% background intercept and slopey_array=pk_qvt(x_array,fit_par);y_array=y_array+ 0.5*rand(size(x_array));lsq(x_array,y_array);