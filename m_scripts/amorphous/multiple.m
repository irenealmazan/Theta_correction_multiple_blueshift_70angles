function imult=multiple(k,coh,inc,energy,s,sample,mu,pol,geo)% MULTIPLE  Calculate multiple scattering for unpolarized x-ray beam%% imult=multiple(k,coh,inc,energy,s,sample,mu,pol,geo)%% Calculates the ratio of double scattering to single scattering,% of x-rays, for either an unpolarized incident beam or a linearly % polarized beam.%% Input:%        k      = array of scattering vector magnitudes (inv. Angstroms)%        coh    = independent coherent scattering intensity%        inc    = incoherent scattering intensity%        energy = incident x-ray energy (eV)%        s      = composition of sample (see CALSCAT for details)%        sample = [thickness density atden] where%                 thickness = sample thickness (cm)%                 density = sample density (g/cm^3)%                 atden = atomic density (atoms/�^3)%        mu     = linear absorption coefficient of sample (cm^-1)%        pol    = polarization constant (see CORINT). Two cases:%                 pol>=0 unpolarized beam (Warren's method)%                 pol<0  linearly polarized beam%        geo    = scattering geometry string:%                 'refl' = symmetric reflection, infinite thickness%                 'trans' = symmetric transmission%% Output:%        imult  = i2/i1 = ratio of double- to single- scattering%% Requires: leasqr.m, mult.mat, jfact.m, atomdata.m, strtoz.m, absorb.m% Todd Hufnagel (hufnagel@jhu.edu), Johns Hopkins University% 14-May-03 Combined old MULTIPLE and MULTSCAT into single routineglobal PROBEif isempty(PROBE)    warning('Global PROBE note set, assuming 0 (for x-rays)')    PROBE=0;endif PROBE~=0    warning('Global PROBE~=0! Are you sure you know what you''re doing?')end% Work with column vectorsk=k(:);coh=coh(:);inc=inc(:);if pol>=0       % Unpolarized beam    imult=warrenmozzi(k,energy,s,sample,mu,coh,inc,geo);else            % polarized beam    imult=multscat(k,energy,coh,sample,mu,geo);end    plot(k,imult)title('Ratio of double- to single- scattering')xlabel('q (inv. Angstroms)');ylabel('i2/i1');grid%--------------------------------------------------------------------------function i2i1=warrenmozzi(k,energy,s,sample,mu,coh,inc,geo)%% WARRENMOZZI  Calculate multiple scattering for unpolarized x-ray beam%% i2i1=warrenmozzi(k,energy,s,sample,mu,coh,inc,geo)%% Calculates the ratio of double scattering to single scattering,% of x-rays, for two scattering geometries. Based on the method% of Warren and Mozzi (see internal comments for references). Note% that these calculations are for an unpolarized incident beam.%% Input:%        k      = array of scattering vector magnitudes (inv. Angstroms)%        energy = incident x-ray energy (eV)%        s      = composition of sample, as in CALSCAT%        sample = [thickness density atden] where%                 thickness = sample thickness (cm)%                 density = sample density (g/cm^3)%                 atden = atomic density (atoms/�^3)%        mu     = linear absorption coefficient of sample (cm^-1)%        coh    = independent coherent scattering intensity%        inc    = incoherent scattering intensity%        geo    = scattering geometry string:%                 'refl' = symmetric reflection, infinite thickness%                 'trans' = symmetric transmission%% Output:%        i2i1  = i2/i1 = ratio of double- to single- scattering%% Copyright 1999, Todd Hufnagel, Johns Hopkins University (hufnagel@jhu.edu)% 08-25-99 TCH First version% 09-02-99 TCH Fixed bug in interpolation of transmission Q, added some error checking% 05-02-00 TCH Changed units on E from keV to eV, replace call to MYPRHO with call to AVGATOM% 06-10-02 TCH Minor changes to comments to reflect changes to other routines.% 12-May-03 TCH Renamed from MULTIPLE.M. Incorporated into new MULTIPLE that does both % polarized and unpolarized.%% Based on these references:% B.E. Warren and R.L. Mozzi, Acta Cryst. 21, 459-461 (1966)% C.W. Dwiggins, Jr. and D.A. Park, Acta Cryst. A27, 264-272 (1971)% C.W. Dwiggins, Jr., Acta Cryst. A28, 155-163 (1972)%% Requires: leasqr.m, mult.mat, jfact.m, atomdata.m, strtoz.m, absorb.m%% Note that this routine does not exactly replicate Warren's result (e.g. Fig. 10.12% in his book) but it's fairly close. Differences include whether or not you count% the incoherent scattering (he doesn't) and what polarization factor you put in% (here, none). But even taking these things into account, I could not exactly% replicate Warren's figure. global DISK_PATH DISK_DELIM AMORPHOUS_FILES% Work with column vectorsk=k(:);coh=coh(:);inc=inc(:);% Load the parameter tables from mult.matload([DISK_PATH AMORPHOUS_FILES DISK_DELIM 'mult.mat'])btrans=btrans(:);twothetatrans=twothetatrans(:);qtrans=qtrans(:);muttrans=muttrans(:);%% Parse the string containing the chemistry information%	[z,w] = strtoz(s);	nz = max(size(z));%% ww is the number of atoms in the formula unit%	ww=sum(w);%% Calculate two terms we'll need later%	asum=0;	sumz2=0;	for j=1:nz		[mr,rho,ifail]=absorb(z(j),energy);		m=atomdata(z(j)); 		asum=asum+w(j).*mr.*m.*10^4;		sumz2=sumz2+w(j).*z(j).^2; 	end%% Calculate theta%	theta=ktotheta(k,energy,'xe');%% Calculate the total independent scattering (elastic+inelastic)%	indscat=ww*(coh+inc);%% Fit the total independent scattering to approximate form%	pin=[50 0.1 sumz2 energy];	dp=[0.001 0.001 0 0];	[fit,p]=leasqr(k,indscat,pin,'jfact',0.001,20,1,dp);%% Plot the results for reality check; give warnings or halt as needed.%		plot(k,indscat,k,fit)	title('Total independent scattering and approximation')	xlabel('k');ylabel('Intensity')	disp(sprintf('Approximation constants = %g, %g; area ratio = %g.',p(1),p(2),trapz(k,fit)./trapz(k,indscat)))        if p(2)<0            warning('Compton component from fit (q) negative; setting to zero.')            p(2)=0;    endif lower(geo(1:4))=='refl'    if (p(1)<min(brefl))|(p(1)>max(brefl))            error(['Fit parameter b=' num2str(p(1)) ' outside interpolation table range (' num2str(min(brefl)) '<b<' num2str(max(brefl)) ' ).'])    end    if (p(2)<min(qrefl))|(p(2)>max(qrefl))            error(['Fit parameter q=' num2str(p(2)) ' outside interpolation table range (0' num2str(min(qrefl)) '<q<' num2str(max(qrefl)) ' ).'])    endelseif lower(geo(1:5))=='trans'    if (p(1)<min(btrans))|(p(1)>max(btrans))            error(['Fit parameter b=' num2str(p(1)) ' outside interpolation table range (' num2str(min(btrans)) '<b<' num2str(max(btrans)) ' ).'])    end    if (p(2)<min(btrans))|(p(2)>max(btrans))            error(['Fit parameter q=' num2str(p(2)) ' outside interpolation table range (' num2str(min(qtrans)) '<q<' num2str(max(qtrans)) ' ).'])    endend            disp('Hit any key to continue.')	pause%% Give warnings%	if p(2)<0;p(2)=0;end	if p(1)<10;p(1)=20;end%% Interpolate the appropriate Q table to get approximate value of Q%if lower(geo(1:4))=='refl'		q=interp3(brefl,twothetarefl,qrefl,bigqrefl,p(1),2.*theta,p(2),'cubic');	elseif lower(geo(1:5))=='trans'		mut=mu*sample(1);		q=interpn(twothetatrans,btrans,qtrans,muttrans,bigqtrans,2.*theta,p(1),p(2),mut,'cubic');		q=q(:);%   		plot(k,q,'o');%   		pause	else		error('Unknown scattering geometry')	end%	find(isnan(q))		if ~isempty(find(isnan(q)))		warning('Warning: Some NaN values in i2/i1 - extrapolation range exceeded.')	end		i2i1=sumz2.^2.*q./(indscat.*asum);%--------------------------------------------------------------------------function i2i1=multscat(q,energy,coh,sample,mu,geo);% MULTSCAT%% Usage: i2i1=multscat(q,energy,coh,sample,geo);%% Multiple scattering calculation for a linearly polarized source (e. g.% a synchrotron) with the scattering plane plane perpendicular to the % direction of polarization (i. e. a vertical scattering plane at a synchrotron).%% Input:%         q      = experimental q scale (q=4*pi*sin(theta)/lambda)%         energy = x-ray energy (eV)%         coh    = coherent indenpendent scattering intensity (electron units)%         sample = [thickness density atden] where%                 thickness = sample thickness (cm)%                 density = sample density (g/cm^3)%                 atden = atomic density (atoms/�^3)%         mu     = linear absorption coefficient of sample (cm^-1)%         geo    = scattering geometry:%                  geo(1:4)='refl' for symmetric reflection%                  geo(1:4)='tran' for symmetric transmission%% Output:%         i2i1   = ratio of double- to single-scattering intensity%% References:% 	Malet, et. al., J. Appl. Cryst. 6, 139-44 (1973)%	Derivation adjusted for linear polarization perpendicular to the scattering plane%	from Todd Hufnagel, PhD Thesis, Stanford Univeristy (1995) (Appendix C)% 9/13/2000 Hope Ishii% 	Notes: Currently, the code outputs k vs. i2i1 where the k vector does not% 	correspond to the experimental k vector. Use an interpolation to get i2i1% 	at experimental points. Also, the experimental coherent intensity is not% 	used in the integral for two reasons: 1) the correction is really small in% 	any case (<0.1%) and 2) using the experimental coherent intensity% 	requires that the intensity be normalized to a per electron basis and would% 	make the calculation much longer. Finally, the code is not optimized to take% 	full advantage of Matlab's matrix capabilities - yet - but still runs quickly.% 7/3/02 Todd Hufnagel (hufnagel@jhu.edu)%	Modified to make the experimental q scale an input; the calculation is%	done over a coarser grid, and then interpolated onto the input q scale.%	Also changed from input scattering factors and composition to input%	coherent scattering intensity. (This is calculated anyhow in the data analysis%	process, and doing it this way opens the door to a later modification to use%	actual experimental coherent scattering.) Changed input E to eV%	(from keV), to conform to other code. Permit calculation for either%	reflection or transmission geometry.% 14-May-03 TCH Incorporated into MULTIPLE.M; use INTERP1 in place of% INTERPOL; results in much faster execution.% 15-May-03 TCH Fixed extrapolation bug.warning offticdisp('This may take a few minutes. Why don''t you check your email?')t=sample(1);atden=sample(3);lambda = etolambda(energy,'x'); 			%wavelengthtwoth=[0:5:180]; th=twoth./2; 		% index for theta is ik=4.*pi./lambda.*sind(th);			% k scale (�^-1)coh=spline(q,coh,k);				% coherent intensity per atom as a function of kcohn=coh./coh(1);				% coherent intensity normalized to a per electron basisdphi=2.5;						% step in phi (�)dpsi=2.5;						% step in psi (�)phi=[0:dphi:180];				% index for phi is jpsi=[-90:dpsi:90];				% index for psi is hnelectrons=atden.*sqrt(coh(1)).*(10.^24);	% electrons / cm^3sige=7.377e-26;					% sigma_e = e^4/m^2/c^4  [cm^2]								%     (scattering cross-sxn of an electron)for i=1:length(th)				% cycle through theta	ct=cosd(th(i));				%   & calculate trig functions for current theta	st=sind(th(i));	for j=1:length(phi)			% cycle through phi		cph=cosd(phi(j));		%   & calculate trig functions for current phi		sph=sind(phi(j));		for h=1:length(psi)			% cycle through psi			cps(h)=cosd(psi(h));		%    % calculate trig functions...			sps(h)=sind(psi(h));			if psi(h)==0						% part of absorption correction				A(h)=((ct-sps(h)).^-2).* ...		%    (some terms factored out)					(mu.*t.*(ct-sps(h)));			elseif psi(h)==90-th(i)				% special case psi=(pi/2)-theta				A(h)=mu.*mu.*t.*t./2./(ct.^2);				else				A(h)=((ct-sps(h)).^-2).* ...					(mu.*t.*(ct-sps(h))-ct.*abs(sps(h)).*  ...					(1-exp(-mu.*t.*(ct-sps(h))./(ct.*abs(sps(h))))));			end			P2(h)=1-((cps(h)^2).*(sph.^2)).*  ...	% polarization correction term				(1-((ct.*sps(h)-st.*cps(h).*cph).^2));			thone(h)=acosd(ct.*sps(h)+st.*cps(h).*cph)./2;	% scattering angles for			thtwo(h)=acosd(ct.*sps(h)-st.*cps(h).*cph)./2;	% double scattering			% interpolate to get cohn at thone and thtwo			cohone(h)=interp1(th,real(cohn),real(thone(h)));	% theoretical coherent intensity			cohtwo(h)=interp1(th,real(cohn),real(thtwo(h)));	% for those two scattering events			% full terms to integrate - I2 is the special case where psi~(pi/2)-theta			I1(h)=A(h).*P2(h).*cohone(h).*cohtwo(h).*cps(h);			I2(h)=(mu.*mu.*t.*t./2./ct).*P2(h).*cohone(h).*cohtwo(h).*cps(h);		end;		% carry out integral over psi (there are three ranges of integration)		if i==1			jM(j)=sum(dpsi.*pi./180.*I1(1:length(psi)-i))+dpsi.*pi./180.*I2(length(psi));		elseif i==2			jM(j)=sum(dpsi.*pi./180.*I1(1:length(psi)-i))+dpsi.*pi./180.*I2(length(psi)-i+1)+ ...				dpsi.*pi./180.*I1(length(psi));		elseif i==length(th)-1			jM(j)=dpsi.*pi./180.*I1(1)+dpsi.*pi./180.*I2(2)+sum(dpsi.*pi./180.*I1(3:length(psi)));		elseif i==length(th)			jM(j)=dpsi.*pi./180.*I2(1)+sum(dpsi.*pi./180.*I1(2:length(psi)));		else			jM(j)=sum(dpsi.*pi./180.*I1(1:length(psi)-i))+dpsi.*pi./180.*I2(length(psi)-i+1)+ ...				sum(dpsi.*pi./180.*I1(length(psi)-i+2:length(psi)));		end;	end;	% calculate ratio of double to single scattering	iM(i)=2.*sige.*nelectrons.*ct./mu./mu./t./cohn(i).*sum(dphi.*pi./180.*jM);end;warning onbeepbeephowlong=toc;minutes=floor(howlong/60);seconds=round(rem(howlong,60));if minutes>0    disp(['That only took ' num2str(minutes,'%0.5g') ' minutes, ' num2str(seconds,'%0.5g') ' seconds.'])else    disp(['That only took ' num2str(seconds,'%0.5g') ' seconds.'])endplot(k,iM);disp('Hit any key to continue');pause% Interpolate onto input q scale.% Last point of iM is INF, so discard it and do an extrapolation.lengthk=length(k);k(1:lengthk);iM(1:lengthk);i2i1=interp1(k(1:lengthk-1),iM(1:lengthk-1),q,'spline');%-------------------------------------------------------------------------% Here follow some useful functionsfunction y=cosd(x)y=cos(x.*pi./180);function y=sind(x)y=sin(x.*pi./180);function y=acosd(x)y=180./pi.*acos(x);function y=asind(x)y=180./pi.*asin(x);