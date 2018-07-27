function w = Ludwigmunrow(xa,xb,fa0,fa1,fa2,fb0,fb1,fb2)% w = Ludwigmunrow(xa,xb,fa0,fa1,fa2,fb0,fb1,fb2)%  Calculates the (wrong) Munro matrix, w, for the equation I=wS in which we seek%  the partial structure factors, Sij for a binary alloy. (Differences in%  I must be calculated the same direction as the differences in fa1, etc!)%% Input:%	xa,xb	= mole fractions of elements a and b in sample%	fa0,fb0 = column vectors of the total scattering factor for elements%				a and b for an independent measurement (size mx1 where the%				m direction contains f as a function of k)%	fa1,fb1	= matrices of total scattering factors for elements a and b%				below the a absorption edge for different energies (size%				mxn where n is the # of energies and m is the length of k)%	fa2,fb2	= mxn matrices (or mx1 vectors) of total scattering factors%				below b edge (or a third independent measurement)% Output:%	w		= matrix of weighting factors (M in I=MS) -- typically a 3x3%			   matrix for 2 intensity differences and one independent%			   measurement, but >1 difference can be used below each edge.%% Originally by ?? -- LOTS of mistakes fixed Jan 2001% by Hope 'I'm-SO-sick-of-bugs' Ishii% The 1st column of I (actually I-<f^2>) contains a total intensity curve% The corresponding 1st row of the Munro matrix is...w(1,1) = xa*fa0.*conj(fa0);w(1,2) = 2*xa*real(fa0.*conj(fb0));w(1,3) = xb*fb0.*conj(fb0);% The columns 2 ... n1 of I contain differences at edge a:%	I(E1)-<f(E1)^2>-I(E2)+<f(E2)^2> where E1 and E2 are below edge a% The corresponding 2nd row of the Munro matrix is...[m1,n1]=size(fa1);	j = 1;for i = 2:n1,   w(i,1) = xa*(fa1(j).*conj(fa1(j))-fa1(i).*conj(fa1(i)));   % Note that we're subtracting |f(E1)|^2-|f(E2)|^2 - Make sure E1 and E2   % have the same definitions in calculating the intensity difference!   w(i,2) = xa*(real(fa1(j).*conj(fb1(j)))-real(fa1(i).*conj(fb1(j))));   w(i,3) = 0;end;% The columns n1+1 ... of I contain data at edge b:%	I(E3)-<f(E3)^2>-I(E4)+<f(E4)^2> where E3 and E4 are below edge b% or another total intensity curve (if n2 = 1):  I-<f^2>[d1,d2]=size(w);	[m2,n2]=size(fa2);if n2 == 1,   w(d1+1,1) = xa*fa2.*conj(fa2);   w(d1+1,2) = 2*xa*real(fa2.*conj(fb2));   w(d1+1,3) = xb*fb2.*conj(fb2);end;if (n2 >= 2),	j = 1;   for i = 2:n2,      w(i+d1-1,1) = 0;      w(i+d1-1,2) = xb*(real(fa2(j).*conj(fb2(j)))-real(fa2(j).* ...	  		conj(fb2(i))));      w(i+d1-1,3) = xb*(fb2(j).*conj(fb2(j))-fb2(i).*conj(fb2(i)));   end;end;w=real(w);