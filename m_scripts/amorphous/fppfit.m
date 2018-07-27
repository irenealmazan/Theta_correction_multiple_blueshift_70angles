function [fpp,den]=fppfit(energy,xas,atom,edgename,erange,nfit);% FPPFIT - Determine f'' from absorption data %% Usage: [fpp,den]=fppfit(energy,xas,atom,{edgename,erange,nfit})%% Determines the anomalous scattering factor f'' from measured% absorption data for a single element. The absorption from% all other elements and the detector response function are% removed, and the absorption is normalized to a parameterized% f'' (e.g. Cromer-Liberman) far from the edge.%% Input:%   erg       = energy array for absorption data (eV)%   xas       = array of absorption data (ln(Io/I))%   atom      = atomic name (string) or number for element of interest%   edgename  = edge of interest. Can be 'K', 'LI', 'LII', or 'LIII' %               (optional, defaults to 'K').%   erange    = [emin emax] lower and upper bounds on %               range of energies near edge excluded from fit%               (optional; defaults to 50 eV below and 100 eV above)%   nfit      = integer for maximum order of fit (can be 1-6)%               (optional; default is 2)%% Output:%   fpp     = experimental f'' for element of interest%   den     = number density * thickness for this element (1/AA^2)%% See also: FPPFITG, KRAMKRON% Reference: Lane Wilson Ph.D. thesis, Stanford University, 1990.% Based on FPPFIT.FOR written by Lane Wilson 3/3/86% originally based on FPPCL.FOR written by Karl Ludwig% Matlab version 7/21/2000 Hope Ishii - hope.ishii@stanford.edu% Batch mode 23-Jul-00 SMB Bren@slac.stanford.edu% 16-May-03 TCH Error checking, better comments, and code clean-up%           Do any edge (not only K).% Requires: element.m, atomdata.m, anomal.m% Use column vectorsenergy=energy(:);xas=xas(:);% Verify size of energy and XAS arrays are the sameif length(energy)~=length(xas),        error('Energy and XAS data are not of same length')end% Verify Z_edge is an atomic numberatom=element(atom);if nargin <6,        nmax=2;else        if nfit>1 & nfit < 6,                nmax= nfit;        else                nmax=2;        endend% Obtain edge energies for element of interest[dum1,dum2, alledge] = atomdata(atom);% Pick out the edge that we wantif nargin<4    edge=alledge(1);               % Default to K edgeelseif isempty(edgename)    edge=alledge(1);else        edgename=lower(edgename);    switch edgename        case 'k'                    % K edge            edge=alledge(1);        case 'li'            edge=alledge(2);        case 'lii'            edge=alledge(3);        case 'liii'            edge=alledge(4);        otherwise            error('Unknown absorption edge name (must be K, LI, LII, or LIII)')    endendedge=edge.*1000;                    % Convert from keV to eV% Error checking on edgeif (edge<energy(1)|edge>energy(end))    error(['Edge at E=' num2str(edge) ' eV not within given range of energies. Do you have the right element and edge?'])end% Determine boundaries of excluded range of energies near edge% Do Emin firstif nargin < 5,   Emin= edge-50;                   % Emin defaults to 50 eV below edge   if Emin< energy(1),              % unless that exceeds range of data        Emin= mean([energy(1) edge]);   endelse                                % Sanity check on user-specified Emin   Emin= erange(1);   if Emin>edge      error('Specified emin for excluded range is above edge!')   end    end% Now do Emaxif nargin < 5 | length(erange)==1,  % Emax defaults to 100 eV below   Emax= energy(end)- 100;          % maximum in data range   if Emax < edge,                  % unless that is below the edge        Emax= mean([edge energy(end)]);    endelse        Emax= erange(2);        if Emax<edge          error('Specified emax for excluded range is below edge!')        end    end% Calculate the parameterized f', f'' (not necessarily Cromer-Liberman;% depends on setting of global ANOMAL_FLAG. See ANOMAL for details.[CLfp,CLfpp]=anomal(atom,energy);% d1 and d2 are indices to energies other than those in the excluded% region (below and above edge, respectively)d1=find(energy<Emin);d2=find(energy>Emax); erglim=[energy(d1); energy(d2)];% Get parameterized f'' values for energies outside edge regionCLfpplim=[CLfpp(d1); CLfpp(d2)];% Set up matrix for solution to% ln(Io/I)= constant*CLfpp/E + B(E)%     where B(E)= C1+C2/E+C3/E^2+C4/E^3+ ...V=zeros(length(erglim), nmax+2);V(:,1)=1;V(:,2)=1./erglim;for j=3:nmax+1, V(:,j)=V(:,j-1).*V(:,2); endV(:,nmax+2)=CLfpplim.*V(:,2);% Solution:   ln(Io/I) for the energies outside the edge regionW=[xas(d1); xas(d2)];% Solve linear set of equations for the coefficients%C=V\W(:);C=pinv(V)*W;B=C(1)*ones(size(energy));% B(E)=backgroundfor j=2:nmax+1    B=B+C(j).*energy.^(1-j);endyfit=B+C(nmax+2).*CLfpp./energy;    % constant*CLfpp/E + B(E)% Calculate (number density * thickness)den=C(nmax+2)./0.69886;      % 0.69886=(2he^2)/(mc) [eV*AA^2]%disp(sprintf(' Number density * thickness =  %6.2f [atoms/AA^2]', den));% Experimental XAS - background function, B(E)xasnobg=xas-B;% Experimental XAS without B(E) and normalized to Cromer-Liberman f''% values far from the edgefpp=xasnobg.*energy./C(nmax+2);% Visual check of resultsfigure(1);plot(energy,xas,'r',energy, yfit,'b',energy,B,'m',energy(d1),xas(d1),'g',energy(d2),xas(d2),'g')title('Raw data, fit and background function')xlabel('Energy (ev)'); ylabel('ln(Io/I)');legend('Raw data', 'Fit to data','Background','Fitted regions')figure(2), plot(energy, fpp, energy, CLfpp)title(['Parameterized f" and the experimental f"'])%(nt)_{\alpha} = ', ...num2str(den,5),' atoms/area [�^{-2}]']);xlabel('Energy (eV)'), ylabel('f" (electrons)');legend('Experimental f"','Parameterized f"')