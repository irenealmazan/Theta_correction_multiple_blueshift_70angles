function [x,y]= loadspec(name,scan)% LOADSPEC - Load data from a SPEC file (or a SPEC-compatible SUPER file)%            and extract one scan   %% Usage: [x,y]= loadspec(name{,scan})%% Input: name = file name (string)%        scan = scan number (scalar, optional for SUPER files)%% Output:   x = sqrt(h^2+k^2+l^2) for hkl scans; for all other types of scans,
%               it is whatever is in the first column of SPEC data (typically a%               motor position.%%           y = intensity data (presumed to be last column) for SPEC files%             = matrix of all columns (other than x) for SUPER files% 05-Feb-96 Sean Brennan SSRL Bren@slac.stanford.edu% 08-Jan-99 Todd Hufnagel JHU hufnagel@jhu.edu%  Changes:%   (1) Fixed bug that truncated first character of data array%   (2) Fixed bug that would result in wrong scan loading if (say)%		scans numbered 1 and 10 existed in the same file%	(3) Fixed bug that prevented extraction of last scan in file%	(4) For hklscans, made x=sqrt(h^2+k^2+l^2)%   (5) Cleaned up help comments and added error checking% 08-Mar-02 Todd Hufnagel JHU hufnagel@jhu.edu%	(1) Added checking for scan numbers padded with zeros%	(2) If SUPER file, return all data as y matrix, not just last column% 24-Apr-04 Todd Hufnagel JHU hufnagel@jhu.edu%   (1) Added explicit SUPER file check at beginning%   (2) Fixed bug that caused error in reading SUPER files numbered >=100%   (3) Verified to work with both SSRL and JHU SPEC files%   (4) Made specification of scan number optional for SUPER files% Notes: Entire file is read in at the beginning and the% parsed internally.  Much faster for reasonably small files,% but memory hungry for larger files.fid=fopen(name,'r');[dch,count]=fread(fid,inf,'uchar');fclose(fid);chdat=[dch(:)',setstr(13)];% replace lf's with spaces or cr's, depending on platformid10= find(chdat==setstr(10));comp= computer;if strcmp(comp(1:3),'PCW')|strcmp(comp(1:3),'VAX')|strcmp(comp(1:3),'ALP'),	chdat(id10)= setstr(' '*ones(1,length(id10)));else	chdat(id10)= setstr(13*ones(1,length(id10)));end% figure out whether this is a SUPER file in disguise% SUPER files will always have a comment line with the ring currentsuperfile=0;searchstring=['#C Current='];if ~isempty(findstr(chdat,searchstring))    superfile=1;    if nargin<2 % if not given, determine the scan number        scan=str2num(name((length(name)-2):length(name)));    else        % if given, check for consistency        if scan~=str2num(name((length(name)-2):length(name)))            error('Scan number does not match SUPER file name.')        end        end    end% For SPEC file, check to ensure that scan number is specifiedif (~superfile)&(nargin<2)    error('Scan number must be specified for SPEC files.')end    searchstring=['#S '];% idscns is the indexes to all the scansidscns= findstr(chdat,searchstring);idscns=[idscns count];	% last char is boundary% locate the particular scan that we want. The search string is padded with% spaces to ensure that we don't locate (say) scan 10 when we are looking for% scan 1searchstring= ['#S ',num2str(scan),'  '];idscan= findstr(chdat,searchstring);% If not found on first try, maybe it's because the scan number is padded% with zeros (as in Sean's SUPER files with scan numbers <100)if (length(idscan)==0)	searchstring= ['#S ',num2str(scan,'%1.3d')];	idscan= findstr(chdat,searchstring);end	% Do some error checking on the number of particular scans foundif (length(idscan)==0)	error(['Scan ' num2str(scan) ' not found in file ' name])elseif (length(idscan)>1)	warning([num2str(length(idscan)) ' instances of scan ' num2str(scan) ' found. Loading last instance.'])end% indscan is the index of the last instance of the specific scanindscan= idscan(length(idscan));% next is the index of the start of the next scannext=1+ find(idscns==indscan);% determine whether this is an hkl scanhklscan=strcmp(chdat(indscan+7:indscan+13),'hklscan');% find instance of #N between two instances of #Sindn=findstr(chdat(indscan:idscns(next)),'#N');indn=indscan+indn;% find all of the cr's in the scan of interestindcr=findstr(chdat(indn:idscns(next)),setstr(13));indcr=indcr+indn;% extract the number of columns in the datacolumns= str2num(chdat(indn+2:indcr(1)-1));% find instances of #. If last scan, there is no second # so make second% entry the EOF. If it's a SUPER file, #L defines the last row before the data,% and there is by definition only one scan per file.if superfile	% find instance of #L between two instances of #S	indl=findstr(chdat(indscan:idscns(next)),'#L');	indl=indscan+indl;	indcr=findstr(chdat(indl:idscns(next)),setstr(13));	indcr=indcr+indl;	indl(2)=length(chdat)-1;	chdat(indcr(1):indl(2)-1);	data=sscanf(chdat(indcr(1):indl(2)-1),'%g',[columns inf])';else		indp=findstr(chdat(indn:idscns(next)),'#');	indp=indp+indn;	if (length(indp)==1)		indp(2)=length(chdat)+1;	end	data=sscanf(chdat(indcr(2):indp(2)-1),'%g',[columns inf])';end% assign abscissa (depends on type of scan)if hklscan	x=sqrt(data(:,1).^2+data(:,2).^2+data(:,3).^2);else	x=data(:,1);end% assign ordinate (depends on file type)if superfile	if hklscan		y=data(:,4:columns);	else		y=data(:,2:columns);	endelse % (assumed to be in last column)	y=data(:,columns);end