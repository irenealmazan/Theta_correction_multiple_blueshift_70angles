function [retrphase newobj] = hio3(dp, support, num, savenum, prevobj, beta,plotResults)

dp = abs(dp);

%initial guess as to the phase of the DP comes from a FT of the support
if isempty(prevobj) 
    %A = fftn(support); 
    A = ones(size(support)) .* exp(sqrt(-1) * rand(size(support))*2*pi - pi );
end

if ~isempty(prevobj) A = fftn(prevobj.object); end
%imagesc(fftshift(log(abs((A)))));pause%
%imagesc(fftshift(angle((A))));pause

if plotResults
    h=waitbarpos(0, 'hio phase iterations');
    
    if isempty( findobj('Name', 'phasing'))
        hfig = figure;
        set(gcf, 'Name','phasing');
    else
        hfig = findobj('Name', 'phasing');
    end
    
    if isempty( findobj('Name', 'chi'))
        hchi = figure;
        set(gcf, 'Name','chi');
    else
        hchi = findobj('Name', 'chi');
    end
    
    figure(hchi);clf;
    
end

ssqexpamp = sum(sum(sum(abs(dp))));

if isempty(prevobj) 
    Dold= abs(support); 
    %Dold = 
end
if ~isempty(prevobj) Dold = prevobj.object; end

outsup = find(abs(support) == 0);
insup = find(abs(support) ~= 0);

chi=[];
chilen=0;
if ~isempty(prevobj) 
    chi = prevobj.chi; 
    chilen = length(chi);
end

ii=sqrt(-1);

for i=1:num
    
    
    %we apply the phase guess to the measured dp amplitudes
    B = abs(dp) .* exp(ii * angle(A));
    
    %we bring the new guess to real space via IFT
    C = ifftn(B);

    %we knock off all intensity in real space that falls outside the bounds
    %of the support
    
    if(i>1) Dold = D; end
    
    D=C;
    %D(insup) = C(insup);
    D(outsup) = Dold(outsup) - beta * C(outsup);
        
    if(mod(i,savenum)==0)        
        displayphasing(A, B, C, D, support, 3, hfig, hchi);
    end
    
    %return to dp space and apply the 
    A = fftn(D);

    %find a chi fit.
    chi(chilen+i,1) = sum(sum(sum( (abs(A)-abs(dp)).^2))) / numel(dp) ;%sum(sum(sum( (abs(A)-abs(dp)).^2))) / ssqexpamp ;
    
    if plotResults
        figure(hchi);
        plot([1:length(chi)], log10(chi));
        title(['chi fit is ' num2str(chi(i,1))]);
        ylabel('log chi');
        
        waitbar(i/num,h);
    end
end

if plotResults
    close(h);
end
%imagesc(angle(A));

retrphase = abs(dp).*exp(sqrt(-1)*angle(A));

%find the minimum number of points to describe the support
[ind1 ind2 ind3] = ind2sub(size(support), find(support==1));
suphull = [ind2 ind1 ind3];
K = convhulln(suphull);
suphull = suphull(unique(K),:);

if isempty(prevobj)
    histnum=1;
    newobj.history(histnum).alg = 'hio';
    newobj.history(histnum).sup = inputname(2);
    newobj.history(histnum).suphull = suphull;
    newobj.history(histnum).dp = inputname(1);
    newobj.history(histnum).beta = beta;
    newobj.history(histnum).iter = i;
    newobj.chi = chi;
end

if ~isempty(prevobj)
    
    newobj = prevobj;
    histnum = length(prevobj.history);
    histnum= histnum+1;
    newobj.history(histnum).alg = 'hio';
    newobj.history(histnum).sup = inputname(2);
    newobj.history(histnum).suphull = suphull;
    newobj.history(histnum).dp = inputname(1);
    newobj.history(histnum).beta = beta;
    newobj.history(histnum).iter = i;
    newobj.chi = chi; 
end

if num>0
    newobj.object = D;
    newobj.dp = fftshift(A);
    totiter = 0;
    for i=1:histnum totiter = totiter + newobj.history(i).iter; end
    newobj.totiter = totiter;
end

if num ==0 newobj = prevobj; end
