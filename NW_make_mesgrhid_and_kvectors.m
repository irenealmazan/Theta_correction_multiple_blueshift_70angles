% this script prepares the scattering geometry

[pixsize,lam,Npix,detdist,d2_bragg,depth,defocus,th,del,gam,...
    thscanvals,alphavals,phivals,...
    delta_thscanvals] = InitializeFunctions.NW_scatgeo_2110();

% scattering geometry
NW_diff_vectors_BCDI; % does both the vectors ki and kf and creates the object

% load sample
load('Original_sample_70angles.mat');

if(addNWstrain)
    NW  = img;
else
    NW  = abs(img);
end

if plotResults
    NW_plot_diff_vectors_sample_BCDI;
end

probe = ones(size(X));