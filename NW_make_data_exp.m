% This script prepares the structure data_exp



    
[Niter_rho, Niter_pos,...
    Niter_theta,freq_pos,freq_store,freq_restart,freq_shrink_wrap,tau_backtrack_rho,beta_ini_rho,...
    counter_max_rho,tau_backtrack_theta,beta_ini_theta,counter_max_theta] = ...
    InitializeFunctions.NW_experimental_phretrieval_parameters();


%% Scattering condition:
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

%% Calculate diffraction patterns
NW_calc_dp_BCDI;
NW_add_dp_noise;