classdef InitializeFunctions
    % This library contains all the functions initializing the parameters
    % to be used in the different scripts
    properties(Constant)
    end
    
    
    methods(Static)
        
        function [Niter_rho, Niter_pos,...
                Niter_theta,freq_pos,freq_store,freq_restart,freq_shrink_wrap,tau_backtrack_rho,beta_ini_rho,...
                counter_max_rho,tau_backtrack_theta,beta_ini_theta,counter_max_theta] = NW_experimental_phretrieval_parameters()
            
            % In this script we initialize the values of the experimental set-up and
            % the details concerning the sample:
          
            
            %% Phase retrieval parameters:
            
            % Iteration parameters:
            Niter_rho = 2000;
            Niter_pos = 1;
            Niter_theta = 1;
            freq_pos = 1;
            %freq_rho = 10;
            freq_store = 100;
            freq_shrink_wrap = 500;
            
            % Beta adaptative step parameters and conjugated gradient restart parameter:
            freq_restart = 20;
            
            tau_backtrack_rho = 0.1;
            beta_ini_rho = 0.5;
            counter_max_rho = 10;
            
            tau_backtrack_theta = 0.5;
            beta_ini_theta = 0.1;
            counter_max_theta = 15;
            
           
        end
        
        
        function [pixsize,lam,Npix,detdist,d2_bragg,depth,defocus,th,del,gam,...
                thscanvals,alphavals,phivals,...
                delta_thscanvals] = NW_scatgeo_2110()
            
            % This function contains the value of the experimental parameterse
            % when the [2110] reflection was measured.
            
            % experimental setup details
            pixsize = 55; %microns, for Merlin
            lam = etolambda(10400)*1e-4;
            
            Npix = 128;
            detdist = 0.529e6; % in micrometers
            d2_bragg = detdist * lam /(Npix*pixsize);
            depth = 70;%128;%60;
            defocus = 0;
            
            % expected values of the motors for m-plane
            th = 73.3;
            del = -32.6; %in plane
            gam = 0; %out of plane
            
            
            % rocking curve scans
            if mod(depth,2)
                thscanvals =  [72.8:1.0/depth:73.3-1.0/depth 73.3 73.3+1.0/depth:1.0/depth:73.8];
            else
                thscanvals =  [72.8:1.0/depth:73.8-1.0/depth];
            end
            
            alphavals = zeros(numel(thscanvals),1);
            phivals = zeros(numel(thscanvals),1);
            delta_thscanvals = thscanvals-th;%[thscanvals-th;alphavals';phivals']';
            
            %usesimI = 1;
            
            
        end
        
        function [pixsize,lam,Npix,detdist,d2_bragg,depth,defocus,th,del,gam,...
                thscanvals,alphavals,phivals,...
                delta_thscanvals] = NW_scatgeo_1010()
            
            
            % This script contains the value of the experimental parameterse
            % when the [1010] reflection was measured.
            
            % experimental setup details
            pixsize = 55; %microns, for Merlin
            lam = etolambda(10400)*1e-4;
            
            Npix = 516;
            detdist = 0.33e6; % in micrometers
            d2_bragg = detdist * lam /(Npix*pixsize);
            depth = 100;
            defocus = 0;
            
            % expected values of the motors for m-plane
            th = -9.3;
            del = -18.6; %in plane
            gam = 0; %out of plane
            
            % rocking curve scans
            thscanvals =  [-9.42:.02:-9.13];
            alphavals = zeros(numel(thscanvals),1);
            phivals = zeros(numel(thscanvals),1);
            delta_thscanvals = thscanvals-th;%[thscanvals-th;alphavals';phivals']';
           
          
        end
        
    end   
        
end
    

