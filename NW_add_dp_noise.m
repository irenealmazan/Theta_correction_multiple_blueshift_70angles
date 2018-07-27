% In this script we add Poison noise to the simulated data and we set the
% signal to background ratio with mncntrate; we also overwrite the
% experimental data with the simulated data.

mn = mean(data_exp(middpind).simI(:));
noise_data = zeros(size(X));
clean_data = zeros(size(X));
only_noise= zeros(size(X));

for ii=1:numel(data_exp)

    if (noiseflag)
        if(mod(ii,50)==0) display(['adding noise to sim dp ' num2str(ii) ' of ' num2str(numel(data_exp))]); end
        data_exp(ii).noiseI = poisrnd( data_exp(ii).simI /mn * mncntrate);
    else
        data_exp(ii).noiseI = ( data_exp(ii).simI /mn * mncntrate);
    end

    clean_data(:,:,ii) = data_exp(ii).simI /mn * mncntrate;
    noise_data(:,:,ii) = data_exp(ii).noiseI;
    only_noise(:,:,ii) = noise_data(:,:,ii) - clean_data(:,:,ii);
    
end

% signal_to_noise ratio:
snr = sum(sum(sum(clean_data.^2)))/sum(sum(sum(only_noise.^2)));

phase_NW = angle(fftshift(fftn(fftshift(NW))));
noise_NW = fftshift(ifftn(fftshift(sqrt(noise_data).*exp(1i*phase_NW))));

if plotResults
fig_num = 401;
dimension = '2';
DisplayResults.compare_two_objects(NW,noise_NW,'Original object','Compatible object',[50 80 1 70],64,dimension,fig_num);
end
    


%%

%%
%if use simulated data, overwrite the experimental I with the noisy simI
%and calculates the noisy rocking curve

rock_curve_noise = zeros(numel(data_exp),1);

if(usesimI)
    display(['overwriting experimental I with sim I']);
    for ii=1:numel(data_exp)
        if (noiseflag)
            data_exp(ii).I = data_exp(ii).noiseI;            
        else
            data_exp(ii).I = data_exp(ii).simI /mn * mncntrate;            
        end
        rock_curve_noise(ii) = sum(sum(data_exp(ii).I));
        
    end
end

if plotResults
    figure(7);
    clf; setfigsize(gcf, 800,400);
    colormap jetvar;
    
    ca = [0 2];
    for ii=1:numel(data_exp)
        subplot(121);
        imagesc(data_exp(ii).I);axis image;colorbar;
        
        drawnow;
        subplot(122);
        imagesc( (data_exp(ii).simI)); axis image;colorbar;
        title(num2str(ii));
        
        display(['dth: ' num2str(data_exp(ii).dth_nominal)]);
        pause(.5);
    end
end

    
% plot rocking curve: 
if plotResults
    if newSample
        DisplayResults.show_rocking_curve(delta_thscanvals,rock_curve_3D/mn * mncntrate,'new',8,'k','true using 3D FT');
        DisplayResults.show_rocking_curve(delta_thscanvals+dth_disp',rock_curve_noise,'hold',8,'r','true using 2D FT');
    else
        DisplayResults.show_rocking_curve(delta_thscanvals+dth_disp',rock_curve_noise,'new',8,'r','true using 2D FT');
    end
end


