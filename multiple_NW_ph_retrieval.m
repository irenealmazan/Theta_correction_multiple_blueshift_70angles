% this script shows the results

addpath(genpath('./m_scripts/'));
addpath(genpath('./calc_functions'));


jitterlevel = [0 5 10 20 40];
mncrate_array = [2e2 2e2 1e3];
noiseflag_array = [0 1 1];
noiselevel_array = [1 2 3];
%%%% no noise



for mm = 1:numel(mncrate_array)
    
    noiseflag = noiseflag_array(mm);
   
    mncntrate = mncrate_array(mm);
    
    noiselevel_str = num2str(noiselevel_array(mm));
    
    for jjj = 1:numel(jitterlevel)
        
       percent = jitterlevel(jjj);
       
        savefolder = ['jitter_' num2str(jitterlevel(jjj)) '_noiselevel_' noiselevel_str '_70angles'];
        mkdir(savefolder);
       
       load(['data_initial_alljitter/data_cnt_' noiselevel_str '_jitter_' num2str(percent)]);
   
       load(['data_ERHIO/struct_ERHIO_ini' noiselevel_str '_jitter_' num2str(percent)]);

       support_new = struct_best_ERHIO.support_new;
       newobj.dp = struct_best_ERHIO.dp;
       
        %{
        figure(1);
        clf;
        for kk = 1:numel(data_exp)
            imagesc(data_exp(kk).I);
            colorbar;
            axis image;
            title([' noiselevel ' noiselevel_str ' noiseflag= ' num2str(noiseflag) 'jitter =' num2str(percent) 'mncnrate = ' num2str(mncntrate) 'kk = ' num2str(kk)]);
            pause(.1);
        end
        %}
        
        NW_ph_retrieval_BCDI;
        
       save([savefolder '/results.mat']);

    end
    
end