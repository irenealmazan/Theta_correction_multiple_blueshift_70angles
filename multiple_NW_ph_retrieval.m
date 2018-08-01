% this script shows the results

addpath(genpath('../m_scripts/'));
addpath(genpath('../calc_functions'));


jitterlevel_1 = [0 5 10 20 40];%0 10 20 40];
mncrate_array_1 = [1e3];%[2e2 2e2 1e3];
noiseflag_array_1 = [1];%[0 1 1];
noiselevel_array_1 = [3];%[1 2 3];
%%%% no noise



for mm = 1:numel(mncrate_array_1)
    
    noiseflag = noiseflag_array_1(mm);
   
    mncntrate = mncrate_array_1(mm);
    
    noiselevel_str = num2str(noiselevel_array_1(mm));
    
    for jjj = 1:numel(jitterlevel_1)
        
       percent = jitterlevel_1(jjj);
       
        savefolder = ['jitter_' num2str(jitterlevel_1(jjj)) '_noiselevel_' noiselevel_str '_70angles'];
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