addpath(genpath('./m_scripts/'));
addpath(genpath('./calc_functions'));


jitterlevel = [0 5 10 20 40];%0 10 20 40];
mncrate_array = [2e2 2e2 1e3];
noiseflag_array = [0 1 1];
noiselevel_array = [1 2 3]
%%%% no noise



for mm = 1:numel(mncrate_array)
    
    noiseflag = noiseflag_array(mm);
    if(noiseflag)
        display('ADDING NOISE')
    else
        display('NO NOISE')
    end
    
    mncntrate = mncrate_array(mm);
    
    noiselevel_str = num2str(noiselevel_array(mm));
    
    for jjj = 1:numel(jitterlevel)
        
        percent = jitterlevel(jjj);
       
        load(['data_initial_alljitter/data_cnt_' noiselevel_str '_jitter_' num2str(percent)]);

        min_chi = 1;
        for kk = 1:5
            ER_HIO;
           struct_ER_HIO(kk).chi = newobj.chi;     
           struct_ER_HIO(kk).dp = newobj.dp;
            struct_ER_HIO(kk).support_new = support_new;
            
           if min_chi <newobj.chi(end)
               min_kk = kk;
               min_chi = newobj.chi(end);
           end
           
        end
        
        struct_best_ERHIO = struct_ERHIO(min_kk);
        
        save(['data_ERHIO/struct_ERHIO_ini' noiselevel_str '_jitter_' num2str(percent)],'struct_best_ERHIO');
        
        
    end
    
end