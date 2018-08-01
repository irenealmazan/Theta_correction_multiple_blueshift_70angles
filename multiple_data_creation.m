addpath(genpath('../m_scripts/'));
addpath(genpath('../calc_functions'));


jitterlevel = [70 90];%[0 5 10 20 40];%0 10 20 40];
mncrate_array = [2e2 2e2 1e3];
noiseflag_array = [0 1 1];
noiselevel_array = [1 2 3];%[1 2 3];
%%%% no noise

NW_flags;

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
       
        NW_make_data_exp;
        
                
        save(['data_initial_alljitter/data_cnt_' noiselevel_str '_jitter_' num2str(percent)]);
        
        
    end
    
end