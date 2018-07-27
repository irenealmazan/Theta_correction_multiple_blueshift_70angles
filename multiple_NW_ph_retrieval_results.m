% this script shows the results

addpath(genpath('./m_scripts/'));
addpath(genpath('./calc_functions'));


jitterlevel_1 = [0 5 10 20 40];%0 10 20 40];
mncrate_array_1 = [2e2 2e2 1e3];
noiseflag_array_1 = [0 1 1];
noiselevel_array_1 = [1 2 3];

% flip is required for jitter_0_noiselevel_0, jitter_5_noiselevel_0, 
flipflag_array = [1 1 0 1 0];

%%%% no noise



for mm = 1:numel(mncrate_array_1)
    
    noiseflag = noiseflag_array_1(mm);
   
    mncntrate = mncrate_array_1(mm);
    
    noiselevel_str = num2str(noiselevel_array_1(mm));
    
    for jjj = 1:numel(jitterlevel_1)
        
       percent = jitterlevel_1(jjj);
       
       flipflag = flipflag_array(jjj);
       
        savefolder = ['jitter_' num2str(jitterlevel_1(jjj)) '_noiselevel_' noiselevel_str '_70angles'];
        load([savefolder '/results.mat']);

        load(['data_ERHIO/struct_ERHIO_ini' noiselevel_str '_jitter_' num2str(percent)]);

        
        newobj.chi = struct_best_ERHIO.chi;
        
        NW_figures_single_jitter;
        

    end
    
end