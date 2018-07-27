% jitterlevel = [10 20 40];%0 10 20 40];
% 
% 
% for jjj = 1:numel(jitterlevel)
%     
%     percent = jitterlevel(jjj);
%       
%      mncntrate = 1e3;
%      
%      noiseflag = 0;
%      if(noiseflag)
%          display('ADDING NOISE')
%      else
%          display('NO NOISE')
%      end
%      
%     savefolder = ['jitter_' num2str(jitterlevel(jjj)) '_noiselevel_0'];
%     mkdir(savefolder);     
%     
%     NW_masterscript_BCDI;
%     
%     save([savefolder '/results.mat']);
%     
% end

%%%% noise
jitterlevel = [10 20 40];%
for jjj = 1:numel(jitterlevel)
    
    percent = jitterlevel(jjj);
    display(['jitterlevel = ' num2str(jitterlevel(jjj))]); 
    
    
     mncntrate = 1e3;
     
     noiseflag = 1;
     if(noiseflag)
         display('ADDING NOISE')
     else
         display('NO NOISE')
     end
     
    savefolder = ['jitter_' num2str(jitterlevel(jjj)) '_noiselevel_1'];
    mkdir(savefolder);     
    
    NW_masterscript_BCDI;
    
    save([savefolder '/results.mat']);
    
end

