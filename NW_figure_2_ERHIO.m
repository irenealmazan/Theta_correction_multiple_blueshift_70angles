% this script the calculations for ER/HIO simulations are runned with the
% different degrees of jittering employed 

flagER_direct = 0;

jitterlevel_summary = [0 5 10 20 40];
noiselevel_str = '1';
for kk = 1:numel(jitterlevel_summary)
    
    savefolder = ['jitter_' num2str(jitterlevel_summary(kk)) '_noiselevel_' noiselevel_str '_freqsw500'];
    load([savefolder '/results.mat']);
    
    ER_HIO;
    
    er_iter = 2000;
    flagER = 1;
    [retrphase,newobj,err_ERHIO] = Phretrieval_functions.do_ERHIO(err_ERHIO,dp,support_new,newobj,er_iter,original_object,delta_thscanvals,ki_o,kf_o,probe,d2_bragg,X,Y,Z,plotResults,flagER,flagER_direct);
    
    
    ER_HIOstruct(kk).chi = newobj.chi ;
    ER_HIOstruct(kk).dp = newobj.dp;
    ER_HIOstruct(kk).support_new = support_new;
    
    
    
    %%%%%%% prepare figures to plot
    
     flipflag = 0;%flipflag_list(kk);
    
    if flipflag
        rho_plot = ifftn(conj(newobj.dp));
        support_plot = abs(ifftn(conj(fftn( ER_HIOstruct(kk).support_new))));
    else
        rho_plot = ifftn(newobj.dp);
        support_plot = struct_err(kk).support_new;
    end
    
    
    rho_shift = DiffractionPatterns.shift_object(NW*sqrt(mncntrate/mn),rho_plot,delta_thscanvals,ki_o,kf_o,kf_o-ki_o,d2_bragg,X,Y,Z);
    
    
    support_shift = DiffractionPatterns.shift_object(abs(NW*sqrt(mncntrate/mn)),support_plot,delta_thscanvals,ki_o,kf_o,kf_o-ki_o,d2_bragg,X,Y,Z);
    support_shift_abs = abs(support_shift);
    support_shift_fin = (support_shift_abs>0.1*max(support_shift_abs(:)));
    
    
    phase_rho_shift = angle(rho_shift(65,65,65));
    phase_NW = angle(NW(65,65,65));
    
    
    struct_toplot(kk).rho_shift = rho_shift*exp(-1i*phase_rho_shift).*support_shift_fin;
    struct_toplot(kk).support_shift_fin = support_shift_fin;
    struct_toplot(kk).rho_nophase = rho_shift.*conj(NW)*exp(-1i*phase_rho_shift).*support_shift_fin*exp(-1i*phaseoffset(kk));
    
    
    
    
    
    
end




fig_num = 40;
phase_color = [0 2.1];
phase_color_2 = [-0.1 0.1];
FiguresForPaper.figure5_bottompanel(struct_toplot,phase_color,phase_color_2,[40 90 40 90],[65],'3',fig_num);

save('rseults_sim_blueshift/ER_HIO_200iter','ER_HIOstruct','struct_toplot');
