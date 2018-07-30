% In this script the figures for the paper are generated:

%%%%% Figure 1 %%%%%%%%%%%%%%


fig_num = 1;
slice_array = [1:floor(round(numel(angles_list)/4)):numel(angles_list)];
[Psi_mod_matrix] = FiguresForPaper.display_slice_dp(NW,probe,delta_thscanvals,slice_array,ki_o,kf_o,d2_bragg,th,X,Y,Z,fig_num);

figure(2);
clf;

for jj = 1:numel(slice_array)
    %subplot(1,numel(slice_array),jj);
    figure;
    imagesc(Psi_mod_matrix(:,:,slice_array(jj)));
    axis image;
    colorbar;
    title(num2str(slice_array(jj)));
end

fig_num = 3;
DisplayFunctions.display_diff_geom(NW,ki,kf,qbragg,fig_num,X,Y,Z);




%%%%%%%%%%%% Figure 1: low panel %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%% Figure 5: comparison of the quality of the reconstruction for
%%%%%%%% different jittering and different levels of noise

jitterlevel_summary = [0 5 10 20 40];
noiselevel_str = '2';
for kk = 1:numel(jitterlevel_summary)
    
    %savefolder = ['jitter_' num2str(jitterlevel_summary(kk)) '_noiselevel_' noiselevel_str '_freqsw500'];
    savefolder = ['allresults_blueshift/jitter_' num2str(jitterlevel_summary(kk)) '_noiselevel_' noiselevel_str '_70angles'];
    load([savefolder '/results.mat']);
    load(['data_ERHIO/struct_ERHIO_ini' noiselevel_str '_jitter_' num2str(jitterlevel_summary(kk))]);

    struct_err(kk).chi = [struct_best_ERHIO.chi' errlist];
    %struct_err(kk).chid_direct = [err_ERHIO errlist_direct];
    struct_err(kk).rho = rho;
    struct_err(kk).support_iter = support_iter;
    struct_err(kk).rho_3DFT = rho_3DFT;
    struct_err(kk).support_new = support_new;
end

%mkdir results_sim_blueshift;

save(['allresults_blueshift/struct_err_level' noiselevel_str '.mat'],'struct_err');


figure(100);
clf;
hold on;

for kk = 1:numel(struct_err)
   
    chi_final(kk) = struct_err(kk).chi(end);
    %chi_direct_final(kk) = struct_err(kk).chid_direct(end);
    
    plot(log10(struct_err(kk).chi),'LineWidth',3.0);
    legend_str{kk} = [num2str(jitterlevel_summary(kk)) ];
end
legend(legend_str);

figure(100);
savefig(['allresults_blueshift/error_chi_alljiiter_' noiselevel_str '.fig']);



figure(1);
%subplot(121);
plot(jitterlevel_summary,log10(chi_final),'-o');
%title('recip');
xlabel('% of angular jitter');ylabel('log(\epsilon)');
ax = gca;
set(ax,'FontSize',20);

figure(1);
savefig(['allresults_blueshift/chi_vs_jitter_levelnoise_' noiselevel_str '.fig']);


% subplot(122);
% plot(jitterlevel_summary,log10(chi_direct_final));
% title('direct');

% for noiselevel = 1;
%flipflag_list = [1 0 1 0];
%phaseoffset = [1.55 1.58 1.43 1.57];

% for noiselevel = 0;
phaseoffset = [1.48 1.59 1.52 1.59 1.51];%[1.55 1.58 1.58 1.545 1.01];
flipflag_list = [1 0 0 1 0];%[1 1 1 0 1];

for kk = 1:numel(struct_err)
    
    flipflag = flipflag_list(kk);
    
    if flipflag
        rho_plot = ifftn(conj(fftn(struct_err(kk).rho)));
        support_plot = abs(ifftn(conj(fftn(struct_err(kk).support_iter))));
    else
        rho_plot = struct_err(kk).rho;
        support_plot = struct_err(kk).support_iter;
    end
    
    
    rho_shift = DiffractionPatterns.shift_object(NW*sqrt(mncntrate/mn),rho_plot,delta_thscanvals,ki_o,kf_o,kf_o-ki_o,d2_bragg,X,Y,Z);
    
    
    support_shift = DiffractionPatterns.shift_object(abs(NW*sqrt(mncntrate/mn)),support_plot,delta_thscanvals,ki_o,kf_o,kf_o-ki_o,d2_bragg,X,Y,Z);
    support_shift_abs = abs(support_shift);
    support_shift_fin = (support_shift_abs>0.1*max(support_shift_abs(:)));
    
    midpoint_1 = [round(size(rho_shift,1)/2)+1 round(size(rho_shift,2)/2)+1 round(size(rho_shift,3)/2)+1];


    phase_rho_shift = angle(rho_shift(midpoint_1(1),midpoint_1(2),midpoint_1(3)));
    phase_NW = angle(NW(midpoint_1(1),midpoint_1(2),midpoint_1(3)));
    
    
    
    struct_toplot(kk).rho_shift = rho_shift*exp(-1i*phase_rho_shift).*support_shift_fin;
    struct_toplot(kk).support_shift_fin = support_shift_fin;
    struct_toplot(kk).rho_nophase = rho_shift.*conj(NW)*exp(-1i*phase_rho_shift).*support_shift_fin*exp(-1i*phaseoffset(kk));
    
    
    
end

DisplayResults.compare_two_objects(NW*sqrt(mncntrate/mn).*conj(NW),struct_toplot(5).rho_nophase,'','',[40 90],[midpoint_1(2) midpoint_1(3)],'23',42);


fig_num = 30;
phase_color = [0 2.0];
phase_color_2 = [-0.5 0.5];
intenscolor = [0 0.35];
dimension = 3;
FiguresForPaper.figure5_bottompanel(struct_toplot,intenscolor,phase_color,phase_color_2,[40 90 40 90],[midpoint_1(dimension)],num2str(dimension),fig_num);
% check phase offset

DisplayResults.compare_two_objects(NW*sqrt(mncntrate/mn).*conj(NW),struct_toplot(1).rho_nophase,'','',[1 128],[midpoint_1(2) midpoint_1(3)],'23',42);
DisplayResults.compare_two_objects(NW*sqrt(mncntrate/mn).*conj(NW),struct_toplot(2).rho_nophase,'','',[40 90],[midpoint_1(2) midpoint_1(3)],'23',42);
DisplayResults.compare_two_objects(NW*sqrt(mncntrate/mn).*conj(NW),struct_toplot(3).rho_nophase,'','',[40 90],[midpoint_1(2) midpoint_1(3)],'23',42);
DisplayResults.compare_two_objects(NW*sqrt(mncntrate/mn).*conj(NW),struct_toplot(4).rho_nophase,'','',[40 90],[midpoint_1(2) midpoint_1(3)],'23',42);
DisplayResults.compare_two_objects(NW*sqrt(mncntrate/mn).*conj(NW),struct_toplot(5).rho_nophase,'','',[40 90],[midpoint_1(2) midpoint_1(3)],'23',42);


figure(fig_num);
savefig(['allresults_blueshift/rho_alljitter_noiselevel' noiselevel_str '.fig']);

save(['allresults_blueshift/rhostructtoplot_level' noiselevel_str '.mat'],'struct_toplot','flipflag_list','phaseoffset');



%%%%%%%%%%%%%%%% Figure 2: effect of jittering in ER_HIO results


% for noiselevel = 1;
flipflag_list = [0 0 1 1 1];
%phaseoffset = [1.55 1.58 1.43 1.57 1.57];

% for noiselevel = 0;



phaseoffset_ERHIO =  [1.52 1.64 1.68 1.53 1.57];%[1.59 1.65 1.67 1.57 1.40];
%flipflag_list = [1 1 1 0 1];

for kk = 1:numel(struct_err)
    
    flipflag = flipflag_list(kk);
    
    if flipflag
        rho_plot = ifftn(conj(fftn(struct_err(kk).rho_3DFT)));
        support_plot = abs(ifftn(conj(fftn(struct_err(kk).support_new))));
    else
        rho_plot = struct_err(kk).rho_3DFT;
        support_plot = struct_err(kk).support_new;
    end
    
    
    rho_shift = DiffractionPatterns.shift_object(NW*sqrt(mncntrate/mn),rho_plot,delta_thscanvals,ki_o,kf_o,kf_o-ki_o,d2_bragg,X,Y,Z);
    
    
    support_shift = DiffractionPatterns.shift_object(abs(NW*sqrt(mncntrate/mn)),support_plot,delta_thscanvals,ki_o,kf_o,kf_o-ki_o,d2_bragg,X,Y,Z);
    support_shift_abs = abs(support_shift);
    support_shift_fin = (support_shift_abs>0.1*max(support_shift_abs(:)));
    
   midpoint_1 = [round(size(rho_shift,1)/2)+1 round(size(rho_shift,2)/2)+1 round(size(rho_shift,3)/2)+1];

    
     phase_rho_shift = angle(rho_shift(midpoint_1(1),midpoint_1(2),midpoint_1(3)));
    phase_NW = angle(NW(midpoint_1(1),midpoint_1(2),midpoint_1(3)));
           
    
    struct_toplot(kk).rho_shift = rho_shift*exp(-1i*phase_rho_shift).*support_shift_fin;
    struct_toplot(kk).support_shift_fin = support_shift_fin;
    struct_toplot(kk).rho_nophase = rho_shift.*conj(NW)*exp(-1i*phase_rho_shift).*support_shift_fin*exp(-1i*phaseoffset_ERHIO(kk));
    
    %DisplayResults.compare_two_objects(NW*sqrt(mncntrate/mn).*conj(NW),struct_toplot(kk).rho_nophase,'','',[40 90],[65 65],'23',32+kk);

    
end

DisplayResults.compare_two_objects(NW*sqrt(mncntrate/mn).*conj(NW),struct_toplot(5).rho_nophase,'','',[40 90],[midpoint_1(2) midpoint_1(3)],'23',42);


fig_num = 40;
phase_color = [0 2.0];
phase_color_2 = [-0.5 0.5];
intenscolor = [0 0.17];
dimension = 3;
FiguresForPaper.figure5_bottompanel(struct_toplot,intenscolor,phase_color,phase_color_2,[40 90 40 90],[midpoint_1(dimension)],num2str(dimension),fig_num);

figure(40);
savefig(['allresults_blueshift/rhoERHIO_alljitter_noiselevel' noiselevel_str '.fig']);


