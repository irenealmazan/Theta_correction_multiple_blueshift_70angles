% this script checks that the center of mass coincides with the center of
% the array

X_com = sum(abs(NW(:)).*X(:))/sum(abs(NW(:)));
Y_correct = max(Y)-Y;
Y_com = sum(abs(NW(:)).*Y(:))/sum(abs(NW(:)));
Z_com = sum(abs(NW(:)).*Z(:))/sum(abs(NW(:)));


NW_shift = circshift(NW,-round([Y_com,X_com,Z_com]));

[NW_shift] = DiffractionPatterns.shift_object_known_shift(NW,-[X_com,Y_com,Z_com],delta_thscanvals,ki_o,kf_o,kf_o-ki_o,X,d2_bragg);


X_com_center = sum(abs(NW_shift(:)).*X(:))/sum(abs(NW_shift(:)))/d2_bragg;
Y_com_center = sum(abs(NW_shift(:)).*Y(:))/sum(abs(NW_shift(:)))/d2_bragg;
Z_com_center = sum(abs(NW_shift(:)).*Z(:))/sum(abs(NW_shift(:)))/th_pixel_size;