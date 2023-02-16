%plot neurons
%%%%%%%%%%%%%%%%%%%%%%%%%

figure
subplot(3,2,1)
hold on
%plot(sum((allvec_individual_vector5(sst_inds1,:))),'g')
plot_allen_layers(sst_inds1,allvec_individual_vector5,layer_array)

title('SST')

subplot(3,2,2)
hold on
%plot(sum((allvec_individual_vector5(pvalb_inds1,:))),'r')
plot_allen_layers(pvalb_inds1,allvec_individual_vector5,layer_array)

title('Pvalb')

subplot(3,2,3)
hold on
plot(sum((allvec_individual_vector5(htr_inds1,:))),'k')
title('HTR')

subplot(3,2,3)
hold on
plot_allen_layers(vip_inds1,allvec_individual_vector5,layer_array)
title('VIP')

subplot(3,2,4)
hold on
plot_allen_layers(lamp_inds1,allvec_individual_vector5,layer_array)
title('Lamp')

subplot(3,2,5)
hold on
inds=sncg_inds1;
plot_allen_layers(sncg_inds1,allvec_individual_vector5,layer_array)
title('Sncg')

subplot(3,2,6)
hold on
%m1=squeeze((expression_matrix4(num_genes+1,:,:)));
plot_allen_layers(inh_inds1,allvec_individual_vector5,layer_array)
%plot(sum((allvec_individual_vector5(inh_inds1,:))),'Color',[0.2,0.7,0.9])
title('Inhibitory')

%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
%subplot(3,4,1)
hold on
plot_allen_layers(neur_inds1,allvec_individual_vector5,layer_array)
title('Neurons')
