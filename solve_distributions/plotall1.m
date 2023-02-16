%make all of the plots
l12=7;
l23=20;
l34=35;
l45=45;
l56=70;

if(exist('map.mat'))
    load('map')
else
    map=rand(200,3);
    save('map','map')
end

validCellInds=find(goodCells);
inhRatio5=sum(sum(allvec_individual_vector5(inh_inds1,:)))./sum(sum(allvec_individual_vector5(neur_inds1,:)))
astRatio5=sum(sum(allvec_individual_vector5(astro_inds1,:)))./sum(sum(allvec_individual_vector5(neur_inds1,:)))
oliRatio5=sum(sum(allvec_individual_vector5(oligo_inds1,:)))./sum(sum(allvec_individual_vector5(neur_inds1,:)))
nrnRatio5=sum(sum(allvec_individual_vector5(neur_inds1,:)))./sum(sum(allvec_individual_vector5(validCellInds,:)))
sstRatio5=sum(sum(allvec_individual_vector5(sst_inds1,:)))./sum(sum(allvec_individual_vector5(inh_inds1,:)))
pvalbRatio5=sum(sum(allvec_individual_vector5(pvalb_inds1,:)))./sum(sum(allvec_individual_vector5(inh_inds1,:)))
otherRatio5=sum(sum(allvec_individual_vector5([vlmc_inds1;smc_inds1;endo_inds1;peri_inds1;micro_inds1],:)))./sum(sum(allvec_individual_vector5(validCellInds,:)));
err1=(5*(inhRatio5-0.2)./0.2)^2+((astRatio5-0.2)./0.2)^2+((oliRatio5-0.885)./0.885)^2+((nrnRatio5-0.46)./0.46)^2+((sstRatio5-0.3)./0.3)^2+((pvalbRatio5-0.4)./0.4)^2;

matchlist=[];

%%%%%%%%%
plot_inh_and_exc
goodinds3=goodinds2;
load ttype_labels
plot_exc
plot_neurons
plot_non_neuronal
plot_inhibitory_profiles5
plot_marker_fractions
