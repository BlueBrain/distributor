%plot inhibitory marker subtypes

inhSum=sum(sum(allvec_individual_vector5(inh_inds1,:)));
sstFract=sum(sum(allvec_individual_vector5((SSTexprInds),:)))/inhSum

pvalbFract=sum(sum(allvec_individual_vector5((PvalbexprInds),:)))/inhSum

meyerCR=xlsread('meyer_CR.csv');
meyerSST=xlsread('meyer_SST.csv');
meyerPV=xlsread('meyer_PV.csv');

figure
hold on
my_inh=sum(allvec_individual_vector5(inh_inds1,:));
my_htr=sum(allvec_individual_vector5(HTRexprInds,:));
plot(my_htr,'k','LineWidth',2);
my_pv=sum(allvec_individual_vector5(PvalbexprInds,:));
plot(my_pv,'r','LineWidth',2);
my_sst=sum(allvec_individual_vector5(SSTexprInds,:));
plot(my_sst,'g','LineWidth',2);
my_cr=sum(allvec_individual_vector5(CRexprInds,:));
plot(my_cr,'b','LineWidth',2);

scale1=1;%0.05*1e-6;
%figure
hold on
plot(meyerPV(:,1),scale1.*meyerPV(:,2),'r')
plot(meyerSST(:,1),scale1.*meyerSST(:,2),'g')
plot(meyerCR(:,1),scale1.*meyerCR(:,2),'b')
xlim([0 100])
max1=1.1* max(max(my_pv),max(my_sst));


hold on
plot([l12 l12],[0 max1],'--','color',[.7 .7 .7])
plot([l23 l23],[0 max1],'--','color',[.7 .7 .7])
plot([l34 l34],[0 max1],'--','color',[.7 .7 .7])
plot([l45 l45],[0 max1],'--','color',[.7 .7 .7])
plot([l56 l56],[0 max1],'--','color',[.7 .7 .7])

ylim([0 ( max1)])
make_bar_charts2


figure
subplot(2,4,1)
plot(sum((squeeze(sum(allvec_individual6(l2_inds1,:,:),1)))))
title('L2')
set(gca,'ytick',[])

subplot(2,4,2)
plot(sum((squeeze(sum(allvec_individual6(l23_inds1,:,:),1)))))
title('L23')
set(gca,'ytick',[])
subplot(2,4,3)
plot(sum((squeeze(sum(allvec_individual6(l3_inds1,:,:),1)))))
title('L3')
set(gca,'ytick',[])

subplot(2,4,4)
plot(sum((squeeze(sum(allvec_individual6(l4_inds1,:,:),1)))))
title('L4')
set(gca,'ytick',[])

subplot(2,4,5)
plot(sum((squeeze(sum(allvec_individual6(l45_inds1,:,:),1)))))
title('L45')
set(gca,'ytick',[])

subplot(2,4,6)
plot(sum((squeeze(sum(allvec_individual6(l5_inds1,:,:),1)))))
title('L5')
set(gca,'ytick',[])

subplot(2,4,7)
plot(sum((squeeze(sum(allvec_individual6(l6_inds1,:,:),1)))))
title('L6')
set(gca,'ytick',[])

subplot(2,4,8)
plot(sum((squeeze(sum(allvec_individual6(l6b_inds1,:,:),1)))))
title('L6b')
set(gca,'ytick',[])