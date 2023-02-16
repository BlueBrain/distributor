%plot both exc and inh 

inhVec1=sum(allvec_individual_vector5(inh_inds1,:),1);
excVec1=sum(allvec_individual_vector5(exc_inds1,:),1);
figure
hold on
plot(inhVec1+excVec1,'k')
plot(excVec1,'g')
plot(inhVec1,'r')
max1=1.1 * max(inhVec1+excVec1);
plot([l12 l12],[0 max1],'--','color',[.7 .7 .7])
plot([l23 l23],[0 max1],'--','color',[.7 .7 .7])
plot([l34 l34],[0 max1],'--','color',[.7 .7 .7])
plot([l45 l45],[0 max1],'--','color',[.7 .7 .7])
plot([l56 l56],[0 max1],'--','color',[.7 .7 .7])
ylim([0 max1])
%legend({'Excitatory','Inhibitory','All Neurons'},'Location','northwest')
