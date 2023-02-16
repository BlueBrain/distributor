%plot non-neuronal cells

%%%%%%%%%%%
figure
temp_allvec_individual_vector5=allvec_individual_vector5;
temp_allvec_individual_vector5(:,1:7)=nan;
gab_ast=xlsread('Gabbott Astrocytes.csv');
gab_oli=xlsread('Gabbott Oligodendrocytes.csv');
gab_mic=xlsread('Gabbott Microglia.csv');
inds=find(gab_ast(:,2)>1323);
gab_ast(inds,:)=[];
gab_ast(:,2)=gab_ast(:,2)./1323*100;
inds=find(gab_oli(:,2)>1323);
gab_oli(inds,:)=[];
gab_oli(:,2)=gab_oli(:,2)./1323*100;
inds=find(gab_mic(:,2)>1323);
gab_mic(inds,:)=[];
gab_mic(:,2)=gab_mic(:,2)./1323*100;

mcc_ast=xlsread('McCaslin2011.csv');
lw=1;
subplot(1,3,1)
hold on
plot_allen_layers(astro_inds1,temp_allvec_individual_vector5,layer_array)
plot(transpose(temp_allvec_individual_vector5(astro_inds1,:)),'LineWidth',lw)
hold on
plot(gab_ast(:,2),gab_ast(:,1),'k','LineWidth',1.25)
max_ast=floor(max(gab_ast(:,2)));
ast_mean=mean(sum(allvec_individual_vector5(astro_inds1,8:100),1))
ast_mean_gab=mean(interp1(gab_ast(:,2),gab_ast(:,1),[8:max_ast]))
%plot(mcc_ast(:,1)./480*45,mcc_ast(:,2),'g','LineWidth',2)
ylim([0 (1.1*max(gab_ast(:,1)))])
pbaspect([1 1 1])
m1=ttype_labels(astro_inds1);
m1(2:length(m1)+1)=m1;
m1(1)={'Total Astro'};
for count1=1:length(m1)
    m1{count1}=strrep(m1{count1},'_',' ');
end
m1=[m1,"Gabbott"]
%set(gca,'FontSize', 2)
%legend(m1)
title('Astrocytes')

subplot(1,3,2)
hold on
plot_allen_layers(oligo_inds1,temp_allvec_individual_vector5,layer_array)
plot(transpose(temp_allvec_individual_vector5(oligo_inds1,:)),'LineWidth',lw)
pbaspect([1 1 1])
m1=ttype_labels(oligo_inds1);
m1(2:length(m1)+1)=m1;
m1(1)={'Total Oligo'};
for count1=1:length(m1)
    m1{count1}=strrep(m1{count1},'_',' ');
end
m1=[m1,"Gabbott"]
plot(gab_oli(:,2),gab_oli(:,1),'k','LineWidth',1.25)
max_oli=floor(max(gab_oli(:,2)));
oli_mean=mean(sum(allvec_individual_vector5(oligo_inds1,8:100),1))
oli_mean_gab=mean(interp1(gab_oli(:,2)./1500*100,gab_oli(:,1),[8:max_oli]))

ylim([0 (1.1*max(gab_oli(:,1)))])

%set(gca,'FontSize', 2)
%legend(m1)
title('Oligodendrocytes')

subplot(1,3,3)
hold on
inds=micro_inds1;
plot_allen_layers(micro_inds1,temp_allvec_individual_vector5,layer_array)
plot(transpose(temp_allvec_individual_vector5(micro_inds1,:)),'LineWidth',lw)
max_mic=floor(max(gab_mic(:,2)./1500*100));
mic_mean=mean(sum(allvec_individual_vector5(micro_inds1,8:100),1))
mic_max=max(sum(allvec_individual_vector5(micro_inds1,8:100),1));

mic_mean_gab=mean(interp1(gab_mic(:,2)./1500*100,gab_mic(:,1),[8:max_mic]))

pbaspect([1 1 1])
m1=ttype_labels(micro_inds1);
m1(2:length(m1)+1)=m1;
m1(1)={'Total Micro'};
for count1=1:length(m1)
    m1{count1}=strrep(m1{count1},'_',' ');
end
m1=[m1,"Gabbott"]
plot(gab_mic(:,2),gab_mic(:,1),'k','LineWidth',1.25)
ylim([0 (1.1*max(mic_max,max(gab_mic(:,1))))])

%set(gca,'FontSize', 2)
%legend(m1)
title('micro')

figure
subplot(2,3,4)
hold on
inds=smc_inds1;
plot_allen_layers(smc_inds1,temp_allvec_individual_vector5,layer_array)
plot(transpose(temp_allvec_individual_vector5(smc_inds1,:)),'LineWidth',lw)
pbaspect([1 1 1])
m1=ttype_labels(smc_inds1);
m1(2:length(m1)+1)=m1;
m1(1)={'Total SMC'};
for count1=1:length(m1)
    m1{count1}=strrep(m1{count1},'_',' ');
end
%set(gca,'FontSize', 2)
%legend(m1)
title('Smc')


subplot(2,3,5)
hold on
plot_allen_layers(vlmc_inds1,temp_allvec_individual_vector5,layer_array)
plot(transpose(temp_allvec_individual_vector5(vlmc_inds1,:)),'LineWidth',lw)

m1=ttype_labels(vlmc_inds1);
m1(2:length(m1)+1)=m1;
m1(1)={'Total Vlmc'};
for count1=1:length(m1)
    m1{count1}=strrep(m1{count1},'_',' ');
end
%set(gca,'FontSize', 2)
%legend(m1)
pbaspect([1 1 1])
title('Vlmc')

subplot(2,3,6)
hold on
inds=endo_inds1;
plot_allen_layers2(endo_inds1,temp_allvec_individual_vector5,layer_array)
pbaspect([1 1 1])
title('Endo')

drawnow
refresh
clear('temp_allvec_individual_vector5')
