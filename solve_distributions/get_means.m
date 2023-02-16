goodCells=sum(clusters0(goodinds,:)>0);

temp=zeros(1,numCells);
temp(l2_inds1)=1;
mean2=transpose(mean(transpose(clusters0(:,find(goodCells.*temp)))));


temp=zeros(1,numCells);
temp(l23_inds1)=1;
mean23=transpose(mean(transpose(clusters0(:,find(goodCells.*temp)))));
std23=transpose(std(transpose(clusters0(:,find(goodCells.*temp)))));
other=1-temp;
mean23Other=transpose(mean(transpose(clusters0(:,find(goodCells.*other)))));
std23Other=transpose(std(transpose(clusters0(:,find(goodCells.*other)))));
z23=erf( ((mean23-mean23Other)./sqrt(std23.*std23+std23Other.*std23Other)));
z23(isnan(z23))=0;

temp=zeros(1,numCells);
temp(l3_inds1)=1;
mean3=transpose(mean(transpose(clusters0(:,find(goodCells.*temp)))));


temp=zeros(1,numCells);
temp(l4_inds1)=1;
mean4=transpose(mean(transpose(clusters0(:,find(goodCells.*temp)))));
std4=transpose(std(transpose(clusters0(:,find(temp.*goodCells)))));
other=1-temp;
mean4Other=transpose(mean(transpose(clusters0(:,find(goodCells.*other)))));
std4Other=transpose(std(transpose(clusters0(:,find(goodCells.*other)))));
z4=erf( ((mean4-mean4Other)./sqrt(std4.*std4+std4Other.*std4Other)));
z4(isnan(z4))=0;

temp=zeros(1,numCells);
temp(l45_inds1)=1;
mean45=transpose(mean(transpose(clusters0(:,find(goodCells.*temp)))));
std45=transpose(std(transpose(clusters0(:,find(temp.*goodCells)))));
other=1-temp;
mean45Other=transpose(mean(transpose(clusters0(:,find(other.*goodCells)))));
std45Other=transpose(std(transpose(clusters0(:,find(other.*goodCells)))));
z45=erf( ((mean45-mean45Other)./sqrt(std45.*std45+std45Other.*std45Other)));
z45(isnan(z45))=0;

temp=zeros(1,numCells);
temp(l5_inds1)=1;
mean5=transpose(mean(transpose(clusters0(:,find(goodCells.*temp)))));
std5=transpose(std(transpose(clusters0(:,find(temp.*goodCells)))));
other=1-temp;
mean5Other=transpose(mean(transpose(clusters0(:,find(other.*goodCells)))));
std5Other=transpose(std(transpose(clusters0(:,find(other.*goodCells)))));
z5=erf( ((mean5-mean5Other)./sqrt(std5.*std5+std5Other.*std5Other)));
z5(isnan(z5))=0;

temp=zeros(1,numCells);
temp(l6_inds1)=1;
mean6=transpose(mean(transpose(clusters0(:,find(goodCells.*temp)))));
std6=transpose(std(transpose(clusters0(:,find(temp.*goodCells)))));
other=1-temp;
mean6Other=transpose(mean(transpose(clusters0(:,find(other.*goodCells)))));
std6Other=transpose(std(transpose(clusters0(:,find(other.*goodCells)))));
z6=erf( ((mean6-mean6Other)./sqrt(std6.*std6+std6Other.*std6Other)));
z6(isnan(z6))=0;

temp=zeros(1,numCells);
temp(l6b_inds1)=1;
mean6b=transpose(mean(transpose(clusters0(:,find(goodCells.*temp)))));
std6b=transpose(std(transpose(clusters0(:,find(temp.*goodCells)))));
other=1-temp;
mean6bOther=transpose(mean(transpose(clusters0(:,find(other.*goodCells)))));
std6bOther=transpose(std(transpose(clusters0(:,find(other.*goodCells)))));
z6b=erf( ((mean6b-mean6bOther)./sqrt(std6b.*std6b+std6bOther.*std6bOther)));
z6b(isnan(z6b))=0;

temp=zeros(1,numCells);
%astro_inds1=[272,273];
temp(astro_inds1)=1;
meanAstro=transpose(mean(transpose(clusters0(:,find(temp.*goodCells)))));
stdAstro=transpose(std(transpose(clusters0(:,find(temp.*goodCells)))));
other=1-temp;
meanAstroOther=transpose(mean(transpose(clusters0(:,find(other.*goodCells)))));
stdAstroOther=transpose(std(transpose(clusters0(:,find(other.*goodCells)))));
zAstro=erf( ((meanAstro-meanAstroOther)./sqrt(stdAstro.*stdAstro+stdAstroOther.*stdAstroOther)));
zAstro(isnan(zAstro))=0;
astroscalevec=(zAstro);

temp=zeros(1,numCells);
%oligo_inds1=[274:279];
temp(oligo_inds1)=1;
meanOligo=transpose(mean(transpose(clusters0(:,find(temp.*goodCells)))));
stdOligo=transpose(std(transpose(clusters0(:,find(temp.*goodCells)))));
other=1-temp;
meanOligoOther=transpose(mean(transpose(clusters0(:,find(other.*goodCells)))));
stdOligoOther=transpose(std(transpose(clusters0(:,find(other.*goodCells)))));
zOligo=erf( ((meanOligo-meanOligoOther)./sqrt(stdOligo.*stdOligo+stdOligoOther.*stdOligoOther)));
zOligo(isnan(zOligo))=0;
oliscalevec=(zOligo);

temp=zeros(1,numCells);
temp(vlmc_inds1)=1;
meanVlmc=transpose(mean(transpose(clusters0(:,find(temp.*goodCells)))));
stdVlmc=transpose(std(transpose(clusters0(:,find(temp.*goodCells)))));
other=1-temp;
meanVlmcOther=transpose(mean(transpose(clusters0(:,find(other.*goodCells)))));
stdVlmcOther=transpose(std(transpose(clusters0(:,find(other.*goodCells)))));
zVlmc=erf( ((meanVlmc-meanVlmcOther)./sqrt(stdVlmc.*stdVlmc+stdVlmcOther.*stdVlmcOther)));
vlmcscalevec=(zVlmc);

temp=zeros(1,numCells);
temp(endo_inds1)=1;
meanEndo=transpose(mean(transpose(clusters0(:,find(temp.*goodCells)))));
stdEndo=transpose(std(transpose(clusters0(:,find(temp.*goodCells)))));
other=1-temp;
meanEndoOther=transpose(mean(transpose(clusters0(:,find(other.*goodCells)))));
stdEndoOther=transpose(std(transpose(clusters0(:,find(other.*goodCells)))));
zEndo=erf( ((meanEndo-meanEndoOther)./sqrt(stdEndo.*stdEndo+stdEndoOther.*stdEndoOther)));
endoscalevec=(zEndo);

%Smc_ind1=[288];
temp=zeros(1,numCells);
temp(smc_inds1)=1;
meanSmc=transpose(mean(transpose(clusters0(:,find(temp.*goodCells)))));

stdSmc=0.*meanSmc;
temp=zeros(1,numCells);
temp(smc_inds1)=1;
other=1-temp;
meanSmcOther=transpose(mean(transpose(clusters0(:,find(other.*goodCells)))));
stdSmcOther=transpose(std(transpose(clusters0(:,find(other.*goodCells)))));
zSmc=erf( ((meanSmc-meanSmcOther)./sqrt(stdSmc.*stdSmc+stdSmcOther.*stdSmcOther+1e-9)));
Smcscalevec=(zSmc);

temp=zeros(1,numCells);
%micro_inds1=[289,numCells];
temp(micro_inds1)=1;
meanmicro=transpose(mean(transpose(clusters0(:,find(temp.*goodCells)))));
meanmicro=transpose(mean(transpose(clusters0(:,find(temp.*goodCells)))));
stdmicro=transpose(std(transpose(clusters0(:,find(temp.*goodCells)))));
other=1-temp;
meanmicroOther=transpose(mean(transpose(clusters0(:,find(other.*goodCells)))));
stdmicroOther=transpose(std(transpose(clusters0(:,find(other.*goodCells)))));
zmicro=erf( ((meanmicro-meanmicroOther)./sqrt(stdmicro.*stdmicro+stdmicroOther.*stdmicroOther)));
microscalevec=(zmicro);

otherscalevec=transpose(min(transpose(erf([zmicro,zSmc,zEndo,zVlmc]))));
%  scaleVec=(ones(length(clusters00),1)+i1*abs(inhscalevec-thr)+e1*abs(excscalevec-thr)+a1*abs(astroscalevec-thr)+o1*abs(oliscalevec-thr)+oth1.*abs(otherscalevec-thr));

%inhibitory subtypes

temp=zeros(1,numCells);
temp(pvalb_inds1)=1;
meanpvalb=transpose(mean(transpose(clusters0(:,find(goodCells.*temp)))));

temp=zeros(1,numCells);
temp(vip_inds1)=1;
meanvip=transpose(mean(transpose(clusters0(:,find(goodCells.*temp)))));

temp=zeros(1,numCells);
temp(lamp_inds1)=1;
meanlamp=transpose(mean(transpose(clusters0(:,find(goodCells.*temp)))));

temp=zeros(1,numCells);
temp(sncg_inds1)=1;
meansncg=transpose(mean(transpose(clusters0(:,find(goodCells.*temp)))));
%%%%%%%

temp=zeros(1,numCells);
temp(sst_inds1)=1;
meansst=transpose(mean(transpose(clusters0(:,find(goodCells.*temp)))));
stdsst=transpose(std(transpose(clusters0(:,find(temp.*goodCells)))));
other=zeros(1,numCells);
other(inh_inds1)=1;
other=other-temp;
meansstOther=transpose(mean(transpose(clusters0(:,find(other.*goodCells)))));
stdsstOther=transpose(std(transpose(clusters0(:,find(other.*goodCells)))));
zsst=erf( ((meansst-meansstOther)./sqrt(stdsst.*stdsst+stdsstOther.*stdsstOther)));
zsst(isnan(zsst))=0;

temp=zeros(1,numCells);
temp(inh_inds1)=1;
meanInh=transpose(mean(transpose(clusters0(:,find(goodCells.*temp)))));

temp=zeros(1,numCells);
temp(pvalb_inds1)=1;
meanpvalb=transpose(mean(transpose(clusters0(:,find(goodCells.*temp)))));
stdpvalb=transpose(std(transpose(clusters0(:,find(temp.*goodCells)))));
other=zeros(1,numCells);
other(inh_inds1)=1;
other=other-temp;
meanpvalbOther=transpose(mean(transpose(clusters0(:,find(other.*goodCells)))));
stdpvalbOther=transpose(std(transpose(clusters0(:,find(other.*goodCells)))));
zpvalb=erf( ((meanpvalb-meanpvalbOther)./sqrt(stdpvalb.*stdpvalb+stdpvalbOther.*stdpvalbOther)));
zpvalb(isnan(zpvalb))=0;


temp=zeros(1,numCells);
temp(htr_inds1)=1;
meanhtr=transpose(mean(transpose(clusters0(:,find(goodCells.*temp)))));


%now do nozero
meansst2=zeros(num_genes,1);
meanpvalb2=zeros(num_genes,1);
meanhtr2=zeros(num_genes,1);
meanneur2=zeros(num_genes,1);
meaninh2=zeros(num_genes,1);
meanexc2=zeros(num_genes,1);
meanendo2=zeros(num_genes,1);
meansmc2=zeros(num_genes,1);
meanmicro2=zeros(num_genes,1);
meanvlmc2=zeros(num_genes,1);
meanastro2=zeros(num_genes,1);
meanoligo2=zeros(num_genes,1);
meanall2=zeros(num_genes,1);
mean232=zeros(num_genes,1);
mean42=zeros(num_genes,1);
mean452=zeros(num_genes,1);
mean52=zeros(num_genes,1);
mean62=zeros(num_genes,1);
mean6b2=zeros(num_genes,1);

stdsst2=zeros(num_genes,1);
stdpvalb2=zeros(num_genes,1);
stdhtr2=zeros(num_genes,1);
stdneur2=zeros(num_genes,1);
stdinh2=zeros(num_genes,1);
stdexc2=zeros(num_genes,1);
stdendo2=zeros(num_genes,1);
stdSmc2=zeros(num_genes,1);
stdmicro2=zeros(num_genes,1);
stdvlmc2=zeros(num_genes,1);
stdastro2=zeros(num_genes,1);
stdoligo2=zeros(num_genes,1);
stdall2=zeros(num_genes,1);
std232=zeros(num_genes,1);
std42=zeros(num_genes,1);
std452=zeros(num_genes,1);
std52=zeros(num_genes,1);
std62=zeros(num_genes,1);
std6b2=zeros(num_genes,1);

for count1=1:length(goodinds)
    m1=clusters0(goodinds(count1),sst_inds1);
    meansst2(goodinds(count1))=mean(m1(m1>=0));
    stdsst2(goodinds(count1))=std(m1(m1>=0));
    m1=clusters0(goodinds(count1),pvalb_inds1);
    meanpvalb2(goodinds(count1))=mean(m1(m1>=0));
    stdpvalb2(goodinds(count1))=std(m1(m1>=0));
    m1=clusters0(goodinds(count1),htr_inds1);
    meanhtr2(goodinds(count1))=mean(m1(m1>=0));
    stdhtr2(goodinds(count1))=std(m1(m1>=0));
    
    m1=clusters0(goodinds(count1),exc_inds1);
    meanexc2(goodinds(count1))=mean(m1(m1>=0));
    stdexc2(goodinds(count1)) =std(m1(m1>=0));
    m1=clusters0(goodinds(count1),inh_inds1);
    meaninh2(goodinds(count1))=mean(m1(m1>=0));
    stdinh2(goodinds(count1))=std(m1(m1>=0));
    
    m1=clusters0(goodinds(count1),vlmc_inds1);
    meanvlmc2(goodinds(count1))=mean(m1(m1>=0));
    stdvlmc2(goodinds(count1))=std(m1(m1>=0));
    m1=clusters0(goodinds(count1),endo_inds1);
    meanendo2(goodinds(count1))=mean(m1(m1>=0));
    stdendo2(goodinds(count1))=std(m1(m1>=0));
    m1=clusters0(goodinds(count1),smc_inds1);
    meansmc2(goodinds(count1))=mean(m1(m1>=0));
    
    m1=clusters0(goodinds(count1),micro_inds1);
    meanmicro2(goodinds(count1))=mean(m1(m1>=0));
    stdmicro2(goodinds(count1))=std(m1(m1>=0));
    
    m1=clusters0(goodinds(count1),astro_inds1);
    meanastro2(goodinds(count1))=mean(m1(m1>=0));
    stdastro2(goodinds(count1))=std(m1(m1>=0));
    m1=clusters0(goodinds(count1),oligo_inds1);
    meanoligo2(goodinds(count1))=mean(m1(m1>=0));
    stdoligo2(goodinds(count1))=std(m1(m1>=0));
    
    m1=clusters0(goodinds(count1),neur_inds1);
    meanneur2(goodinds(count1))=mean(m1(m1>=0));
    stdneur2(goodinds(count1))=std(m1(m1>=0));
    
    m1=clusters0(goodinds(count1),:);
    meanall2(goodinds(count1))=mean(m1(m1>0));
    stdall2(goodinds(count1))=std(m1(m1>=0));
    
    m1=clusters0(goodinds(count1),l23_inds1);
    mean232(goodinds(count1))=mean(m1(m1>0));
    std232(goodinds(count1))=std(m1(m1>=0));
    m1=clusters0(goodinds(count1),l4_inds1);
    mean42(goodinds(count1))=mean(m1(m1>0));
    std42(goodinds(count1))=std(m1(m1>=0));
    m1=clusters0(goodinds(count1),l45_inds1);
    mean452(goodinds(count1))=mean(m1(m1>0));
    std452(goodinds(count1))=std(m1(m1>=0));
    m1=clusters0(goodinds(count1),l5_inds1);
    mean52(goodinds(count1))=mean(m1(m1>0));
    std52(goodinds(count1))=std(m1(m1>=0));
    
    m1=clusters0(goodinds(count1),l6_inds1);
    mean62(goodinds(count1))=mean(m1(m1>0));
    std62(goodinds(count1))=std(m1(m1>=0));
    
    m1=clusters0(goodinds(count1),l6b_inds1);
    mean6b2(goodinds(count1))=mean(m1(m1>0));
    std6b2(goodinds(count1))=std(m1(m1>=0));
end
%meanhtr2=transpose(mean(transpose(clusters0(:,htr_inds1))));
%meanpv2=transpose(mean(transpose(clusters0(:,pvalb_inds1))));
%meanInh2=0.4.* pv_genes2+0.3.*sst_genes2+0.3.*htr_genes2;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%clusters2=[meanNeur,meanAstro,meanOligo,meanVlmc,meanEndo,meanSmc,meanmicro];
%meanall=transpose(mean(transpose(clusters2)));

ttype_labels={'Neurons','Astro','Oligo','VLMC','Endo','Smc','micro'}  ;
%basemat=[meanall,meanall,meanall,meanall,meanall,meanall  ,meanall ];

mat2=[meanneur2,meanastro2,meanoligo2,meanvlmc2,meanmicro2,meansmc2,meanendo2];
stdmat2=[stdneur2,stdastro2,stdoligo2,stdvlmc2,stdmicro2,stdvlmc2,stdendo2];

meanall2=zeros(num_genes,1);
for count1=1:length(goodinds)
    m1=mat2(goodinds(count1),:);
    meanall2(goodinds(count1))=mean(m1(m1>=0));
end

meanother=transpose(mean(transpose(clusters0(:,[vlmc_inds1;peri_inds1;micro_inds1;endo_inds1;smc_inds1]))));

