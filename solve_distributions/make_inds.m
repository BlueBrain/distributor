load ttype_labels

uncorrected_pvalb_inds1=[108:123];
pvalb_inds1 = correct_inds(uncorrected_pvalb_inds1,ttype_labels,clusters0,goodinds);

uncorrected_sst_inds1=[63:107];
sst_inds1 = correct_inds(uncorrected_sst_inds1,ttype_labels,clusters0,goodinds);

uncorrected_vip_inds1=[40:62];
vip_inds1 = correct_inds(uncorrected_vip_inds1,ttype_labels,clusters0,goodinds);

uncorrected_lamp_inds1=[5:18];
lamp_inds1 = correct_inds(uncorrected_lamp_inds1,ttype_labels,clusters0,goodinds);

uncorrected_sncg_inds1=[19:39];
sncg_inds1 = correct_inds(uncorrected_sncg_inds1,ttype_labels,clusters0,goodinds);
inh_inds1=[pvalb_inds1;sst_inds1;vip_inds1 ;lamp_inds1 ;sncg_inds1];

%uncorrected_exc_inds1=[303:354,258:262,234:257,126:233];
%exc_inds1 = correct_inds(uncorrected_exc_inds1,ttype_labels,clusters0,goodinds);


vlmc_inds1=[];
peri_inds1=[];
oligo_inds1=[];
micro_inds1=[];
endo_inds1=[];
astro_inds1=[];
smc_inds1=[];

l2_inds1=[];
l23_inds1=[];
l3_inds1=[];
l4_inds1=[];
l3_inds1=[];
l45_inds1=[];
l5_inds1=[];
l6_inds1=[];
l6b_inds1=[];

for count1=1:length(ttype_labels)
    if (sum(clusters0(goodinds,count1))>0)
        if(length(strfind(ttype_labels{count1},'VLMC')))
            vlmc_inds1=[vlmc_inds1;count1];
        end
        if(length(strfind(ttype_labels{count1},'Peri')))
            peri_inds1=[peri_inds1;count1];
        end
             if(length(strfind(ttype_labels{count1},'SMC')))
            smc_inds1=[smc_inds1;count1];
        end
        if(length(strfind(ttype_labels{count1},'Oligo')))
            oligo_inds1=[oligo_inds1;count1];
        end
        if(length(strfind(ttype_labels{count1},'Micro')))
            micro_inds1=[micro_inds1;count1];
        end
        if(length(strfind(ttype_labels{count1},'Endo')))
            endo_inds1=[endo_inds1;count1];
        end
        if(length(strfind(ttype_labels{count1},'Astro')))
            astro_inds1=[astro_inds1;count1];
        end
        if(length(strfind(ttype_labels{count1},'L2')))
            if(length(strfind(ttype_labels{count1},'L2_3')))
                l23_inds1=[l23_inds1;count1];
            else
                l2_inds1=[l2_inds1;count1];
            end
        end
        if(length(strfind(ttype_labels{count1},'L3')))
            l3_inds1=[l3_inds1;count1];            
        end
        if(length(strfind(ttype_labels{count1},'L4')))
            
            if(length(strfind(ttype_labels{count1},'L4_5')))
                
                l45_inds1=[l45_inds1;count1];
            else
                l4_inds1=[l4_inds1;count1];
            end
        end
        if(length(strfind(ttype_labels{count1},'L5')))
            l5_inds1=[l5_inds1;count1];
        end
        if(length(strfind(ttype_labels{count1},'L6')))
            if(length(strfind(ttype_labels{count1},'L6b')))
                l6b_inds1=[l6b_inds1;count1];
            else
                l6_inds1=[l6_inds1;count1];
            end
        end
    end
end
exc_inds1=[l2_inds1;l23_inds1;l3_inds1;l4_inds1;l45_inds1;l5_inds1;l6_inds1;l6b_inds1];

htr_inds1=[];
vip_only_inds1=[];
sst_only_inds1=[];
pvalb_only_inds1=[];

%sst_genes=transpose(mean(transpose(clusters0(:,sst_inds1))));
%htr_genes=transpose(mean(transpose(clusters0(:,htr_inds1))));
%pv_genes=transpose(mean(transpose(clusters0(:,pvalb_inds1))));
neur_inds1=[inh_inds1;exc_inds1];

placed_types3=find(sum( clusters0(goodinds,:))>0);%[1:290];
vec1(placed_types3)=1;
vec2=zeros(1,length(vec1));
validCellInds=sort([inh_inds1;exc_inds1;peri_inds1;astro_inds1;oligo_inds1;vlmc_inds1;micro_inds1;smc_inds1;endo_inds1]);
vec2(validCellInds)=1;
validCellInds=find(vec1.*vec2);

[foo numCells]=size(clusters0);
validVec=zeros(1,numCells);
validVec(validCellInds)=1;
goodCells=(validVec);

[a b]=size(clusters0);
temp=zeros(1,b);
temp(inh_inds1)=1;
temp(pvalb_inds1)=0;
temp(sst_inds1)=0;
inhib_minus_pv_sst_inds1=find(temp);

%changed to clusters0
vec1=1.*(clusters0(4358,:)>1);
vec1(inh_inds1)=vec1(inh_inds1)+1;
CRexprInds=find(vec1==2);

%pvalb
vec1=1*(clusters0(38359,:)>=0.65);
vec1(inh_inds1)=vec1(inh_inds1)+1;
%vec1(SSTexprInds)=0;
PvalbexprInds=find(vec1==2);

%sst
vec1=1*(clusters0(41532,:)>=0.69);
vec1(inh_inds1)=vec1(inh_inds1)+1;
vec1(PvalbexprInds)=0;
vec1(CRexprInds)=0;
%vec1(pvalb_inds1)=0;
SSTexprInds=find(vec1==2);

localinds2=[l2_inds1;l23_inds1;l3_inds1;l4_inds1;l45_inds1;l5_inds1;l6_inds1;l6b_inds1;inh_inds1;astro_inds1;oligo_inds1;vlmc_inds1;micro_inds1;smc_inds1;peri_inds1;endo_inds1];

%htra3
%vec1=1.*(sum(clusters00([669,935,5026,7927,7981,8559,8706,11067,14200,16320,18818,19696,19697,20151,24095,27304,29992,30675],:))>0.0);
vec1=1.*((clusters0([24436],:))>2.5);
vec1(inh_inds1)=vec1(inh_inds1)+1;
HTRexprInds=find(vec1==2);

temp1=zeros(1,b);
temp1(localinds2)=1;
temp2=temp1;
temp2(pvalb_inds1)=temp2(pvalb_inds1)+1;
pvalb_inds2=find(temp2==2);
temp2=temp1;
temp2(sst_inds1)=temp2(sst_inds1)+1;
sst_inds2=find(temp2==2);
temp2=temp1;
temp2(vip_inds1)=temp2(vip_inds1)+1;
vip_inds2=find(temp2==2);
temp2=temp1;
temp2(lamp_inds1)=temp2(lamp_inds1)+1;
lamp_inds2=find(temp2==2);
temp2=temp1;
temp2(sncg_inds1)=temp2(sncg_inds1)+1;
sncg_inds2=find(temp2==2);

temp2=temp1;
temp2(inh_inds1)=temp2(inh_inds1)+1;
inh_inds2=find(temp2==2);
temp2=temp1;
temp2(exc_inds1)=temp2(exc_inds1)+1;
exc_inds2=find(temp2==2);
temp2=temp1;
temp2(astro_inds1)=temp2(astro_inds1)+1;
astro_inds2=find(temp2==2);
temp2=temp1;
temp2(oligo_inds1)=temp2(oligo_inds1)+1;
oligo_inds2=find(temp2==2);
temp=zeros(1,400);
temp(inh_inds2)=1;
temp(pvalb_inds2)=0;
inh_minus_pvalb_inds=find(temp);
temp(sst_inds2)=0;
inh_minus_pvalb_sst_inds=find(temp);
other_inds1=[vlmc_inds1;smc_inds1;endo_inds1;peri_inds1;micro_inds1];
temp2=temp1;
temp2(other_inds1)=temp2(other_inds1)+1;
other_inds2=find(temp2==2);

inhvec0=zeros(1,400);
inhvec0(inh_inds1)=1;
inhvec0(CRexprInds)=0;
inhvec0(PvalbexprInds)=0;
inhvec0(SSTexprInds)=0;
other_inh_inds1=find(inhvec0);