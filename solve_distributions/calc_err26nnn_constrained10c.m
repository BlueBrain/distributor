%function [allvec_individual,errvec,error,ratiovec,ishratiovec,baseratiovec]= calc_err23(scaleVec,exc_inds1,goodinds,expression_matrix,expression_matrix0,baseline_Scale,ish_Scale,clusters0,num_ttypes,placed_types,base_mat);

function [allvec_individual,errvec,error,ratiovec,ishratiovec,baseratiovec]= calc_err23(scaleVec,exc_mat,exc_inds1,goodinds,expression_matrix,expression_matrix0,baseline_Scale,ish_Scale,clusters0,num_ttypes,placed_types,base_mat,inh_mat,cr_mat,pv_mat,sst_mat,inh_inds1,sumAstro,sumOligo,sumVlmc,sumMicro,sumSmc,sumPeri,sumEndo,l2_inds1,l23_inds1,l3_inds1,l4_inds1,l45_inds1,l5_inds1,l6_inds1,l6b_inds1,astro_inds1,oligo_inds1,vlmc_inds1,micro_inds1,smc_inds1,peri_inds1,endo_inds1,CRexprInds,PvalbexprInds,SSTexprInds,neur_inds1,neur_mat);
[foo numcells]=size(clusters0);
placed_types=find(sum(clusters0(goodinds,:))>0);
allvec_individual=zeros(numcells,100,40);
%allvec_individual5(localinds,:,:)=allvec_individual_temp5;

placed_types_no_exc=[];
for count1=1:length(placed_types)
    if(sum(placed_types(count1)==exc_inds1)==0)
        placed_types_no_exc=[placed_types_no_exc;placed_types(count1)];
    end
end

num_ttypes=length(placed_types);
%allvec_individual=zeros(num_ttypes,100,40);
errvec=0.*ish_Scale;
ratiovec=0.*ish_Scale;
expr_base=zeros(length(goodinds),1);
pred_err2=zeros(length(goodinds),1);
sum1=0;

%options = optimset('Display','notify','TolX',1e-13);
%%options = optimset('Display','notify','TolX',1e-11);
%options = optimset('Display','notify','TolX',1e-5);
options = optimset('Display','notify','TolX',1e0);

%enforceVal=10000000;
enforceVal= 1000;
error=0;
%goodCellInds=find(sum( clusters0(goodinds(length(goodinds)-3),:)));
pred_genes=zeros(length(goodinds),100,40);
for count1=1:100
    cellinds=placed_types;
    if(count1<=6)
        cellinds=placed_types_no_exc;
    end
    if(1)
    inhvec0=zeros(1,max(cellinds));
    inhvec0(inh_inds1)=1;
    inhvec1=inhvec0(cellinds);
    excvec0=zeros(1,max(cellinds));
    eqz1=find(sum(clusters0(goodinds,cellinds))==0);
    excvec0(exc_inds1)=1;
    excvec0(eqz1)=0;
    excvec1=excvec0(cellinds);

    neurvec0=zeros(1,max(cellinds));
    neurvec0(neur_inds1)=1;
    neurvec0(eqz1)=0;
    neurvec1=neurvec0(cellinds);

    astro_vec=zeros(1,max(cellinds));
    astro_vec(astro_inds1)=1;
    astro_vec(eqz1)=0;
    astro_vec=astro_vec(cellinds);

    oligo_vec=zeros(1,max(cellinds));
    oligo_vec(oligo_inds1)=1;
    oligo_vec(eqz1)=0;
    oligo_vec=oligo_vec(cellinds);

    vlmc_vec=zeros(1,max(cellinds));
    vlmc_vec(vlmc_inds1)=1;
    vlmc_vec(eqz1)=0;
    vlmc_vec=vlmc_vec(cellinds);

    micro_vec=zeros(1,max(cellinds));
    micro_vec(micro_inds1)=1;
    micro_vec(eqz1)=0;
    micro_vec=micro_vec(cellinds);

    smc_vec=zeros(1,max(cellinds));
    smc_vec(smc_inds1)=1;
    smc_vec(eqz1)=0;
    smc_vec=smc_vec(cellinds);

    peri_vec=zeros(1,max(cellinds));
    peri_vec(peri_inds1)=1;
    peri_vec(eqz1)=0;
    peri_vec=peri_vec(cellinds);

    endo_vec=zeros(1,max(cellinds));
    endo_vec(endo_inds1)=1;
    endo_vec(eqz1)=0;
    endo_vec=endo_vec(cellinds);


    cr_vec=zeros(1,max(cellinds));
    cr_vec(CRexprInds)=1;
    cr_vec(eqz1)=0;
    cr_vec=cr_vec(cellinds);



    pv_vec=zeros(1,max(cellinds));
    pv_vec(PvalbexprInds)=1;
    pv_vec(eqz1)=0;
    pv_vec=pv_vec(cellinds);

    sst_vec=zeros(1,max(cellinds));
    sst_vec(SSTexprInds)=1;
    sst_vec(eqz1)=0;
    sst_vec=sst_vec(cellinds);
    end

    for count2=2:40

        inds=find(expression_matrix0(goodinds,count1,count2)>0);
        goodinds2=goodinds(inds);
leftside= scaleVec(goodinds2).*( expression_matrix(goodinds2,count1,count2));%- base_mat(count1,count2) .*baseline_Scale(goodinds2) );
  
        % leftside=leftside.*(leftside>0);
        leftside=[leftside;enforceVal.*exc_mat(count1,count2)];
        leftside=[leftside;enforceVal.*inh_mat(count1,count2)];

        leftside=[leftside;enforceVal.*sumAstro(count1,count2)];
        leftside=[leftside;enforceVal.*sumOligo(count1,count2)];
        leftside=[leftside;enforceVal.*sumVlmc(count1,count2)];
        leftside=[leftside;enforceVal.*sumMicro(count1,count2)];
        leftside=[leftside;enforceVal.*sumSmc(count1,count2)];
        leftside=[leftside;enforceVal.*sumPeri(count1,count2)];
        leftside=[leftside;enforceVal.*sumEndo(count1,count2)];

%         leftside=[leftside;enforceVal.*cr_mat(count1,count2)];
%         leftside=[leftside;enforceVal.*pv_mat(count1,count2)];
%         leftside=[leftside;enforceVal.*sst_mat(count1,count2)];

        %oligo_inds1,vlmc_inds1,micro_inds1,smc_inds1,peri_inds1,endo_inds1);

        %enforceVal.*sumL2(count1,count2)];
          rightside=ish_Scale(goodinds2).* clusters0(goodinds2,cellinds) .*repmat(scaleVec(goodinds2),1,length(cellinds));
       
        %  size(rightside)
        %  size(l2_vec)
        rightside=[rightside;enforceVal.*excvec1];
        rightside=[rightside;enforceVal.*inhvec1];

      rightside=[rightside;enforceVal.*astro_vec];
        rightside=[rightside;enforceVal.*oligo_vec];
        rightside=[rightside;enforceVal.*vlmc_vec];
        rightside=[rightside;enforceVal.*micro_vec];
        rightside=[rightside;enforceVal.*smc_vec];
        rightside=[rightside;enforceVal.*peri_vec];
        rightside=[rightside;enforceVal.*endo_vec];

        %         rightside=[rightside;enforceVal.*cr_vec];
%         rightside=[rightside;enforceVal.*pv_vec];
%         rightside=[rightside;enforceVal.*sst_vec];

        [vec1,resnorm,res1]=lsqnonneg(rightside,leftside,options);

        leftside= ( expression_matrix(goodinds2,count1,count2));%- base_mat(count1,count2) .*baseline_Scale(goodinds2) );
        % leftside=leftside.*(leftside>0);
        sol1=clusters0(goodinds2,cellinds)*vec1;

        errvec(goodinds2)= errvec(goodinds2)+transpose(min([transpose(sol1);transpose(leftside)]));
        %        errvec(goodinds)=errvec(goodinds)+abs(res1(1:length(goodinds)))./(max(scaleVec(goodinds),1e-12));
        error=error+resnorm;
        %       vec1=vec1.*(vec1>0);
        %allvec_individual(:,count1,count2)=vec1;
        allvec_individual(cellinds,count1,count2)=max(0,vec1);
        pred_genes(:,count1,count2)=(transpose(vec1)*transpose(ish_Scale(goodinds).*clusters0(goodinds,cellinds)));
        pred_genes(:,count1,count2)=pred_genes(:,count1,count2).*(expression_matrix(goodinds,count1,count2)>0);
        %pred_genes(:,count1,count2)=transpose(vec1)*transpose(clusters0(goodinds,cellinds));
        if(sum(isnan( pred_genes(:,count1,count2))))
            count1
            count2
            pause
        end
    end
end


%errvec(goodinds)= squeeze(sum(sum( min( expression_matrix(goodinds,:,:),pred_genes),2),3));
temp=expression_matrix(goodinds,:,:)-pred_genes;
%temp=expression_matrix(goodinds,:,:),length(goodinds),1,1) .*(repmat(reshape(baseline_Scale(goodinds),length(goodinds),1,1),1,100,40)) -pred_genes;
%size(temp)
errvec(goodinds)= squeeze(sum(sum(  temp.^2 ,2),3));


%errvec(goodinds)=errvec(goodinds)./scaleVec(goodinds)./(expr_base+1e-9);
%errvec(goodinds)=errvec(goodinds);%./expr_base;


%for count1=1:length(goodinds)
%    errvec(goodinds(count1))=sqrt(  errvec(goodinds(count1)))./sum(sum(expression_matrix(goodinds(count1),:,:) ));
%end
error=sqrt(sum( errvec(goodinds) ));
errvec(goodinds)=sqrt(errvec(goodinds));
ratiovec=ones(1,length(clusters0));
ishratiovec=ones(1,length(clusters0));
baseratiovec=zeros(1,length(clusters0));
if(0)
    for count1=1:length(goodinds)
        if sum(sum( pred_genes(count1,:,2:40) ))
            ratiovec(goodinds(count1))=reshape(squeeze(expression_matrix(goodinds(count1),:,2:40)),1,3900)/reshape(squeeze(pred_genes(count1,:,2:40)),1,3900);
        end
        leftside=reshape(squeeze(expression_matrix(goodinds(count1),:,2:40)),1,3900);
        rightside=[reshape(squeeze(pred_genes((count1),:,2:40)),1,3900)];%;reshape(squeeze(base_mat(:,2:40)),1,3900)];
        %baseratiovec(goodinds(count1))=max(0,reshape(squeeze(expression_matrix(goodinds(count1),:,2:40)-pred_genes((count1),:,2:40)),1,3900))/reshape(squeeze(baseline(:,2:40)),1,3900);
        a =(lsqnonneg(transpose(rightside),transpose(leftside)));

        if(a(2)>0)
            ishratiovec(goodinds(count1))=a(1);
            baseratiovec(goodinds(count1))=a(2);
        end
    end
end