
function [allvec_individual,errvec,error,ratiovec,ishratiovec,baseratiovec]= calc_err23(scaleVec,exc_inds1,goodinds,expression_matrix,expression_matrix0,baseline_Scale,ish_Scale,clusters0,num_ttypes,placed_types,baseline);
[foo numcells]=size(clusters0);

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
options = optimset('Display','notify','TolX',1e-10);
%options = optimset('Display','notify','TolX',1e-7);

error=0;
%goodCellInds=find(sum( clusters0(goodinds(length(goodinds)-3),:)));
pred_genes=zeros(length(goodinds),100,40);
for count1=1:100
    cellinds=placed_types;
    if(count1<=6)
        cellinds=placed_types_no_exc;
    end
    for count2=2:40
        inds=find(expression_matrix0(goodinds,count1,count2)>0);
        goodinds2=goodinds(inds);
leftside= scaleVec(goodinds2).*( expression_matrix(goodinds2,count1,count2));
  %      leftside= scaleVec(goodinds2).*( expression_matrix(goodinds2,count1,count2)- baseline(count1,count2) .*baseline_Scale(goodinds2) );
         %leftside=leftside.*(leftside>0);
%          size(ish_Scale(goodinds2))
%          size( clusters0(goodinds2,cellinds) )
%          size(repmat(scaleVec(goodinds2),1,length(cellinds)))
%       rightside=ish_Scale(goodinds2).* clusters0(goodinds2,cellinds) .*repmat(scaleVec(goodinds2),1,length(cellinds));
%           dedd
        rightside=ish_Scale(goodinds2).* clusters0(goodinds2,cellinds) .*repmat(scaleVec(goodinds2),1,length(cellinds));
        [vec1,resnorm,res1]=lsqnonneg(rightside,leftside,options);
   leftside= ( expression_matrix(goodinds2,count1,count2));
       % leftside= ( expression_matrix(goodinds2,count1,count2)- baseline(count1,count2) .*baseline_Scale(goodinds2) );
        %leftside=leftside.*(leftside>0);
        sol1=clusters0(goodinds2,cellinds)*vec1;

        errvec(goodinds2)= errvec(goodinds2)+transpose(min([transpose(sol1);transpose(leftside)]));
        %        errvec(goodinds)=errvec(goodinds)+abs(res1(1:length(goodinds)))./(max(scaleVec(goodinds),1e-12));
        error=error+resnorm;
        %       vec1=vec1.*(vec1>0);
        %allvec_individual(:,count1,count2)=vec1;
        
        allvec_individual(cellinds,count1,count2)=vec1;
        if(sum(sum(isnan(ish_Scale(goodinds2).* clusters0(goodinds2,cellinds) .*repmat(scaleVec(goodinds2),1,length(cellinds))))))
       (sum(sum(isnan(clusters0(goodinds2,cellinds)))))
        
          disp 'isnan'
          pause
        end
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

if(0)
pred_genes=0.* pred_genes;
[a b c]=size(allvec_individual);
err1=0;
for count1=1:a
    vec1=squeeze(sum(allvec_individual(count1,:,:),2));
    for count2=1:100
        vec2=squeeze(allvec_individual(count1,count2,:));
        if (sum(vec2)==0)
            continue
        end
        coeff=transpose(vec2)/transpose(vec1);
        bestvec=coeff*vec1;        
        allvec_individual(count1,count2,:)=bestvec;

        %err1=err1+sum( (bestvec-vec2).^2 );
        for count3=2:40
            pred_genes(:,count2,count3)=pred_genes(:,count2,count3)+(bestvec(count3))*squeeze(clusters0(goodinds,(count1)));
        end
    end
end
end

%errvec(goodinds)= squeeze(sum(sum( min( expression_matrix(goodinds,:,:),pred_genes),2),3));
%temp=expression_matrix(goodinds,:,:)-pred_genes;
temp=expression_matrix(goodinds,:,:)- repmat(reshape(baseline,1,100,40),length(goodinds),1,1) .*(repmat(reshape(baseline_Scale(goodinds),length(goodinds),1,1),1,100,40)) -pred_genes;
%size(temp)
errvec(goodinds)= squeeze(sum(sum(  temp.^2 ,2),3));


%errvec(goodinds)=errvec(goodinds)./scaleVec(goodinds)./(expr_base+1e-9);
%errvec(goodinds)=errvec(goodinds);%./expr_base;


%for count1=1:length(goodinds)
%    errvec(goodinds(count1))=sqrt(  errvec(goodinds(count1)))./sum(sum(expression_matrix(goodinds(count1),:,:) ));
%end
error=sqrt(sum( errvec(goodinds) ));
%errvec(goodinds)=sqrt(errvec(goodinds));
ratiovec=ones(1,length(clusters0));
ishratiovec=ones(1,length(clusters0));
baseratiovec=zeros(1,length(clusters0));
if(0)
    for count1=1:length(goodinds)
        if sum(sum( pred_genes(count1,:,2:40) ))
            ratiovec(goodinds(count1))=reshape(squeeze(expression_matrix(goodinds(count1),:,2:40)),1,3900)/reshape(squeeze(pred_genes(count1,:,2:40)),1,3900);
        end
        leftside=reshape(squeeze(expression_matrix(goodinds(count1),:,2:40)),1,3900);
        rightside=[reshape(squeeze(pred_genes((count1),:,2:40)),1,3900);reshape(squeeze(baseline(:,2:40)),1,3900)];
        %baseratiovec(goodinds(count1))=max(0,reshape(squeeze(expression_matrix(goodinds(count1),:,2:40)-pred_genes((count1),:,2:40)),1,3900))/reshape(squeeze(baseline(:,2:40)),1,3900);
       a =(lsqnonneg(transpose(rightside),transpose(leftside)));
     
       if(a(2)>0)
        ishratiovec(goodinds(count1))=a(1);
        baseratiovec(goodinds(count1))=a(2);
       end
    end
end