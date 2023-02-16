if(1)

clear
close all
load('allsol1')
load('clusters0')
load('clusters')
load('ttype_labels')
load('goodinds')
load('expression_matrix_processed');
%load('clusters')
load('num_genes')
load('sumexpr0')
load('max_base')
load('base_mat')
%load expression_cv;
end

if(1)
base1=squeeze(sum(expression_matrix_processed(goodinds,10:100,2:40),1));
base1=base1./sum(sum(base1));
cv_exp=zeros(num_genes,1);
for count1=1:length(goodinds)
    target=squeeze(expression_matrix_processed(goodinds(count1),10:100,2:40));
    mult1=reshape(target,1,91*39)/reshape(transpose(base1),1,91*39);
    sum_vec=(squeeze(expression_matrix_processed(goodinds(count1),10:100,2:40))-mult1.*base1);
    %sum_vec=sum_vec(7:100);
    cv_exp(goodinds(count1))=sqrt(sum(sum(sum_vec.*sum_vec)));
end
cv_exp(goodinds)=cv_exp(goodinds)./max(cv_exp);
end

m1=sum(transpose(clusters0(goodinds,:))>0);
goodinds(find(m1==0))=[];

make_inds;

[a b ]=size(clusters);


clusters1=0.*clusters;
clusters1(goodinds,(localinds2))=clusters(goodinds,(localinds2));
clusters=clusters1;

num_ttypes=numCells;
clusters00=clusters;
clusters=0.*clusters;
clusters(goodinds,:)=clusters00(goodinds,:);

ish_Scale=ones(num_genes,1);
scaleVec=ones(num_genes,1);
scaleVec0=ones(num_genes,1);
baseline_Scale=zeros(num_genes,1);
bestScaleVec5=ones(num_genes,1);


cv_genes=transpose(std(transpose(clusters(:,localinds2))));


m1=(sum(clusters(goodinds,:))>0);
localinds=unique(localinds2(find(m1(localinds2))));
localinds2=localinds;

get_means;

means=transpose(mean(transpose(clusters(:,localinds2))));
ss=(clusters(goodinds,localinds)-repmat(means(goodinds),1,length(localinds2))).^2;
ss=max(0,(clusters(goodinds,localinds)-repmat(means(goodinds),1,length(localinds2)))).^2;


ish_Scale=1+0.*ish_Scale;
goodinds2=goodinds;

%%%%%%%%%

disp 'expression matrix processing'

if(1)
    ave_transcriptome=zeros(length(goodinds),100);
    for count1=1:100
        for count2=1:length(goodinds)
            vec1=expression_matrix_processed(goodinds(count2),count1,:);
            inds2=find(vec1>0);
            if(length(inds2))
                ave_transcriptome(count2,count1)=mean(vec1(inds2));
            end
        end
    end
    gene_expression_matrix=[];
    for count1=1:100
        for count2=2:40
            vec1=squeeze(transpose(squeeze(expression_matrix_processed(goodinds,count1,count2) )));
            inds=find(vec1>0);
            if(length(inds)==0)
                continue
            end
            val=vec1(inds)/transpose(ave_transcriptome(inds,count1));
            vec2=vec1;
            % vec2=val .*transpose(ave_transcriptome(:,count1));
            % vec2(inds)=vec1(inds);
            gene_expression_matrix=[gene_expression_matrix;vec2];
            % expression_matrix(goodinds,count1,count2) =vec2;
        end
    end
end


if(1)
expression_matrix_processed2=expression_matrix_processed;
    for count1=1:length(goodinds)
        %count1
        temp=squeeze(expression_matrix_processed(goodinds(count1),:,:));

         sum1=sum(sum(sum( expression_matrix_processed(goodinds(count1),:,2:40) )));
%         sumexpr0(goodinds(count1))=sum1;
%        sum1=sumexpr0(goodinds(count1));
        expression_matrix_processed2(goodinds(count1),:,:)=expression_matrix_processed(goodinds(count1),:,:)./sum1;
    end
end
%sumexpr0=squeeze(sum(sum(expression_matrix_processed2(:,:,:),2),3));

expression_matrix_processed2(:,:,1)=0;

ngtz=sum(gene_expression_matrix>0,1);
clusters1=clusters;



%std_clusters0(std_clusters0==0)=2*max(max(std_clusters0))/100;

%[a b]=size(std_clusters0);
%vec0=transpose(std(transpose(clusters0)));
%std_clusters0=std_clusters0.*repmat(vec0,1,b);

%what percent of expression matrix is occupied?
goodExpr=[];
for count1=1:length(goodinds)
    goodExpr=[goodExpr;sum(sum(sum( expression_matrix_processed(goodinds(count1),:,:)>0 )))];
end

gtz1=sum(clusters0(goodinds,:)>0);
goodClust=[];
for count1=1:length(goodinds)
    inds=find(clusters00(goodinds(count1),:)>0);  
    val=0;
    if(length(inds))
%        val=sum(gtz1(inds));
         val=min(gtz1(inds));
    end
    goodClust=[goodClust;val];
end

if(0)
meanL4=transpose(mean(transpose(clusters(:,l4_inds1))));
meanL2=transpose(mean(transpose(clusters(:,l2_inds1))));
meanL3=transpose(mean(transpose(clusters(:,l3_inds1))));
meanL23=transpose(mean(transpose(clusters(:,l23_inds1))));
meanL5=transpose(mean(transpose(clusters(:,l5_inds1))));
meanL6=transpose(mean(transpose(clusters(:,l6_inds1))));

meanInh=transpose(mean(transpose(clusters(:,inh_inds1))));
meanExc=transpose(mean(transpose(clusters(:,exc_inds1))));
meanAll=transpose(mean(transpose(clusters(:,localinds2))));
end

%scaleVec2=sumexpr0;
scaleVec2=sumexpr0;

geneScale=1./sumexpr0;
geneScale0=geneScale;

%%%%%%%%%%%%%
m1=[clusters(goodinds,localinds2)];
%sol1=(transpose(sumexpr0(goodinds))/transpose(m1));
%sol1=lsqnonneg(transpose(m1),transpose(sumexpr0(goodinds)));
%sol1=transpose(lsqnonneg(m1,sumexpr0(goodinds)));
%plot(mean(transpose(clusters(goodinds,localinds2))),sumexpr0(goodinds),'.')
sol1=allsol1(localinds2);
figure
hold on
maxx=0;
maxy=0;
for count1=1:length(goodinds)
    %    sum1=sum(sum(sum(expression_matrix_processed2(goodinds(count1),:,:))));
    sum1=sumexpr0(goodinds(count1));

    inds2=find(clusters(goodinds(count1),localinds2));
    %  plot(mean(transpose(clusters(goodinds(count1),localinds2(inds2)))),sumexpr0(goodinds(count1)),'.')
    predval=sum(clusters(goodinds(count1),localinds2) .*sol1(1:(length(sol1))));
    plot(predval,sum1,'.')

    %  plot(predval,sumexpr0(goodinds(count1)),'.')
    maxx=max(predval,maxx);
    maxy=max(sumexpr0(goodinds(count1)),maxy);
end
xlim([0 maxx])
ylim([0 maxy])
plot([0 maxx],[0 maxx],'k')
sum2=zeros(1,400);
sum2(localinds2)=sol1(1:length(sol1));
%sum(sum2(inh_inds1))/sum(sum2(neur_inds1))

%%%%%%%%%%%%%
disp 'solving'

clusters1=clusters;

%scaleVec2=ones(num_genes,1);
% scaleVec2=sumexpr0;
% scaleVec2=scaleVec2.*(transpose(sum(transpose(clusters)>0)));
%scaleVec2=sumexpr0;
% scaleVec2(goodinds)=scaleVec2(goodinds).*goodClust(goodinds);
% scaleVec2(goodinds)=scaleVec2(goodinds).*goodExpr;
scaleVec2=scaleVec2./max(scaleVec2);
scaleVec3=scaleVec2;
scaleVec3(goodinds)=scaleVec2(goodinds);%.*vec1;
baseline_Scale5=baseline_Scale;%50.*transpose(baseratiovec5);

clusters2=clusters;
if(1)
    clusters2=consolidate(clusters2,clusters0,goodinds,l2_inds1);
    clusters2=consolidate(clusters2,clusters0,goodinds,l23_inds1);
    clusters2=consolidate(clusters2,clusters0,goodinds,l4_inds1);
    clusters2=consolidate(clusters2,clusters0,goodinds,l45_inds1);
    clusters2=consolidate(clusters2,clusters0,goodinds,l5_inds1);
    clusters2=consolidate(clusters2,clusters0,goodinds,l6_inds1);
    clusters2=consolidate(clusters2,clusters0,goodinds,l6b_inds1);

    %  clusters2=consolidate(clusters2,clusters0,goodinds,inh_inds1);
    %  clusters2=consolidate(clusters2,clusters0,goodinds,neur_inds1);

    idx=kmeans(transpose(clusters(goodinds,inh_inds1)),10);
for count6=1:10
    clusters2=consolidate(clusters2,clusters0,goodinds,inh_inds2(find(idx==count6)));
end
  
    if(0)
    clusters2=consolidate(clusters2,clusters0,goodinds,CRexprInds);
    clusters2=consolidate(clusters2,clusters0,goodinds,PvalbexprInds);
    clusters2=consolidate(clusters2,clusters0,goodinds,SSTexprInds);
    clusters2=consolidate(clusters2,clusters0,goodinds,other_inh_inds1);
end

    %   clusters2=consolidate(clusters2,clusters0,goodinds,astro_inds1);
    %     clusters2=consolidate(clusters2,clusters0,goodinds,oligo_inds1);
    %     clusters2=consolidate(clusters2,clusters0,goodinds,smc_inds1);
    %     clusters2=consolidate(clusters2,clusters0,goodinds,vlmc_inds1);
    %     clusters2=consolidate(clusters2,clusters0,goodinds,micro_inds1);
    %     clusters2=consolidate(clusters2,clusters0,goodinds,endo_inds1);
    %     clusters2=consolidate(clusters2,clusters0,goodinds,peri_inds1);
end

scaled_clusters_fixed=clusters2.*(clusters2>0);
for count1=1:length(goodinds)
    scaled_clusters_fixed(goodinds(count1),:)=scaled_clusters_fixed(goodinds(count1),:).*geneScale(goodinds(count1));%.*ratlist(goodinds(count1));%.*min(2.63,medlist((count1)));
end

[allvec_individual6,errvec6,error,ratiovec6,ishratiovec5,baseratiovec6]= calc_err26nnn(scaleVec2,exc_inds1,goodinds,expression_matrix_processed2,expression_matrix_processed2,baseline_Scale5,ish_Scale,scaled_clusters_fixed,num_ttypes,localinds,base_mat);
allvec_individual_vector5=1e7.* squeeze(sum(allvec_individual6(:,:,1:40),3));
allvec_individual5=allvec_individual6;

neur_mat=sum(allvec_individual6(neur_inds1,:,:),1);

inh_mat=sum(allvec_individual6(inh_inds1,:,:),1);
exc_mat=sum(allvec_individual6(exc_inds1,:,:),1);

cr_mat=sum(allvec_individual6(CRexprInds,:,:),1);
sst_mat=sum(allvec_individual6(SSTexprInds,:,:),1);
pv_mat=sum(allvec_individual6(PvalbexprInds,:,:),1);
astro_mat=squeeze(sum(allvec_individual6(astro_inds1,:,:),1));
oligo_mat=squeeze(sum(allvec_individual6(oligo_inds1,:,:),1));
vlmc_mat=squeeze(sum(allvec_individual6(vlmc_inds1,:,:),1));
micro_mat=squeeze(sum(allvec_individual6(micro_inds1,:,:),1));
smc_mat=squeeze(sum(allvec_individual6(smc_inds1,:,:),1));
peri_mat=squeeze(sum(allvec_individual6(peri_inds1,:,:),1));
endo_mat=squeeze(sum(allvec_individual6(endo_inds1,:,:),1));
other_inds3=[astro_inds1;oligo_inds1;vlmc_inds1;micro_inds1;smc_inds1;peri_inds1;endo_inds1];
other_mat=squeeze(sum(allvec_individual6(other_inds3,:,:),1));
cr_mat=sum(allvec_individual6(CRexprInds,:,:),1);
sst_mat=sum(allvec_individual6(SSTexprInds,:,:),1);
pv_mat=sum(allvec_individual6(PvalbexprInds,:,:),1);
%sum(sum(sum(oligo_mat)))/sum(sum(sum(astro_mat)))


figure
plot(sum(transpose(squeeze(exc_mat))))
hold on
plot(sum(transpose(squeeze(inh_mat))),'r')
%sum(sum(squeeze(inh_mat)))/sum(sum(squeeze(exc_mat+inh_mat)))

clusters5=clusters;
scaled_clusters_fixed=clusters5.*(clusters5>0);
for count1=1:length(goodinds)
    scaled_clusters_fixed(goodinds(count1),:)=scaled_clusters_fixed(goodinds(count1),:).*geneScale(goodinds(count1));%.*ratlist(goodinds(count1));%.*min(2.63,medlist((count1)));
end
inds1=find(sum(clusters5(goodinds,:)));
[allvec_individual6,errvec5,error,ratiovec5,ishratiovec5,baseratiovec5]= calc_err26nnn_constrained10c(scaleVec2,squeeze(exc_mat),exc_inds1,goodinds,expression_matrix_processed2,expression_matrix_processed2,baseline_Scale5,ish_Scale,scaled_clusters_fixed,num_ttypes,inds1,base_mat,squeeze(inh_mat),squeeze(cr_mat),squeeze(pv_mat),squeeze(sst_mat),inh_inds1,astro_mat,oligo_mat,vlmc_mat,micro_mat,smc_mat,peri_mat,endo_mat,l2_inds1,l23_inds1,l3_inds1,l4_inds1,l45_inds1,l5_inds1,l6_inds1,l6b_inds1,astro_inds1,oligo_inds1,vlmc_inds1,micro_inds1,smc_inds1,peri_inds1,endo_inds1,CRexprInds,PvalbexprInds,SSTexprInds,neur_inds1,squeeze(neur_mat));

sum1=mean(sum(sum(allvec_individual6(neur_inds1,10:100,1:40),3),1),2);
allvec_individual_vector5=110/sum1.* squeeze(sum(allvec_individual6(:,:,1:40),3));
plotall1

figure
m1=squeeze(sum(allvec_individual6(l2_inds1,:,1:40),1));
m1= m1./max(max(m1));
plot_sol(m1)

figure
m1=squeeze(sum(allvec_individual6(l23_inds1,:,1:40),1));
m1= m1./max(max(m1));
plot_sol(m1)

figure
m1=squeeze(sum(allvec_individual6(l3_inds1,:,1:40),1));
m1= m1./max(max(m1));
plot_sol(m1)

figure
m1=squeeze(sum(allvec_individual6(l4_inds1,:,1:40),1));
m1= m1./max(max(m1));
plot_sol(m1)

figure
m1=squeeze(sum(allvec_individual6(l45_inds1,:,1:40),1));
m1= m1./max(max(m1));
plot_sol(m1)

figure
m1=squeeze(sum(allvec_individual6(l5_inds1,:,1:40),1));
m1= m1./max(max(m1));
plot_sol(m1)

figure
m1=squeeze(sum(allvec_individual6(l6_inds1,:,1:40),1));
m1= m1./max(max(m1));
plot_sol(m1)

figure
m1=squeeze(sum(allvec_individual6(l6b_inds1,:,1:40),1));
m1= m1./max(max(m1));
plot_sol(m1)

figure
m1=squeeze(sum(allvec_individual6(l2_inds1,:,1:40),1));
m1= m1./max(max(m1));
plot_sol(m1)

figure
m1=squeeze(sum(allvec_individual6(l2_inds1,:,1:40),1));
m1= m1./max(max(m1));
plot_sol(m1)

other_inds2=[micro_inds1;astro_inds1;smc_inds1;peri_inds1;oligo_inds1;endo_inds1;vlmc_inds1];
tot_other=sum(sum(allvec_individual_vector5(other_inds2,7:100)));
oligoastro=sum(sum(allvec_individual_vector5(oligo_inds1,7:100)))./sum(sum(allvec_individual_vector5(astro_inds1,7:100)));
tot_neur=sum(sum(allvec_individual_vector5(neur_inds1,7:100)));
if(0)
%https://reader.elsevier.com/reader/sd/pii/S0092867422001349?token=971BFA10F5F2B10F0E5289D10A53A5DA84EBF8CBD143E73E342F2BC4E91029C15EC0DF308EFD9403C5882C085F0E14AA&originRegion=eu-west-1&originCreation=20221114130024
%lit_total=71+44+31+15+2+1+5;
lit_neur=417+12+4+2+1+1+14;
endoprop=sum(sum(allvec_individual_vector5(endo_inds1,7:100)))./tot_neur
71/lit_neur
%71/lit_total
astroprop=sum(sum(allvec_individual_vector5(astro_inds1,7:100)))./tot_neur
28.7/60
%(28.7/60+44/lit_neur)/2
%44/lit_neur
microprop=sum(sum(allvec_individual_vector5(micro_inds1,7:100)))./tot_neur
4.5/60
%(4.5/60+31/lit_neur)/2
%31/lit_neur
oligoprop=sum(sum(allvec_individual_vector5(oligo_inds1,7:100)))./tot_neur
14.7/60
%(14.7/60+17/lit_neur)/2
%17/lit_neur
periprop=sum(sum(allvec_individual_vector5(peri_inds1,7:100)))./tot_neur
%1/lit_total
1/lit_neur
end