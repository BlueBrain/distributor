%perform the normization of the data distributions 

clear
load('dot_expression_matrix')
load('expression_matrix')
load ('goodinds')
load ('clusters0')
load ('gtot2')
clusters=clusters0;

if(exist('map.mat'))
    load('map')
else
    map=rand(200,3);
    save('map','map')
end

m1=squeeze(mean(expression_matrix(goodinds,:,:),1));
m0=squeeze(mean(dot_expression_matrix(goodinds,:,:),1));
meanRatio=mean(mean(m1./(m0+1e-9)));
%normalize to dot expression matrix
expression_matrix_processed=expression_matrix;

figure
hold on
count0=1;
for count2=40:-1:2
    
    x1=map(count0,:);
    count0=count0+1;
    for count1=1:100
        vec1=expression_matrix_processed(goodinds,count1,count2)./dot_expression_matrix(goodinds,count1,count2);
        vec2=vec1;
        testsum=sum(vec1>0);
        vec2(isnan(vec2))=0;
        %        plot([0;sort(vec2)],'color',x1)
        plot(count2.*ones(length(vec2),1),vec2,'k.')
    end

end


ratio_matrix=zeros(100,40);
for count1=1:100
    for count2=1:40

        vec1=expression_matrix_processed(goodinds,count1,count2)./dot_expression_matrix(goodinds,count1,count2);
        inds=find(vec1>0);
        if(length(inds))
            median1=median(vec1(inds));
            ratio_matrix(count1,count2)=mean(vec1(inds))/median1;
            expression_matrix_processed(goodinds,count1,count2)=expression_matrix_processed(goodinds,count1,count2)./median1;
        end
    end
end

figure
plot(transpose(ratio_matrix),'b.')
hold on
plot(1:40,ones(1,40),'k.')

expression_matrix(goodinds,:,1)=0;
dot_expression_matrix(goodinds,:,1)=0;


m1=squeeze(mean(expression_matrix(goodinds,:,:),1));
m0=squeeze(mean(dot_expression_matrix(goodinds,:,:),1));
meanRatio=mean(mean(m1./(m0+1e-9)));


%normalize to dot expression matrix
figure
hold on
count0=1;
for count2=40:-1:2
    %    x1=rand(1,3);
    x1=map(count0,:);
    count0=count0+1;
    for count1=1:100
        vec1=expression_matrix_processed(goodinds,count1,count2)./dot_expression_matrix(goodinds,count1,count2);
        vec2=vec1;%(vec1>0);
        testsum=sum(vec1>0);
        vec2(isnan(vec2))=0;
        %        plot([0;sort(vec2)],'color',x1)
        plot(count2.*ones(length(vec2),1),vec2,'k.')

    end

end

localinds2=find(sum((clusters0(goodinds,:))>0));
scaleVec=transpose(sum(transpose(clusters)>0));

scaleVec=scaleVec./max(scaleVec);

sumexpr0=squeeze(sum(sum(expression_matrix_processed,3),2));
save('sumexpr0','sumexpr0')

base_mat=sum(expression_matrix_processed(goodinds,:,:),1);
base_mat=base_mat./sum(sum(sum(base_mat)));

if(1)
    max_base=[];
    for count1=1:length(goodinds)
        count1
        val=sumexpr0(goodinds(count1));
        while(1)
            test=expression_matrix_processed(goodinds(count1),:,:)-val.*base_mat;
            pos1=sum(sum(max(test,0)));
            neg1=sum(sum(max(-test,0)));
            rat1=neg1/pos1;
            if (rat1>0.1)
                val=val.*0.95;
            else
                break
            end
        end
        max_base=[max_base;val];
    end
    plot(sort(max_base))
end
max_base=0.*max_base;

make_inds

[a b]=size(clusters);
rank_clusters=zeros(length(goodinds),b);
for count1=1:b
    [sort1,i1]=sort(clusters(goodinds,(count1)));
    rank_clusters(i1,(count1))=1:length(i1);
    inds=find(clusters(goodinds,(count1))==0);
    rank_clusters(inds,(count1))=0;
end
[astro_rank] = get_medians(goodinds,rank_clusters,astro_inds1);
[oligo_rank] = get_medians(goodinds,rank_clusters,oligo_inds1);
[vlmc_rank] = get_medians(goodinds,rank_clusters,vlmc_inds1);
[smc_rank] = get_medians(goodinds,rank_clusters,smc_inds1);
[peri_rank] = get_medians(goodinds,rank_clusters,peri_inds1);
[micro_rank] = get_medians(goodinds,rank_clusters,micro_inds1);
[endo_rank] = get_medians(goodinds,rank_clusters,endo_inds1);
[l2_rank] = get_medians(goodinds,rank_clusters,l2_inds1);
[l23_rank] = get_medians(goodinds,rank_clusters,l23_inds1);
[l3_rank] = get_medians(goodinds,rank_clusters,l3_inds1);
[l4_rank] = get_medians(goodinds,rank_clusters,l4_inds1);
[l45_rank] = get_medians(goodinds,rank_clusters,l45_inds1);
[l5_rank] = get_medians(goodinds,rank_clusters,l5_inds1);
[l6_rank] = get_medians(goodinds,rank_clusters,l6_inds1);
[l6b_rank] = get_medians(goodinds,rank_clusters,l6b_inds1);
[exc_rank] = get_medians(goodinds,rank_clusters,exc_inds1);
[inh_rank] = get_medians(goodinds,rank_clusters,inh_inds1);
[neur_rank] = get_medians(goodinds,rank_clusters,neur_inds1);
[cell_rank] = get_medians(goodinds,rank_clusters,localinds2);
[all_cell_rank] = get_medians(goodinds,rank_clusters,1:b);

[pvalb_rank] = get_medians(goodinds,rank_clusters,pvalb_inds1);
[sst_rank] = get_medians(goodinds,rank_clusters,sst_inds1);
[vip_rank] = get_medians(goodinds,rank_clusters,vip_inds1);
[lamp_rank] = get_medians(goodinds,rank_clusters,lamp_inds1);
[sncg_rank] = get_medians(goodinds,rank_clusters,sncg_inds1);

[rank_clusters] = fill_in_zeros(rank_clusters,astro_inds1,astro_rank);
[rank_clusters] = fill_in_zeros(rank_clusters,oligo_inds1,oligo_rank);
[rank_clusters] = fill_in_zeros(rank_clusters,vlmc_inds1,vlmc_rank);
[rank_clusters] = fill_in_zeros(rank_clusters,smc_inds1,smc_rank);
[rank_clusters] = fill_in_zeros(rank_clusters,peri_inds1,peri_rank);
[rank_clusters] = fill_in_zeros(rank_clusters,micro_inds1,micro_rank);
[rank_clusters] = fill_in_zeros(rank_clusters,endo_inds1,endo_rank);
[rank_clusters] = fill_in_zeros(rank_clusters,l2_inds1,l2_rank);
[rank_clusters] = fill_in_zeros(rank_clusters,l23_inds1,l23_rank);
[rank_clusters] = fill_in_zeros(rank_clusters,l3_inds1,l3_rank);
[rank_clusters] = fill_in_zeros(rank_clusters,l4_inds1,l4_rank);
[rank_clusters] = fill_in_zeros(rank_clusters,l45_inds1,l45_rank);
[rank_clusters] = fill_in_zeros(rank_clusters,l5_inds1,l5_rank);
[rank_clusters] = fill_in_zeros(rank_clusters,l6_inds1,l6_rank);
[rank_clusters] = fill_in_zeros(rank_clusters,l6b_inds1,l6b_rank);
[rank_clusters] = fill_in_zeros(rank_clusters,exc_inds1,exc_rank);

[rank_clusters] = fill_in_zeros(rank_clusters,pvalb_inds1,pvalb_rank);
[rank_clusters] = fill_in_zeros(rank_clusters,sst_inds1,sst_rank);
[rank_clusters] = fill_in_zeros(rank_clusters,vip_inds1,vip_rank);
[rank_clusters] = fill_in_zeros(rank_clusters,lamp_inds1,lamp_rank);
[rank_clusters] = fill_in_zeros(rank_clusters,sncg_inds1,sncg_rank);
[rank_clusters] = fill_in_zeros(rank_clusters,inh_inds1,inh_rank);

[rank_clusters] = fill_in_zeros(rank_clusters,neur_inds1,neur_rank);
[rank_clusters] = fill_in_zeros(rank_clusters,localinds2,cell_rank);
[rank_clusters] = fill_in_zeros(rank_clusters,1:b,all_cell_rank);

best_rank_clusters=rank_clusters;
minindlength=1e9;
for count1=1:length(localinds2)
    inds=find(clusters(goodinds,localinds2(count1))==0);
    if (length(inds)<minindlength)
        minindlength=length(inds);
        refdist=sort(clusters(goodinds,localinds2(count1)));
        refind=count1;
    end
end

refdist(1:minindlength)=max(refdist)./1000.*[1:minindlength]./minindlength;

refdist=transpose([1:length(goodinds)]./length(goodinds));

clusters0=clusters;

for count3=1:(length(localinds2))
    clusters0(goodinds,localinds2(count3))=refdist(best_rank_clusters(:,localinds2(count3)));
end

figure
hold on
for count1=1:length(localinds2)
    plot(sort(clusters(goodinds,localinds2(count1))))
end

figure
hold on
for count1=1:length(localinds2)
     plot(sort(clusters0(goodinds,localinds2(count1))))
end

best_clusters=clusters0;
bestrefdist=refdist;

options = optimset('Display','notify','TolX',1e-10);
sol1=transpose(lsqnonneg(repmat(scaleVec(goodinds),1,length(localinds2)).*best_clusters(goodinds,localinds2),scaleVec(goodinds).*sumexpr0(goodinds),options));
figure
hold on
maxx=0;
maxy=0;
p1=[];
s1=[];
for count1=1:length(goodinds)
    %    sum1=sum(sum(sum(expression_matrix_processed(goodinds(count1),:,:))));
    sum1=sumexpr0(goodinds(count1));

    inds2=find(clusters0(goodinds(count1),localinds2));
    %  plot(mean(transpose(clusters0(goodinds(count1),localinds2(inds2)))),sumexpr(goodinds(count1)),'.')
    predval=sum(best_clusters(goodinds(count1),localinds2) .*sol1(1:(length(sol1))));
    plot(predval,sum1,'.')

    p1=[p1,predval];
    s1=[s1,sum1];
    %  plot(predval,sumexpr(goodinds(count1)),'.')
    maxx=max(predval,maxx);
    maxy=max(sumexpr0(goodinds(count1)),maxy);
end
xlim([0 maxx])
ylim([0 maxy])
plot([0 maxx],[0 maxx],'k')

converge_queue=ones(5,1);
constant1=1/1000; %to scale random numbers
bad_num=0;

for countY=1:10000 %number of cycles to test
    currrefdist=bestrefdist;
    for count0=1:10000
        %count0
        rank_clusters=best_rank_clusters;

        test_clusters=best_clusters;

        testrefdist=bestrefdist;
        if(count0==1)
            besterr=1e9;
        else

            %  testrefdist=sort(max(0,testrefdist+max(testrefdist).*constant1.*randn(length(testrefdist),1)));
            testrefdist=sort((testrefdist+max(testrefdist).*constant1.*randn(length(testrefdist),1)));
            testrefdist=testrefdist-min(testrefdist);
            testrefdist=testrefdist./max(testrefdist);

            for count3=1:(length(localinds2))
                test_clusters(goodinds,localinds2(count3))=testrefdist(rank_clusters(:,localinds2(count3)));
            end


        end

        m1=[test_clusters(goodinds,localinds2)];

        sol1=transpose(lsqnonneg(repmat(scaleVec(goodinds),1,length(localinds2)).*m1,scaleVec(goodinds).*sumexpr0(goodinds)));

        predvec1=(test_clusters(goodinds,localinds2) *transpose(sol1));
        testerr=sum( scaleVec(goodinds) .* (predvec1-(sumexpr0(goodinds))).^2);

        if(testerr<besterr)
            bad_num=0;
            if (besterr<1e9)
                testerr
            end
            besterr=testerr;
            best_rank_clusters=rank_clusters;
            bestrefdist=testrefdist;
            best_clusters=test_clusters;
            converge_queue=[converge_queue;besterr];
            converge_queue(1)=[];
        else
            bad_num=bad_num+1;
            if (bad_num>100)
                constant1=constant1*0.98
                bad_num=0;
            end
        end

    end

    if(abs(converge_queue(1)-converge_queue(5))/converge_queue(1)<1e-6)
        disp 'Converged'
        save dheck1
        break
    end

    finalrefdist=bestrefdist;
    for multiplier=0.1:0.2:120
        testrefdist=sort(finalrefdist+multiplier *(finalrefdist-currrefdist));
        testrefdist=testrefdist./max(testrefdist);
        testrefdist=min(1,testrefdist);
        testrefdist=max(0,testrefdist);
        for count3=1:(length(localinds2))
            test_clusters(goodinds,localinds2(count3))=testrefdist(best_rank_clusters(:,localinds2(count3)));
        end
        m1=[test_clusters(goodinds,localinds2)];
        %sol1=transpose(lsqnonneg((m1),(sumexpr0(goodinds))));
        sol1=transpose(lsqnonneg(repmat(scaleVec(goodinds),1,length(localinds2)).*m1,scaleVec(goodinds).*sumexpr0(goodinds)));

        predvec1=(test_clusters(goodinds,localinds2) *transpose(sol1));
        testerr=sum(scaleVec(goodinds) .*(predvec1-(sumexpr0(goodinds))).^2);
        if(testerr<besterr)
            multiplier
            disp 'assigning'
            besterr=testerr
            best_rank_clusters=rank_clusters;
            bestrefdist=testrefdist;
            best_clusters=test_clusters;
        else
            break
        end

    end
    
    hold on
    plot(bestrefdist)
    drawnow
    save dheck1
end

m1=[best_clusters(goodinds,localinds2)];
sol1=transpose(lsqnonneg(repmat(scaleVec(goodinds),1,length(localinds2)).*m1,scaleVec(goodinds).*sumexpr0(goodinds)));

for count0=1:33333
    testrefdist=bestrefdist;
    testrefdist(end)=0.9999.*testrefdist(end);
    %    testrefdist(1)=1.001.*testrefdist(1);
    testrefdist=testrefdist./max(testrefdist);
    testrefdist=sort(testrefdist);
    predvec1=(best_clusters(goodinds,localinds2) *transpose(sol1));
    testerr=sum(scaleVec(goodinds) .*(predvec1-(sumexpr0(goodinds))).^2);
    if(testerr<besterr)
        [count1 testerr]
        % pause
        besterr=testerr;
        best_rank_clusters=rank_clusters;
        bestrefdist=testrefdist;
        best_clusters=test_clusters;
    else
        break
    end
end

[a b]=size(clusters0);
optim_table=zeros(b,length(goodinds));
m1=[best_clusters(goodinds,localinds2)];
for count1=1:100
    count1
    for count2=2:40
        if (sum(squeeze(expression_matrix_processed(goodinds,count1,count2))))
            sol1=transpose(lsqnonneg((m1),squeeze(expression_matrix_processed(goodinds,count1,count2)),options));
            inds=(find(sol1>0));
            % length(inds)
            for count3=1:length(inds)
                test=optim_table(localinds2(inds(count3)),:)+transpose(squeeze(expression_matrix_processed(goodinds,count1,count2)));
                if(sum(optim_table(localinds2(inds(count3)),:))==0)
                    optim_table(localinds2(inds(count3)),:)=test;
                    continue
                end
                sol1=transpose(lsqnonneg((m1),squeeze(expression_matrix_processed(goodinds,count1,count2)),options));
                if(sol1((inds(count3))))
                    optim_table(localinds2(inds(count3)),:)=test;
                end
            end
        end
    end
end

%now improve astrocytes to match oligodendrocytes
eqz1=sum((clusters(goodinds,:)==0));

m1=[best_clusters(goodinds,localinds2)];

for countZ=1:4
    switch countZ
        case 0
            help_inds1=pvalb_inds1;
            start_opt=1; %1 3747
        case 4
            help_inds1=astro_inds1;
            start_opt=1850; % 1950
        case 2
            help_inds1=oligo_inds1;
            start_opt=1200;%1000
        case 3
            help_inds1=micro_inds1;
            start_opt=3080;%2980
        case 1
            help_inds1=endo_inds1;
            start_opt=1;%485
        case 5
            help_inds1=smc_inds1;
            start_opt=16820;
        case 6
            help_inds1=peri_inds1;
            start_opt=16820;
        case 7
            help_inds1=vlmc_inds1;
            start_opt=16820;

    end
    assigned=0;
    for countX=1:99
        countX
        for count11=1:length(help_inds1)
            % sol1=transpose(lsqnonneg((m1),transpose(optim_table(astro_inds1(count11),:)),options));
            count1=help_inds1(count11);

            rank_clusters(:,(count1))=best_rank_clusters(:,(count1));
            test_clusters(goodinds,(count1))=best_clusters(goodinds,(count1));

            [s1,i1]=sort( best_rank_clusters(:,(count1)));
            %inds0=find((i1>=start_opt).*(i1<=eqz1(count1)));
            inds0=i1(start_opt:eqz1(count1));

            rank_clusters=best_rank_clusters;
            test_clusters=best_clusters;


            %count1
            if (sum(optim_table(count1,:))==0)
                %     disp 'no solution possible'
                continue
            end
            m1=[best_clusters(goodinds,localinds2)];
            sol1=transpose(lsqnonneg((m1),transpose(optim_table(count1,:)),options));

            predvec1=(best_clusters(goodinds,localinds2) *transpose(sol1));
            origerr=sum((predvec1-transpose(optim_table(count1,:))).^2);
            % inds=find(clusters(goodinds,localinds2(count1))==0);
            sol2=zeros(1,b);
            sol2(localinds2)=sol1;
            remainderplus=sol2(count1)*best_clusters(goodinds,(count1))+ transpose(optim_table(count1,:))-predvec1;

             [s1, i1]=sort(remainderplus(inds0));

             rank_clusters(inds0(i1),count1)=start_opt-1+(1:length(i1));

            test_clusters(goodinds,(count1))=bestrefdist(rank_clusters(:,(count1)));

            m1=[test_clusters(goodinds,localinds2)];

            sol3=transpose(lsqnonneg((m1),transpose(optim_table(count1,:)),options));

            predvec3=(test_clusters(goodinds,localinds2) *transpose(sol3));
            testerr=sum((predvec3-transpose(optim_table(count1,:))).^2);


            if(testerr<(origerr))
                [ (1-testerr/origerr)]
                assigned=assigned+1;
                disp 'assigning'
              
                best_rank_clusters(:,(count1))=rank_clusters(:,(count1));
                best_clusters(goodinds,(count1))=test_clusters(goodinds,(count1));
                assigned=1;
                continue
            else
                disp 'miss'
            end


        end
        if(assigned==0) %none assigned for all inds
            break
        end
    end
end

%%%%%%%%%%
sol1=transpose(lsqnonneg(repmat(scaleVec(goodinds),1,length(localinds2)).*best_clusters(goodinds,localinds2),scaleVec(goodinds).*sumexpr0(goodinds),options));

figure
hold on
maxx=0;
maxy=0;
p1=[];
s1=[];
for count1=1:length(goodinds)
    sum1=sumexpr0(goodinds(count1));

    inds2=find(clusters0(goodinds(count1),localinds2));
    predval=sum(best_clusters(goodinds(count1),localinds2) .*sol1(1:(length(sol1))));
    plot(predval,sum1,'.')

    p1=[p1,predval];
    s1=[s1,sum1];
    maxx=max(predval,maxx);
    maxy=max(sumexpr0(goodinds(count1)),maxy);
end
xlim([0 maxx])
ylim([0 maxy])
plot([0 maxx],[0 maxx],'k')


clusters0=best_clusters;

predvec1=(best_clusters(goodinds,localinds2) *transpose(sol1(1:(length(sol1)))));
vec1=predvec1./(sumexpr0(goodinds));

figure
hold on
for count1=1:length(goodinds)
    sum1=sum(sum(sum(expression_matrix_processed(goodinds(count1),:,:))));
    predval=sum(clusters0(goodinds(count1),localinds2) .*sol1(1:(length(sol1))));
    plot(predval,sum1,'.')
end


%now make expression matrix match nissl sum
nissl_sum=sum(gtot2(:,3:49),2);
expr_mat_sum=sum(squeeze(sum(expression_matrix_processed(goodinds,:,:),1)),2);
rat1=nissl_sum./expr_mat_sum;
rat1=rat1./mean(rat1);
for count1=1:100
    expression_matrix_processed(goodinds,count1,:)=expression_matrix_processed(goodinds,count1,:).*rat1(count1);
end
figure
hold on
plot(nissl_sum./mean(nissl_sum))
plot(expr_mat_sum./mean(expr_mat_sum))

clusters=clusters0;
save('clusters','clusters')
save('goodinds','goodinds')
save('expression_matrix_processed','expression_matrix_processed');
save('max_base','max_base')
save('base_mat','base_mat')

%solve_densities
