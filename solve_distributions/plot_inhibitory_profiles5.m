%plot the inhibitory profiles
load ttype_labels
load layer_array

totnum=0;
for count0=1:length(inh_inds1)
    count1=inh_inds1(count0);
    if(sum(allvec_individual_vector5(count1,:)))
        totnum=totnum+1;
    end
end

numcols=24;
numrows=floor(totnum/numcols);
if(rem(totnum,numcols)>0)
    numrows=numrows+1;
end

numcols=5;
numrows=6;

sumInh=sum(transpose(allvec_individual_vector5(inh_inds1,:)));
[foo,i1]=sort(sumInh,'Descend');

figure
totnum=1;
max1=0;
for count0=1:30%length(inh_inds1)
    count1=inh_inds1(i1(count0));
    if(sum(allvec_individual_vector5(count1,:)))
        subplot(numcols,numrows,totnum)
        totnum=totnum+1;
        box off
        plot(allvec_individual_vector5(count1,:))
        hold on
        max1=max(max1,max(allvec_individual_vector5(count1,:)));
        sum1=sum(allvec_individual_vector5(count1,:));
        vec1=zeros(1,100);
        vec1(1:l12)=layer_array(count1,1);%./length(vec1(1:l12));
        vec1(l12+1:l23)=layer_array(count1,2);%./length(vec1(l12+1:l23));
        vec1(l23+1:l34)=layer_array(count1,3);%./length(vec1(l23+1:l34));
        vec1(l34+1:l45)=layer_array(count1,4);%./length(vec1(l34+1:l45));
        vec1(l45+1:l56)=layer_array(count1,5);%./length(vec1(l45+1:l56));
        vec1(l56+1:100)=layer_array(count1,6);%./length(vec1(l56+1:100));

        normsum=sum(vec1);

        
        ylim([0.0 1.1].*max1);
        plot([l12 l12],[0 max1],'--','color',[.7 .7 .7])
        plot([l23 l23],[0 max1],'--','color',[.7 .7 .7])
        plot([l34 l34],[0 max1],'--','color',[.7 .7 .7])
        plot([l45 l45],[0 max1],'--','color',[.7 .7 .7])
        plot([l56 l56],[0 max1],'--','color',[.7 .7 .7])
        
        title(replace(ttype_labels(count1),'_','-'),'FontSize', 4);
    end
end


figure
totnum=1;
max1=0;
for count0=1:30%length(inh_inds1)
    count1=inh_inds1(i1(count0));

    vec1=squeeze(sum( allvec_individual6(count1,:,:),2));
    if(sum(allvec_individual_vector5(count1,:)))
if(sum(vec1)==0)
    count1
    pause
end
        subplot(numcols,numrows,totnum)
        totnum=totnum+1;
        box off
        plot(vec1./sum(vec1))

        title(replace(ttype_labels(count1),'_','-'),'FontSize', 4);
    end
end