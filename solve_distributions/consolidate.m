function [clusters2] = consolidate(clusters2,clusters00,goodinds,inds0)
vec1=(sum(clusters2(goodinds,:))>0);
inds1=inds0(find(vec1(inds0)));
if(length(inds1)==0)
    return;
end
vec1=[];
for count2=1:length(goodinds)
    inds=find(clusters00(goodinds(count2),inds1)>0);
 %  inds=1:length(inds1);
    med1=0;
    if(length(inds))
        med1=median(clusters2(goodinds(count2),inds1(inds)));
    end
    vec1=[vec1;med1];
end
clusters2(goodinds,inds1)=0;
clusters2(goodinds,inds0(1))=vec1;
end