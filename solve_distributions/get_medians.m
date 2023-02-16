function [cell_rank] = get_medians(goodinds,rank_clusters,cell_inds1)
cell_rank=zeros(length(goodinds),1);
for count1=1:length(goodinds)
    inds=find(rank_clusters(count1,cell_inds1)>0);
    if(length(inds))
        cell_rank(count1)=median(rank_clusters(count1,inds));
   
    end
end
end