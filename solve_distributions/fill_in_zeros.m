function [rank_clusters2] = fill_in_zeros(rank_clusters,astro_inds1,astro_rank)
rank_clusters2=rank_clusters;
for count1=1:length(astro_inds1)
    vec1=rank_clusters2(:,astro_inds1(count1));
    inds=find(vec1==0);
    [s1,i1]=sort(astro_rank(inds));
    vec1(inds(i1))=1:length(i1);
zerinds=find(astro_rank(inds)==0);
        vec1(inds(zerinds))=0;
    rank_clusters2(:,astro_inds1(count1))=vec1;    
end
end