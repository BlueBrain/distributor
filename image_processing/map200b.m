%perform the mapping

function biglist=map200b(newim,newsurflist,innervec,dir1,baseline)
[r s]=size(newsurflist);

biglist=zeros(100,1);
biglistcount=zeros(100,1);

[nrow ncol]=size(newim);
for count1=1:r
    startpoint=newsurflist(count1,:);
    finalpoint=innervec(count1,:);
    veclength=norm(finalpoint-startpoint );
    for count2=1:veclength
        currpos=(count2/veclength)*finalpoint+(1-count2/veclength)*startpoint;
        val=double(newim( max(1,min( round(currpos(2)),nrow)) , max(1,min(ncol,round(currpos(1))))   ));
        
        index=min(100,round(floor(1+ (count2/veclength)*100)));
        biglist(index)=biglist(index)+val;
        biglistcount(index)=biglistcount(index)+1;
    end
end
biglist=biglist./biglistcount;
