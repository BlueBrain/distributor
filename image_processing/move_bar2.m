%identify the surface
function surflist=move_bar2(newim,surflist,dir1,dir2,maxlength)

allvec=[];
[s r]=size(surflist);
lenvec=norm(surflist(1,:)-surflist(s,:));

newsurflist=[];
for count0=1:s
    start1=surflist(count0,:);
    
    %find the lowest point in the range of-50 to +100 of the input length
    margin1=10;%was 40
    extract1=[];
    basestart=start1-(dir1).*(margin1);
    extract1=[];
    
    for count2=0:2*margin1
        currstart=basestart+dir1.*(count2);
        myvec=newim(round(max(1,currstart(2))),round(currstart(1)) );
        extract1=[extract1;(myvec)]; %was mean
    end
    
    
    margin2=5;
    
    extract2=[];
    for count1=(margin1-margin2):2*margin1-margin2
        prior=mean(extract1(count1-margin2+1:count1-1));
        next=mean(extract1(count1:count1+margin2-1));
        extract2=[extract2;next-prior];
    end
    %plot(extract2)
    maxind=find(extract2==max(extract2));
    bestind=maxind(1)-margin2;
    newpos=start1+bestind.*dir1;
    newsurflist=[newsurflist;newpos];
end

surflist=newsurflist;


[length_surfacelist r]=size(surflist);
%remove stragglers
for count2=2:length_surfacelist-1
    curr=surflist(count2,:);
    prev=surflist(count2-1,:);
    next=surflist(count2+1,:);
    curr=max(curr,[1,1]);
    prev=max(prev,[1,1]);
    next=max(next,[1,1]);
    if (( sum((prev-curr).*dir1)>10  )*( sum((next-curr).*dir1)>10  ))
        surflist(count2,:)=(prev +next)./2;
    end
end

[dir1,dir2]=fit_vector(surflist);
surflist=adjust_surface(surflist,dir1,100);

%remove forerunners
xmedian=(median(surflist(:,1)));
ymedian=(median(surflist(:,2)));
for cutoff=25.0:-0.25:1.0
    for count2=2:length_surfacelist
        curr=surflist(count2,:);
        prev=surflist(count2-1,:);
        
        testdist1=sum((curr-prev).*dir1);
        if (( testdist1>cutoff  ) )
            prevdist=sum((prev-[xmedian ymedian]).*dir1);
            currdist=sum((curr-[xmedian ymedian]).*dir1);
            if (abs(currdist)>abs(prevdist))
                surflist(count2,:)=surflist(count2,:)-dir1 .*testdist1/2;
            end
            
        end
    end
    
    for count2=length_surfacelist-1:-1:1
        curr=surflist(count2,:);
        next=surflist(count2+1,:);
        
        testdist1=sum((curr-next).*dir1);
        if (( testdist1>cutoff  ) )
            nextdist=sum((next-[xmedian ymedian]).*dir1);
            currdist=sum((curr-[xmedian ymedian]).*dir1);
            
            if (abs(currdist)>abs(nextdist))
                surflist(count2,:)=surflist(count2,:)-dir1 .*testdist1/2;
            end
        end
        
    end
end



[newdir,newdir2]=fit_vector(surflist);

%now a final cull of surface
%surflist=adjust_surface(surflist,newdir,30);

[length_surfacelist col]=size(surflist);

cutoff=2;
xmedian=(median(surflist(:,1)));
ymedian=(median(surflist(:,2)));
for count0=1:20
    for count1=2:length_surfacelist
        prev=surflist(count1-1,:);
        curr=surflist(count1,:);
        testvec=surflist(count1,:)-prev;
        testdist1=norm(testvec);
        if(testdist1>cutoff)
            surflist(count1,:)=(2.*curr+prev)./3;
            surflist(count1-1,:)=(curr+2.*prev)./3;
        end
    end
    
    for count1=length_surfacelist-1:-1:1
        next=surflist(count1+1,:);
        curr=surflist(count1,:);
        testvec=surflist(count1,:)-next;
        testdist1=norm(testvec);
        if(testdist1>cutoff)
            surflist(count1,:)=(2.*curr+next)./3;
            surflist(count1+1,:)=(curr+2.*next)./3;
        end
    end
    
    
end

%final remove stragglers in the right direction
[newdir1,newdir2]=fit_vector(surflist);
for count3=1:50
    for count2=2:length_surfacelist-1
        curr=surflist(count2,:);
        prev=surflist(count2-1,:);
        next=surflist(count2+1,:);
        curr=max(curr,[1,1]);
        prev=max(prev,[1,1]);
        next=max(next,[1,1]);
        if ( sum((prev-curr).*newdir1)>2  )
            surflist(count2,:)=(2 *prev +curr)./3;
        else
            if ( sum((next-curr).*newdir1)>2  )
                surflist(count2,:)=(2 *next +curr)./3;
            end
        end
    end
    for count2=length_surfacelist-1:-1:2
        curr=surflist(count2,:);
        prev=surflist(count2-1,:);
        next=surflist(count2+1,:);
        curr=max(curr,[1,1]);
        prev=max(prev,[1,1]);
        next=max(next,[1,1]);
        if ( sum((prev-curr).*newdir1)>2  )
            surflist(count2,:)=(2 *prev +curr)./3;
        else
            if ( sum((next-curr).*newdir1)>2  )
                surflist(count2,:)=(2 *next +curr)./3;
            end
        end
    end
end


if (0)
    [s, r]=size(surflist);
    for count2=1:s
        start1=surflist(count2,:);
        newim(start1(2):start1(2)+0,start1(1):start1(1)+0)=0;
        
    end
    imshow(newim)
    %pause
end



