%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%if a point is beyond a certain limit from the projection onto the normal dir, adjust it

function surflist=adjust_surface(surflist,dir1,cutoff)
xmedian=round(median(surflist(:,1)));
ymedian=round(median(surflist(:,2)));

%cutoff=100;
for count2=1:length(surflist)
    curr=surflist(count2,:);
    vec1=curr-[xmedian,ymedian];
    projlength=abs(sum(- vec1 .* dir1));
    if (projlength>cutoff)
        surflist(count2,:)=surflist(count2,:)+dir1 .* (projlength-cutoff);
    end
end
