
%make the depth index used to bin the map
function [mapindex]=make_mapindex(newim,newsurflist,innervec,dir1,min1,min2,pointsL4)

[r s]=size(newsurflist);

[nrow ncol]=size(newim);
mapindex=sparse(nrow,ncol);

for count0=0:1 % a switch for half step
for count1=1:r-count0
if (count0==0)
    startpoint=newsurflist(count1,:);
    finalpoint=innervec(count1,:);
else
  startpoint=(newsurflist(count1,:)+newsurflist(count1+1,:))./2;
  finalpoint=(innervec(count1,:)+innervec(count1+1,:))./2;
end
    veclength=norm(finalpoint-startpoint );
    for count2=1:0.5:veclength
        currpos=(count2/veclength)*finalpoint+(1-count2/veclength)*startpoint;
        testx=max(1,min( round(currpos(2)),nrow));
        testy=max(1,min(ncol,round(currpos(1))));
        if (newim(testx,testy)>=0)
            if (length(pointsL4)>0)
                xa1=pointsL4(1,1);
                ya1=pointsL4(1,2);
                xa2=pointsL4(2,1);
                ya2=pointsL4(2,2);
                
                xb1=startpoint(1);
                yb1=startpoint(2);
                xb2=finalpoint(1);
                yb2=finalpoint(2);
                
                xsol=((xa1* (-ya1 + ya2))/(-xa1 + xa2) - ( xb1* (-yb1 + yb2))/(-xb1 + xb2))/((-ya1 + ya2)/(-xa1 + xa2) - (-yb1 +   yb2)/(-xb1 + xb2));
                ysol=(ya2-ya1)/(xa2-xa1)*(xsol-xa1)+ya1;
                
                ysol=(ya2-ya1)/(xa2-xa1)*(xsol-xa1)+ya1;
                veclength45= sqrt((xsol-xb1)*(xsol-xb1)+(ysol-yb1)*(ysol-yb1));
                if (count2<veclength45)
                    index=round(floor(1+ (count2/veclength45)*45));
                else
                    index=round(floor(46+ (count2-veclength45)/(veclength-veclength45)  *54));
                    index=min(100,index);
                end
                            
            else
                index=min(100,round(floor(1+ (count2/veclength)*100)));
            end
	    mapindex(testx,testy)=index;
            newim(testx,testy)=-1;            
        end
    end
end
end
