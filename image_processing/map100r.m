
%perform the mapping
function [biglist1,map4,map7,heightMask]=map100r(newim0,newim,mapindex,snap_area,min1,min2,areaScale,increment)
%account for -1 baseline
%[r s]=size(newsurflist);

heightMask=0.* snap_area-1;
biglist1=zeros(100,1);
%biglist2=zeros(100,1);
biglistcount=zeros(100,40);
allbiglistcount=zeros(100,1);
map4=zeros(100,40);
map7=zeros(100,40);

[nrow ncol]=size(newim);
for testx=1:nrow
    for testy=1:ncol
        if (mapindex(testx,testy)>0)
            index=mapindex(testx,testy);
            
            val=double(newim0(testx,testy));
            biglist1(index)=biglist1(index)+val;
            if (newim(testx,testy)>0)
                %test the areas
                area=double(snap_area(testx,testy));
                
                index1=area/increment;
                
                if (index1>=40)
                    index1=39;
                end
                
                index1=round(floor(index1));
                heightMask(testx,testy)=index1;
                map4(index,index1+1)=map4(index,index1+1)+val;%./rad1;%volume;
                map7(index,index1+1)=map7(index,index1+1)+1;
            end
            allbiglistcount(index)=allbiglistcount(index)+1;
        end
    end
end

biglist1=biglist1./(areaScale .* allbiglistcount);

for count1=1:100
    map4(count1,:)=map4(count1,:)./(areaScale .* allbiglistcount(count1)+1e-9);
    map7(count1,:)=map7(count1,:)./(areaScale .* allbiglistcount(count1)+1e-9);
end


