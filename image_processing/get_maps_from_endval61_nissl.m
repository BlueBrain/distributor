%make the maps using image processing
function [map2,map4,map7,heightMask,cumulativeImageExpression,allDotsArray]=get_maps_from_endval61(endval2,snap0,baseline,surflist,innervec,dir1,snapMin1,snapMin2,snapMax1,snapMax2,areaScale,pointsL4,r1,r2)

increment=20; %for areas

disp 'making masks'
bestnissl=dlmread('NisslAvg61.dat');
bestnisslsumsum=sum(sum(transpose(bestnissl)));
sumbestnissl0=sum(bestnissl);

snap1=snap0.*(snap0<endval2)+ baseline.*(snap0>=endval2);
newimage=baseline-snap1;
[s1 s2]=size(newimage);

disp 'making index'

%restore this
[mapindex]=make_mapindex(sparse(r1,r2),surflist,innervec,dir1,snapMin1,snapMin2,pointsL4);
mapindex=full(mapindex(snapMin2:snapMax2,snapMin1:snapMax1));

disp 'making masks'

cells1=reflector6(snap0,endval2,mapindex);

%find a new snap_Area
s = regionprops(cells1>0,'Centroid','Area');

snap_area=0.*cells1;
for count1=1:length(s)
    pt=max(round(s(count1).Centroid),1);
    snap_area(pt(2),pt(1))=s(count1).Area;
end


newimage=baseline-snap0;
cc = bwconncomp(cells1>0);
labeled = labelmatrix(cc);
densityMask=0.*cells1;
max1=max(max(labeled));
for count1=1:max1
    im1=(labeled==count1);
    densityMask=densityMask+ im1.*  (sum(sum(im1.*newimage)) ./sum(sum( im1))) ;
end

%now expand the densityMask
densityMask2=densityMask;
[a b]=size(densityMask);
for count0=1:7
    for count1=2:a-1
        for count2=2:b-1
            if (densityMask2(count1,count2 )==0 )
                v1=reshape(densityMask(count1-1:count1+1,count2-1:count2+1 ),9,1);
                if (sum(v1)>1)
                    v2=v1(v1>0);
                    m1=mode(v2);
                    densityMask2(count1,count2 )=m1(1);
                end
            end
        end
    end
    densityMask=densityMask2;
end
densityMask=densityMask.*(snap_area>0);


disp 'mapping'

[map2,map4,map7,heightMask]=map100r(newimage,densityMask,mapindex,snap_area,snapMin1,snapMin2,areaScale,increment);
sumsummap7=sum(sum(transpose(map7)));



disp 'done'

maxsnap=max(max(snap0));
%write the dots mask
%fraction density  maxval scaled_total_dots_at_this_fraction scaled_area
scaled_dot_fractions=[];
for count1=1:100
    scaled_dot_fractions=[scaled_dot_fractions;areaScale .* sum(sum( mapindex==count1 ))];
end
allDotsArray=[];
[s1 s2]=size(densityMask);

numdots=0;
for count1=1:s1
    for count2=1:s2
        if (densityMask(count1,count2)>0 )           
               % snap0(count1,count2)=2.*maxsnap;
                numdots=numdots+1;
                position=mapindex(count1,count2);
                if (position>0)
                    allDotsArray=[allDotsArray;position,densityMask(count1,count2), 0,areaScale * (snap_area(count1,count2)), scaled_dot_fractions(position)];
                end
            
        end
        
    end
end



cumulativeDots=sort((allDotsArray(:,2)-1).*allDotsArray(:,4) );

[a b]=size(densityMask);
d3=reshape(densityMask,a*b,1);
d3=d3(d3>0);
d3=sort(d3);
if (length(d3)>0)
cumulativeImageExpression=d3(max(1,round([1:100]./100 *length(d3))));
else
    cumulativeImageExpression=0.*[1:100];
end



