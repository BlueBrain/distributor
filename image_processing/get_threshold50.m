%returns threshold and also corrects image for baseline

function [center1,outside_maxd,snap1,insideMask]=get_threshold28(snap0,p1,p2,p3,p4);
x1=p1(1);
y1=p1(2);
x2=p2(1);
y2=p2(2);
x3=p3(1);
y3=p3(2);
x4=p4(1);
y4=p4(2);


snapoutside=snap0;
snap1=snap0;

slope1=(y4-y1)/(x4-x1);
[a0 b0]=size(snap0);
for countx=-10:b0+10;
    xstart= round(x1+countx);
    ystart=max(round(y1+slope1*(countx)),1);
    if (xstart>0)
        
        if (xstart<=b0)
            if (ystart<=a0)
                for countdown=ystart:-1:1
                    snap0(countdown,xstart)=-1;
                end
                
                for countup=ystart:1:a0
                    snapoutside(countup,xstart)=-1;
                end
                
                
                
            end
        end
    end
end

slope1=(y3-y2)/(x3-x2);
[a0 b0]=size(snap0);
for countx=-10:b0+10;
    xstart= round(x2+countx);
    ystart=round(y2+slope1*(countx));
    if (xstart>0)
        if (ystart>0)
            if (xstart<=b0)
                if (ystart<=a0)
                    for countdown=ystart:1:a0
                        snap0(countdown,xstart)=-1;
                    end
                end
            end
        end
    end
end

%find inflection point
[a0 b0]=size(snap0);
z=reshape(snap0,a0*b0,1);
z=z(find(z>0));
sorted1=sort(z);

z=reshape(snapoutside,a0*b0,1);
z=z(find(z>0));
sorted2=sort(z);

nlist=[];
dlist=[];
oldval=0;
outside_dlist=[];
outside_oldval=0;
for count1=1:255
    currval=sum(sorted1<=count1);
    nlist=[nlist;currval];
    dlist=[dlist;currval-oldval];
    oldval=currval;
    
    currval=sum(sorted2<=count1);
    outside_dlist=[outside_dlist;currval-outside_oldval];
    outside_oldval=currval;
end
%dlist(255)=0; %saturation can occur
%outside_dlist(255)=0; %saturation can occur
maxd=find(dlist==max(dlist(1:254)));
maxd=maxd(1);


maxd0=maxd;

outside_maxd=find(outside_dlist==max(outside_dlist(maxd:254)));
outside_maxd=outside_maxd(find(outside_maxd>=maxd ));
outside_maxd=outside_maxd(1);

dlist2=dlist;

%identify secondpeak

%are there any places where it increases?
for count2=maxd+1:length(dlist)
    if (dlist(count2)>dlist(count2-1))
        secondPeakVal=max(dlist(count2+1:length(dlist)));
        
        secondPeakInd=find(dlist(count2+1:length(dlist))==secondPeakVal);
if (length(secondPeakInd)==0)
    continue;
end
        secondPeakInd=secondPeakInd(1);
        
        if (secondPeakInd+count2-maxd>10) %if is sufficiently offset from the maxd
            maxdSecond=count2+secondPeakInd;
            
if (maxdSecond>=length(dlist2))
  continue
end
    %find the real symmetry point
            %by the average of the two adjacent neighbors
            vallist1=sort(dlist2(maxdSecond-1:maxdSecond+1));
            vallist1=vallist1(2); %choose the middle value
            
            if(dlist2(maxdSecond-1)==vallist1)
disp 'here1'
                xval=(-vallist1+dlist2(maxdSecond))/(dlist2(maxdSecond)-dlist2(maxdSecond+1)+1e-9 )+maxdSecond;
                center1=(xval+(maxdSecond-1))/2;
                else
disp 'here2'
                xval=(vallist1-dlist2(maxdSecond-1))/(dlist2(maxdSecond)-dlist2(maxdSecond-1)+1e-9 )+maxdSecond-1;
                center1=(xval+(maxdSecond+1))/2;
            end
            

            %center1 is now the real symmetry point
            % center1
            
            rightside2=dlist(secondPeakInd+count2:length(dlist));
            rightsidex2=(secondPeakInd+count2:length(dlist))-center1;
            

            dlist2(  round(floor(center1)): length(dlist2) )=0;
            
            for count3=round(floor(center1)):-1: round(floor(center1-max(rightsidex2)))+1;
                testDist=center1-count3;
                dlist2(count3)=dlist(count3)-interp1(rightsidex2,rightside2,testDist);
            end
        end
        
        %now let dlist2 only be monotonically decreasing from secondPeakInd
        for count3=count2-1:length(dlist2)-1
            slopes=(dlist2(count3+1:length(dlist2))-dlist2(count3));
            slopes=slopes./transpose([1:length(slopes)]);
            min1=min(slopes);
            dlist2(count3+1)=dlist2(count3)+min1;
        end

     %   figure
     %   plot(dlist2,'g')
        break
    end
end


%repeatFlag=1;

%now identify new center
vallist1=sort(dlist2(maxd-1:maxd+1));
vallist1=vallist1(2); %choose the middle value

if(dlist2(maxd-1)==vallist1)
    xval=(-vallist1+dlist2(maxd))/(dlist2(maxd)-dlist2(maxd+1) )+maxd;
    center0=(xval+(maxd-1))/2;
else
    xval=(vallist1-dlist2(maxd-1))/(dlist2(maxd)-dlist2(maxd-1) )+maxd-1;
    center0=(xval+(maxd+1))/2;
end

%center0 is now the real symmetry point


olddlist2=dlist2;
if(1)
for testCount=1:40
    testdlist2=subtract_center(maxd,center0,dlist2);
    if (sum(testdlist2(maxd-10:maxd-2)<0)>0) %if there is a dip
        center0=center0+0.1;
    else
        break
    end
end
end

dlist2=testdlist2;
center0

center1=round(floor(center0));



dlist2=dlist2.*(dlist2>=0);
dlist3=dlist2;
maxd2=center0;

%identify the minimum place to start kernels

dlist0=dlist;
%now do the mapping from dlist to dlist3
%and adjust the snaps
newdlist=[];
for count1=1:length(dlist)
    for count2=1:dlist(count1)
        newdlist=[newdlist;count1];
    end
end

%figure
%plot(dlist3)
newdlist3=[];
for count1=1:length(dlist3)
    for count2=1:round(dlist3(count1))
        newdlist3=[newdlist3;count1];
    end
end

for count1=length(newdlist3)+1:length(newdlist)
    newdlist3=[newdlist3;center1+1];
end

[a b]=size(snap0);

sum1=1;
index1=[];
for count1=1:length(dlist)
  index1=[index1;sum1];
  sum1=sum1+dlist(count1);
end

%replace values in images
for count1=1:a
  %  count1
    for count2=1:b
        val=snap0(count1,count2);
        if (val<1)
            continue
        end
        ind1=round(floor(val));
        for count3= index1(ind1): length(newdlist)          
            if ( newdlist(count3)==ind1 )
                snap0(count1,count2)=newdlist3(count3);
                snap1(count1,count2)=newdlist3(count3);
                newdlist(count3)=-222;     
                index1(ind1)=index1(ind1)+1;
                break
            end
        end
    end
end


disp 'done dlist mapping'

snap1=snap1.*(snap1>=0);


%%%%%%%%%%%%%%%%%%%%%%5
%remove later
