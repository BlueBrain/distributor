function [ finalIm ] = reflector5( snap0,baseline,mapindex)
newsnap=[];
snapTemp=im2bw(snap0<=baseline);
%find additional cells by eroding the dark spots
%imshow(imerode(snap1==0,1))

snapStart = bwpropfilt(snapTemp,'Area',[10 Inf]);
snap1=snapStart;

cells1=0.*snap1;
for inc1=0:5:35
    %do convolution ethod for finding dark centers
    kern3=ones(9,9);
    kern3(4:6,4:6)=0;
    test1=conv2(1.0 .* snap1,1.0 .* kern3);
    [a b]=size(snap1);
    test2=test1(5:(4+a),5:(4+b));
    darkconv3=(test2>20+inc1); %60 35
    
    kern2=ones(6,6);
    kern2(3:4,3:4)=0;
    test1=conv2(1.0 .* snap1,1.0 .* kern2);
    [a b]=size(snap1);
    test2=test1(3:(2+a),3:(2+b));
    darkconv2=(test2>15+inc1); %30 25
    
    if (inc1<30)
        darkCenters = bwpropfilt(snap1==0,'Area',[1 30]);
        erode1=imerode(snap1==0,1);
        darkCenters0 = bwpropfilt(erode1,'Area',[1 30]);
        dC= darkconv2+darkconv3+darkCenters+darkCenters0;
    else
        dC= darkconv2+darkconv3;
    end
    if (inc1>30)
        dC=0.*dC;
    end
    %imshow(darkconv2+darkconv3+(baseline-snap0)./200)
    
    snap1=snap1+dC;
        
    
    %now pull out cell-sized objects
    cells0 = bwpropfilt(snap1>0,'Area',[1 300]);
    cc = bwconncomp(cells0);
    labeled = labelmatrix(cc);
    max1=max(max(labeled));
    
    for count1=1:max1
        im1=labeled==count1;
        xproj=sum(im1)>0;
        yproj=sum(transpose(im1))>0;
        xdim=sum(xproj);
        ydim=sum(yproj);
        if (xdim<=15)
            if (ydim<=15)
                cells1=cells1+im1;
            end
        end
    end
    
        
    snap1=snapStart.*(cells1==0);
end


%now start increasing the threshold

for minthresh=15:3:27
    minthresh
    for inc1=0:10:160
        snap1=(snap0<=(baseline-inc1));
        
        %now pull out cell-sized objects
        cells0 = bwpropfilt(snap1>0,'Area',[30  (0.9 * minthresh.^2) ]);
        cc = bwconncomp(cells0);
        labeled = labelmatrix(cc);
        max1=max(max(labeled));
        
        for count1=1:max1
            im1=labeled==count1;
            xproj=sum(im1)>0;
            yproj=sum(transpose(im1))>0;
            xdim=sum(xproj);
            ydim=sum(yproj);
            if (xdim<=minthresh)
                if (ydim<=minthresh)
                    cells1=cells1+im1;
                end
            end
        end
    end
end
figure
imshow(cells1+(baseline-snap0)./200);
title('cells1')

finalIm=cells1;

if(0)
    snap1=(snap0<=(baseline-120)) .* (cells1==0);
    
    
    
    if(0)
        N=1;
        se = ones(2*N + 1, 2*N + 1);
        bw2 = imdilate(dC,se);
        bw2=bw2-dC;
    end
    %imshow(bw2+(baseline-snap0)./100);
    %pause
    
    %restore
    %snap1=snap1.*(bw2==0);
    %snap1=dC;
    
    %imshow(snap1+(baseline-snap0)./100);
    %pause
    %imshow(snap1+(baseline-snap0)./100);
    %pause
    
    %decide what total density we need
    % by calculating volume and desired cell num
    % then erode objects
    % when they reach a certain size save and delete from erode image
    % also save the erode number
    % repeat till desired number achieved
    
    % repopulate center point with circle proportional to erode number and radius
    % then expand till space filled again
    
    % tally the number allocated to each spot
    % or, alternatively count the lines to create a signature
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %what is the total number of cells we need?
    maxcelldensity=200000; %cells per mm tissue
    
    %what is the volume of the area of interest?
    
    allowedAreas=0.* snap1;
    %erode image by four to avoid assigning centers at periphery
    if(0)
        [a b]=size(snap0);
        allowedAreas=snap1;
        tempsnap1=snap1;
        for count0=1:2 %was 2
            for count1=2:1:a-1
                for count2=2:b-1
                    if (tempsnap1(count1,count2))
                        if (sum(sum( tempsnap1(count1-1:count1+1,count2-1:count2+1)>0   ))==1)
                            continue
                        end
                        if (count1-1>0)
                            if (allowedAreas(count1-1,count2)==0)
                                tempsnap1(count1,count2)=0;
                            end
                        end
                        if (count2-1>0)
                            if (allowedAreas(count1,count2-1)==0)
                                tempsnap1(count1,count2)=0;
                            end
                        end
                        if (count1+1<=a)
                            if (allowedAreas(count1+1,count2)==0)
                                tempsnap1(count1,count2)=0;
                            end
                        end
                        if (count2+1<=b)
                            if (allowedAreas(count1,count2+1)==0)
                                tempsnap1(count1,count2)=0;
                            end
                        end
                    end
                end
            end
            allowedAreas=tempsnap1;
        end
    end
    
    allowedAreas=allowedAreas+dC;
    %figure
    %imshow(allowedAreas .* 0.5 + snap1.*0.25 +(baseline-snap0)./100)
    
    cc = bwconncomp(allowedAreas,4);
    labeled = labelmatrix(cc);
    
    cc2 = bwconncomp(snap1,4);
    labeled2 = labelmatrix(cc2);
    
    allowed_to_base=zeros(1,max(max(labeled)) );
    [a b]=size(labeled);
    for count1=1:a
        for count2=1:b
            if ( labeled(count1,count2) >0)
                
                if ( allowed_to_base( labeled(count1,count2)  )==0  )
                    allowed_to_base( labeled(count1,count2)  ) = labeled2(count1,count2) ;
                    continue
                end
                if ( allowed_to_base( labeled(count1,count2)  )<0  )
                    continue
                end
                if (allowed_to_base( labeled(count1,count2)  )~= labeled2(count1,count2) )
                    allowed_to_base( labeled(count1,count2)  ) = -1 ;
                end
            end
        end
    end
    
    %now do the reverse mapping
    allowed_to_base2=zeros(1,max(max(labeled2)) );
    [a b]=size(labeled2);
    for count1=1:a
        for count2=1:b
            if ( labeled2(count1,count2) >0)
                if ( labeled(count1,count2) >0)
                    
                    if ( allowed_to_base2( labeled2(count1,count2)  )==0  )
                        allowed_to_base2( labeled2(count1,count2)  ) = labeled(count1,count2) ;
                        continue
                    end
                    if ( allowed_to_base2( labeled2(count1,count2)  )<0  )
                        continue
                    end
                    if (allowed_to_base2( labeled2(count1,count2)  )~= labeled(count1,count2) )
                        allowed_to_base2( labeled2(count1,count2)  ) = -1 ;
                    end
                end
            end
        end
    end
    
    
    finalim=0.*snap1;
    %now only keep ones that are allowed
    for count1=1:length(allowed_to_base)
        if (allowed_to_base(count1)>0)
            if (allowed_to_base2(allowed_to_base(count1))>0)
                finalim=finalim+(labeled2==allowed_to_base(count1));
            end
        end
    end
    %imshow(finalim+(baseline-snap0)./100)
    
    
    
    %save fooo
    %title('allowed areas after eroding')
    %pause
    
    
    if(1)
        %add in local minima
        for count1=2:1:a-1
            for count2=2:b-1
                
                max1=max(max( snap0(count1-1:count1+1,count2-1:count2+1) ));
                if (sum(sum( snap0(count1-1:count1+1,count2-1:count2+1) ==max1 ))==1)
                    if (snap0(count1,count2)==max1)
                        allowedAreas(count1,count2)=1;
                    end
                else
                    allowedAreas(count1,count2)=0;
                end
            end
        end
    end
    
    if(0)
        %add in local maxima
        for count1=2:1:a-1
            for count2=2:b-1
                
                min1=min(min( snap0(count1-1:count1+1,count2-1:count2+1) ));
                if (sum(sum( snap0(count1-1:count1+1,count2-1:count2+1) ==min1 ))==1)
                    if (snap0(count1,count2)==min1)
                        allowedAreas(count1,count2)=1;
                    end
                else
                    allowedAreas(count1,count2)=0;
                end
            end
        end
        allowedAreas=allowedAreas+dC;
        
        %allowedAreas=dC;
        
        figure
        imshow(allowedAreas .* 1 + (baseline-snap0)./100)
        title('allowed areas added local maxima')
        % pause
    end
    
    %allowedAreas=bwmorph(snap1,'skel',Inf);
    
    angIm=(0.* snap1);
    
    [a b]=size(snap0);
    
    for count1=3:1:a-3
        for count2=3:b-3
            
            if (allowedAreas(count1,count2))
                if (snap1(count1,count2))
                    %	      disp 'here'
                    %		pause
                    vecMat=[];
                    angradlist=[];
                    for ang=-pi/2:pi/6:pi/2
                        
                        vec1=[cos(ang) sin(ang)];
                        
                        dist1=0;
                        ampvec1=[];
                        for inc1=1:10
                            x=max(min(count1+round(inc1 .* vec1(1)),a),1);
                            y=max(min(count2+round(inc1 .* vec1(2)),b),1);
                            dist1=inc1;
                            ampvec1=[ ampvec1,baseline-snap0(x,y)];
                            if (snap1(x,y)==0)
                                
                                break
                            end
                        end
                        ampvec2=[];
                        dist2=0;
                        for inc1=1:10
                            x=max(min(count1-round(inc1 .* vec1(1)),a),1);
                            y=max(min(count2-round(inc1 .* vec1(2)),b),1);
                            ampvec2=[ ampvec2,baseline-snap0(x,y)];
                            dist2=inc1;
                            if (snap1(x,y)==0)
                                
                                break
                            end
                        end
                        %                    angradlist=[angradlist;min(dist1,dist2)];
                        distmin=max(dist1,dist2);
                        if (dist1>0)
                            if (dist2>0)
                                distmin=min(dist1,dist2);
                            end
                        end
                        if(distmin>0)
                            sub1=abs(ampvec2(1:distmin)-ampvec1(1:distmin));
                            max1=max([ampvec2(1:distmin);ampvec1(1:distmin)]);
                            angradlist=[angradlist;sum( 1- sub1./(max1+1e-9))];
                            
                            vec1=zeros(1,100);
                            vec1(1:length(ampvec1))=vec1(1:length(ampvec1))+ampvec1;
                            vec1(1:length(ampvec2))=vec1(1:length(ampvec2))+ampvec2;
                            vec1(1:distmin)=vec1(1:distmin)./2;
                            vec1=[baseline-snap0(count1,count2),vec1];
                            vec1=max(vec1,0);
                            if (max(vec1)>vec1(1))
                                vecMat=[vecMat;1];
                            else
                                vecMat=[vecMat;0];
                            end
                            %       vecMat=[vecMat;vec1];
                        else
                            angradlist=[angradlist;0];
                        end
                    end
                    %           figure
                    %		    vecMat
                    %                plot(transpose(vecMat))
                    %                pause
                    %  angIm(count1,count2)=mean(angradlist);
                    %angIm(count1,count2)=median((angradlist));
                    angradlist=sort(angradlist);
                    
                    %only add in if the it is mostly montonically decreasing
                    %               if (sum(vecMat)<0.25 .*length(vecMat))
                    %  disp 'here'
                    %  pause
                    angIm(count1,count2)=angradlist(round(0.5 * length( angradlist )));
                    %                end
                    %  angIm(count1,count2)=1;
                end
            end
        end
    end
    
    
    figure
    % imshow((angIm~=0) .* 1 + (baseline-snap0)./200)
    imshow((angIm~=0) .* 1+ (baseline-snap0)./100 )
    %	      sum(sum(angIm))
    %    pause
    title('Angles')
    
    %now assign some dots
    sortList=zeros(a*b,3);
    count0=1;
    for count1=3:1:a-3
        for count2=3:b-3
            if (angIm(count1,count2)>1)
                %sortList=[sortList;angIm(count1,count2),count1,count2];
                sortList(count0,:)=[angIm(count1,count2),count1,count2];
                count0=count0+1;
            end
        end
    end
    sortList=sortList(1:count0-1,:);
    sortList=sortrows(sortList,[-1 2]);
    
    doneList=zeros(a,b);
    [len1 foo]=size(sortList);
    
    finalIm=0.*snap0;
    max1=max(max(snap0));
    doneCount=0;
    for count0=1:len1
        count1=sortList(count0,2);
        count2=sortList(count0,3);
        if(doneList(count1,count2)==0)
            %check to make sure there are no pixels within 4
            cutoff1=5;
            start1=max(count1-cutoff1,1);
            start2=max(count2-cutoff1,1);
            end1=min(count1+cutoff1,a );
            end2=min(count2+cutoff1,b );
            
            if( sum(sum( ( finalIm(start1:end1,start2:end2)>0 )   )) >0 )
                %then double check in more detail
                continueflag=0;
                for count3=start1:end1
                    for count4=start2:end2
                        if (finalIm(count3,count4)>0)
                            
                            dist1=sqrt(sum([count3-count1,count4-count2].^2));
                            if (dist1>0)
                                if (dist1<=cutoff1)
                                    % disp 'found'
                                    % pause
                                    continueflag=1;
                                    break
                                end
                            end
                        end
                    end
                    
                    if (continueflag)
                        break
                    end
                end
                
                if (continueflag)
                    continue
                end
            end
            
            finalIm(count1,count2)=2.*max1;
            %angInc=1/(20*angIm(count1,count2));
            angInc=atan(1/angIm(count1,count2));
            
            for ang=-pi/2:angInc:pi/2+angInc
                
                vec1=[cos(ang) sin(ang)];
                
                dist1=0;
                for inc1=1:13
                    x=max(min(count1+round(inc1 .* vec1(1)),a),1);
                    y=max(min(count2+round(inc1 .* vec1(2)),b),1);
                    
                    
                    if (snap1(x,y)==0)
                        break
                    end
                    dist1=inc1;
                end
                
                %if(dist1==0)
                %    disp 'dist1 0'
                %end
                
                dist2=0;
                for inc1=1:13
                    x=max(min(count1-round(inc1 .* vec1(1)),a),1);
                    y=max(min(count2-round(inc1 .* vec1(2)),b),1);
                    
                    
                    if (snap1(x,y)==0)
                        
                        break
                    end
                    dist2=inc1;
                end
                
                rad1=1.7 .* min(dist1,dist2);
                for inc1=-rad1:rad1
                    x=max(min(count1+round(inc1 .* vec1(1)),a),1);
                    y=max(min(count2+round(inc1 .* vec1(2)),b),1);
                    doneList(x,y)=1;
                end
                
            end
        else
            doneCount=doneCount+1;
            %        disp 'already done'
        end
        %break
    end
    figure
    imshow(finalIm .* 1 + (baseline-snap0)./200)
    pause
    %imshow(doneList .* max1 + (baseline-snap0)./100)
    %title('final')
end