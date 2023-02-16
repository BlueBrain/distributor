%perform image processing on coronal slices
function annotate_translate_coronal61(startdiroffset)
dirlist=dir ('./annotate_coronal_left/a*.dat'); %all neur marker

%check if directory exists, otherwise create
mkdir('result_vectors_coronal_left61')
mkdir('annotate_coronal_L4')
mkdir('baselines61')
mkdir('cumulative61')
mkdir('alldots61')
mkdir('heightdots61')

for count00=startdiroffset:min(startdiroffset+20000, length(dirlist))
    %for count00=startdiroffset:-1:1
    
    name1=dirlist(count00).name
    str1=name1(10:length(name1)-4);
    
    filestr1=strcat(strcat('images_coronal/',str1),'.jpg');
    if (~(exist(filestr1,'file')))
        disp 'no jpeg found'
        continue
    end
    
    %take this part out if you don't want to skip
    mapname=strcat(strcat('result_vectors_coronal_left61/',str1),'.dat');
    if ((exist(mapname,'file')))
        disp 'skipping1'
        continue
    end

    %now read in the marked points
    pointsL4=[];    

    filestr2=strcat('annotate_coronal_L4/',strcat('L4_',name1));
    if (exist(filestr2))
      disp 'reading L4'
      pointsL4=dlmread(filestr2);
    end


    %now read in the marked points
    filestr2=strcat('annotate_coronal_left/',name1);
    points=dlmread(filestr2);
    
    layerThickness=get_layer_thickness(points);
    areaScale=(1299/layerThickness).^2;
    
    p1=(points(1,:));
    p2=(points(2,:));
    p3=(points(3,:));
    p4=(points(4,:));
    %now do a simple markup 
    len1=round(norm(points(1,:)-points(4,:)));
    surflist=[];
    innervec=[];
    for count1=0:len1
        surflist=[surflist;count1/len1 .*p4+(1.0-count1/(len1)) .*p1];
        innervec=[innervec;count1/len1 .*p3+(1.0-count1/(len1)) .*p2];
    end
    dir1=p2+p3-p1-p4;
    width=norm(dir1);
    dir1=dir1./norm(dir1);
    
    dir2=p4+p3-p1-p2;
    dir2=dir2./norm(dir2);
    
    base1=(p1+p2)./2-dir1 .*50;
    %now call recognition on this
    
    try
        newim=imread(filestr1);
    catch
        filestr1
        disp 'skipping2'
        continue
    end
    
    newim=uint16(newim(:,:,3));
    
    if(sum(sum(newim))==0 )
        disp 'is zero'
    end
    
    maxlength=5;%pixel width
    surflist=move_bar2(newim,surflist,dir1,dir2,maxlength);
    
    length1=width;
        
    [s r]=size(newim);
   
    %make a snapshot and run tests on it to determine threshold
    snapMin1=round(min([p1(1),p2(1),p3(1),p4(1)]));
    snapMin2=round(min([p1(2),p2(2),p3(2),p4(2)]));
    snapMax1=round(max([p1(1),p2(1),p3(1),p4(1)]));
    snapMax2=round(max([p1(2),p2(2),p3(2),p4(2)]));
    snap1=newim(snapMin2:snapMax2,snapMin1:snapMax1);
    
    %take next few lines out    
    snap1=double(newim(snapMin2:snapMax2,snapMin1:snapMax1));

    [m n]=size(snap1);
    radius1=sparse(m,n);
    
    %find an appropriate threshold and adjust the image
    
    base=[snapMin1-1,snapMin2-1];

    [baseline,background,snap1]=get_threshold50(snap1,p1-base,p2-base,p3-base,p4-base);
    newim(snapMin2:snapMax2,snapMin1:snapMax1)=snap1 .*(snap1<=baseline)+ baseline .*(snap1>baseline) ;

    map1=baseline-map200b(newim ,surflist,innervec,dir1,baseline);

    mapname0=strcat(strcat('baselines61/base_',str1),'.dat');
    dlmwrite(mapname0,[baseline,background]);
    
    newim(snapMin2:snapMax2,snapMin1:snapMax1)=snap1 .*(snap1>=0);
    map1=baseline-map200b(newim ,surflist,innervec,dir1,baseline);

    newim(snapMin2:snapMax2,snapMin1:snapMax1)=snap1;
    [a b]=size(snap1);
    endval2=baseline;
    snap0=snap1;
    [r1 r2]=size(newim);

    [map2,map4,map7,heightMask,cumulativeImageExpression,allDotsArray]=get_maps_from_endval61_nissl(baseline,snap0,baseline+1,surflist,innervec,dir1,snapMin1,snapMin2,snapMax1,snapMax2,areaScale,pointsL4,r1,r2);

    dlmwrite(mapname,[map1,map4,map7]);
    
    mapname2=strcat(strcat('cumulative61/cumul_',str1),'.dat');
    dlmwrite(mapname2,cumulativeImageExpression);
    
    mapname3=strcat(strcat('alldots61/alldots_',str1),'.dat');
    dlmwrite(mapname3,allDotsArray);

    
    
end

exit

