% Daniel Keller
% copyright 2011-2014 EPFL
% mark cortex bottom of L4 boundary

clear
newdirlist=[];

dirlist=dir ('annotate_sagittal/annotate_*.dat');
%loop through all the text files in the directory

startoffset=1;
for count00=startoffset:length(dirlist)
    count00
    name1=dirlist(count00).name;
    
    k = strfind(name1, '.dat');
    newname='                    ';
    newname(1:k(1)-10)=name1(10:k(1)-1);
    newname=strcat(newname,'.jpg');
    
    inname=strcat('annotate_sagittal/',name1);
    outname=strcat('annotate_sagittal_L4/',strcat('L4_',name1));
    if (~(exist(outname)))
        name1
        
        g1=dlmread(inname);
        minx =(min(g1(:,1)));
        maxx =(max(g1(:,1)));
        miny =(min(g1(:,2)));
        maxy =(max(g1(:,2)));
        
        %calculate old layer 4 boundary
        
        left= 0.55 .* g1(1,:)+0.45 .* g1(2,:);
        
        right= 0.55 .* g1(4,:)+0.45 .* g1(3,:);
               
        leftx=round(left(1));
        lefty=round(left(2));
        rightx=round(right(1));
        righty=round(right(2));
        
        middle=(left+right)./2;
        middlex=round(middle(1));
        middley=round(middle(2));
        
        imname=strcat('./images_sagittal/',newname);
        if (~(exist(imname)))
            continue;
        end
        %load image
        im0=imread(imname);
        
        im0(lefty-5:lefty+5,leftx-5:leftx+5,1)=0;
        im0(righty-5:righty+5,rightx-5:rightx+5,1)=0;
        im0(middley-5:middley+5,middlex-5:middlex+5,1)=0;
        
        imshow(im0(round(miny):round(maxy),round(minx):round(maxx),:))
        [x1,y1,key] = ginput(1);
        
        if (length(x1)==0)
            continue
        end
        
        im0(round(miny)+y1-15:round(miny)+y1+15,round(minx)+x1-15:round(minx)+x1+15,2)=50;
        
        imshow(im0(round(miny):round(maxy),round(minx):round(maxx),:))
        [x2,y2] = ginput(1);
        
        if (length(x2)==0)
            continue
        end
        
        %now save the coordinates
        
        dlmwrite(outname,[round(minx)+x1,round(miny)+y1;round(minx)+x2,round(miny)+y2]);
        
        
    end
end
