% Daniel Keller
% copyright 2011-2023 EPFL
% mark sagittal slices

  clear
  newdirlist=[];
  mkdir('annotate_sagittal')
  dirlist=dir ('./images_sagittal/*.jpg');
  %loop through all the text files in the directory
  

  for count00=length(dirlist):-1:1
    name1=dirlist(count00).name;
    name1=strrep(name1,'?','*');
    k = strfind(name1, '.jpg');
    newname='                    ';
    newname(1:k(1)-1)=name1(1:k(1)-1);
    newname=strcat(newname,'.dat');
    outname=strcat('annotate_sagittal/annotate_',newname);
    if (~(exist(outname)))
        name1
        
        imname=strcat('./images_sagittal/',name1);
        imname=strrep(imname,'?','*');
         if (~(exist(imname)))
	   disp 'image does not exist'
       imname
     %  pause
             continue;
         end
      %load image
	 im0=imread(imname);
         [s r t]=size(im0);
         im1=im0(1:round(s/2),1:round(r/2),3);
         
 %        im1=im0(1:round(2*s/3),1:round(2*r/3),3);
   %      im1=im0;
         imshow(im1)
         [x1,y1] = ginput(1)
         y1=max(y1,1);
      	 im1(y1-15:y1+15,x1-15:x1+15)=0;
         imshow(im1)
         [x2,y2] = ginput(1)
         y2=max(y2,1);
         im1(y2-15:y2+15,x2-15:x2+15)=0;
         imshow(im1)
	 [x3,y3] = ginput(1)
         y3=max(y3,1);
      	 im1(y3-15:y3+15,x3-15:x3+15)=0;
         imshow(im1)
	 [x4,y4] = ginput(1)
         y4=max(y4,1);
      	 im1(y4-15:y4+15,x4-15:x4+15)=0;
         imshow(im1)
	 dlmwrite(outname,[x1,y1;x2,y2;x3,y3;x4,y4]);
	 
    end
  end
