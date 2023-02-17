% Daniel Keller
% copyright 2011-2014 EPFL
% mark nissl stains

  clear
  newdirlist=[];
  mkdir('annotate_sagittal_nissl')
  dirlist=dir ('./images_sagittal_nissl/*.jpg');
  %loop through all the text files in the directory
  

  for count00=1:length(dirlist)
    name1=dirlist(count00).name;

    k = strfind(name1, '.jpg');
    newname='                    ';
    newname(1:k(1)-1)=name1(1:k(1)-1);
    newname=strcat(newname,'.dat');
    outname=strcat('annotate_sagittal_nissl/annotate_',newname);
    if (~(exist(outname)))
        name1
      %load image
	 im0=imread(strcat('./images_sagittal_nissl/',name1));
         [s r t]=size(im0);
         im1=im0(1:round(s/2),1:round(r/2),3);
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
