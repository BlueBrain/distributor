
%fit a vector to newsurflist
function [newdir,newdir2]=fit_vector(newsurflist)

p=polyfit(newsurflist(:,1),newsurflist(:,2),1);
angle=pi/2;
rotmatrix=[cos(angle),-sin(angle);sin(angle),cos(angle)];
point1=[1 ;p(1)];
newdir2=point1./norm(point1);
newpoint =transpose(rotmatrix * point1);
newdir=newpoint./norm(newpoint) ;

