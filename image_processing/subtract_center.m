%subtract the center from dlist
function [ dlist2 ] = subtract_center(maxd,center0,dlist2)
rightside=dlist2(maxd :length(dlist2));
rightsidex=(maxd:length(dlist2))-center0;


for count3=round(floor(center0)):-1: round(floor(center0-max(rightsidex)))+1;
    testDist=center0-count3;
    dlist2(count3)=dlist2(count3)-interp1(rightsidex,rightside,testDist);
end
dlist2(  round(floor(center0))+1: length(dlist2) )=0;

end

