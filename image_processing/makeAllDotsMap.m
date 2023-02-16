%make alldotsmap
function [allDotsMap]=makeAllDotsMap(allDots)
allDotsMap=zeros(1,100);
[s1 s2]=size(allDots);
for count1=1:s1
    allDotsMap((allDots(count1,1)))=allDots(count1,5);
end
%if it was not present then assume a value based on median
median1=median(allDotsMap(allDotsMap>0));
allDotsMap(allDotsMap==0)=median1;
end