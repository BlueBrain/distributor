function [allDots]=replaceAllDotsMap(allDots,allDotsMap)
[s1 s2]=size(allDots);
for count1=1:s1
    allDots(count1,5)=allDotsMap(allDots(count1,1));
end
end