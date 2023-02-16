
% Daniel Keller
% copyright 2011-2023 EPFL
% combine results of all of the individual image vectors

mkdir('vectors_collected_manual61')
mkdir('cumulative_collected_manual61')
mkdir('alldots_collected_manual61')
clear
cvlist=[];
cvlistintensity=[];
semlist=[];
newdirlist=[];
basedir='./result_vectors_sagittal61/';
dirlist=dir (strcat(basedir,'/*.dat'));
%loop through all the text files in the directory

for count00=1:length(dirlist)
    name1=dirlist(count00).name;
    k = strfind(name1, '.');
    newname='                    ';
    newname(1:k(1)-1)=name1(1:k(1)-1);
    newdirlist=[newdirlist;newname];
end

basedir2='./result_vectors_coronal_left61/';
dirlist2=dir (strcat(basedir2,'/*.dat'));
%loop through all the text files in the directory

for count00=1:length(dirlist2)
    name1=dirlist2(count00).name;
    k = strfind(name1, '.');
    newname='                    ';
    newname(1:k(1)-1)=name1(1:k(1)-1);
    
    newdirlist=[newdirlist;newname];
end

backgroundCompensationHigh =[0.0522 ,  11.3067]; %if less than 7 use this

allbaselines=[];
allmeans=[];
allmeansadjusted=[];

%go through unique entries
newdirlist=unique(newdirlist,'rows');
s0=size(newdirlist);
correction=1;
for count1=1:s0(1)
 
    count1
    totCellsList=[];
    totIntensityList=[];
    localcumulative=[];
    localAllDotsMap=zeros(1,100);
    localAllDots=[];
    searchfield=strtrim(newdirlist(count1,:));
    
    bigvector=zeros(100,81);
    
    outputname=strcat('./vectors_collected_manual61/',strcat(searchfield,'.dat'));
    
    outputnameCumulative=strcat('./cumulative_collected_manual61/cumul_',strcat(searchfield,'.dat'));
    outputnameAllDots=strcat('./alldots_collected_manual61/alldots_',strcat(searchfield,'.dat'));
    %skip the processing if it has been processed before
    %to really redo just remove the directory and start anew
    %    if exist(outputname)
    %        continue
    %end
    
    bigcount=0;
    dirlist= dir(strcat(strcat(basedir,searchfield),'*.dat'));
    for count00=1:length(dirlist)
        name1=dirlist(count00).name;
        target1=strcat(basedir,name1);
        if (exist(target1))
            f1=load(strcat(basedir,name1));
            
            f1(:,1)=f1(:,1).*(f1(:,1)>=0);
            
            target2=strcat('cumulative/cumul_',name1);
            if (exist(target2))
                cumulative=transpose(dlmread(target2));
                localcumulative=[localcumulative;cumulative];
            end
	    
            %alldots
            target2=strcat('alldots61/alldots_',name1);
            if (exist(target2))
                allDots=(dlmread(target2));
                allDotsMap=makeAllDotsMap(allDots);
                localAllDotsMap=localAllDotsMap+allDotsMap;           
                localAllDots=[localAllDots;allDots];
                localAllDots=replaceAllDotsMap(localAllDots,localAllDotsMap);
            end
            
            %load up the baselines
            target2=strcat('baselines61/base_',name1);
            if (exist(target2))
                
                baselines=dlmread(target2);
                baselineDiff=baselines(2)-baselines(1);
                
                %    if (baselineDiff<-1)
                %        name1
                %        pause
                %    end
                
                %if the image is saturated and low density observed then do not use
                if (baselineDiff>-5)
                    allbaselines=[allbaselines;(baselines(2)-baselines(1))];
                    allmeans=[allmeans;mean(f1(:,1))];
                    
                    correction=1;%backgroundCompensationHigh(2)+backgroundCompensationHigh(1).*baselineDiff;
                    if (correction<0)                        
                        [correction baselineDiff]
                        pause
                    end
                    f1(:,1)=f1(:,1)./correction;
                    %f1(:,1)=f1(:,1).*(f1(:,1)>=0);
                    %now correct cell counts
                    %f1(:,2:41)=(f1(:,2:41) ./(f1(:,42:81)+1e-9)-correction).*f1(:,42:81);
                    f1(:,2:41)=f1(:,2:41) ./(correction);
                    allmeansadjusted=[allmeansadjusted;mean(f1(:,1))];
                    f1=f1.*(f1>0);
                else
                    disp 'here'
                    
                    continue
                end
            else
                target2
                disp 'baseline not found'
                pause
            end
        end
        
        totCellsList=[totCellsList;sum(sum( f1(:,42:81)))];
        totIntensityList=[totIntensityList;sum(sum( f1(:,1)))];
        bigvector=bigvector+f1;
        bigcount=bigcount+1;
    end
    
    dirlist= dir(strcat(strcat(basedir2,searchfield),'*.dat'));
    for count00=1:length(dirlist)
        name1=dirlist(count00).name;
        target1=strcat(basedir2,name1);
        if (exist(target1))
            f1=load(strcat(basedir2,name1));
            f1(:,1)=f1(:,1).*(f1(:,1)>=0);
            
            
            target2=strcat('cumulative/cumul_',name1);
            if (exist(target2))
                cumulative=transpose(dlmread(target2));
                localcumulative=[localcumulative;cumulative];
            end
            
            target2=strcat('alldots61/alldots_',name1);
            if (exist(target2))
                allDots=(dlmread(target2));
                allDotsMap=makeAllDotsMap(allDots);
                localAllDotsMap=localAllDotsMap+allDotsMap;           
                localAllDots=[localAllDots;allDots];
                localAllDots=replaceAllDotsMap(localAllDots,localAllDotsMap);
            end
           
            
            %load up the baselines
            target2=strcat('baselines61/base_',name1);
            if (exist(target2))
                baselines=dlmread(target2);
                baselineDiff=baselines(2)-baselines(1);
                
                %if the image is saturated and low density observed then do not use
                if (baselineDiff>-5)
                    allbaselines=[allbaselines;(baselines(2)-baselines(1))];
                    allmeans=[allmeans;mean(f1(:,1))];
                    
                    if (correction<0)
                        [correction baselineDiff]
                        pause
                    end
                    f1(:,1)=f1(:,1)./correction;
                    %f1(:,1)=f1(:,1).*(f1(:,1)>=0);
                    %now correct cell counts
                    %f1(:,2:41)=(f1(:,2:41) ./(f1(:,42:81)+1e-9)-correction).*f1(:,42:81);
                    f1(:,2:41)=f1(:,2:41) ./(correction);
                    allmeansadjusted=[allmeansadjusted;mean(f1(:,1))];
                    f1=f1.*(f1>0);
                else
                    
                    continue
                end
            else
                disp 'baseline not found'
                pause
            end
        end
        
        totCellsList=[totCellsList;sum(sum( f1(:,42:81)))];
        totIntensityList=[totIntensityList;sum(sum( f1(:,1)))];
        bigvector=bigvector+f1;
        bigcount=bigcount+1;
    end
    
    
    val=0;
    val2=0;
    val3=0;
    if (length(totCellsList)>1)
        val=std(totCellsList)/(mean(totCellsList)+1e-9);
        val2=std(totCellsList)/(mean(totCellsList)+1e-9)/sqrt(length(totCellsList)-1 );
        val3=std(totIntensityList)/(mean(totIntensityList)+1e-9)/sqrt(length(totIntensityList)-1 );
    end
    
    cvlist=[cvlist;val];
    semlist=[semlist;val2,val3];
    %if (val<1.2) %only take images that are not crazy different
    
    if (bigcount>0)
        dlmwrite(outputname,bigvector./bigcount);
    end
    
 
    [s1 s2]=size(localcumulative);
    if (s1>1)
        localcumulative=mean(localcumulative);
    end
    dlmwrite(outputnameCumulative,localcumulative);
    dlmwrite(outputnameAllDots,localAllDots);
    
end

figure
plot(allbaselines,allmeans,'.')
hold on
newallbaselines=[];
newallmeans=[];
for count1=min(allbaselines):max(allbaselines)
    newallbaselines=[newallbaselines;count1];
    newallmeans=[newallmeans;mean(allmeans(find( (allbaselines>=count1).*(allbaselines<count1+1)     )))];
end
plot(newallbaselines,newallmeans,'r');
xlabel('slice background - image background')
ylabel('average image intensity')
dlmwrite('backgroundCompensation.dat',[newallbaselines,newallmeans]);

y=polyfit(newallbaselines( find(newallbaselines<40 )),newallmeans( find(newallbaselines<40)) ,1)
plot([2:7],y(2)+y(1).*[2:7],'y')
%figure
%plot(allbaselines,allmeansadjusted,'.')
%xlabel('slice background - image background')
%ylabel('adj. average image intensity')

%pause
%coeffs = polyfit(allbaselines,allmeans, 1);

%write the cvs
aveCV=mean(cvlist(find(cvlist>0)));
%cvlist=cvlist+aveCV .*(cvlist==0);
semlist(:,1)=semlist(:,1)+aveCV .*(semlist(:,1)==0);
semlist(:,2)=semlist(:,2)+aveCV .*(semlist(:,2)==0);

for count1=1:s0(1)
    searchfield=strtrim(newdirlist(count1,:));
    outputname=strcat('./vectors_collected_manual61/',strcat(searchfield,'.cv'));
    dlmwrite(outputname,semlist(count1,:));
end

if(0)
 cd ../SST_use_Nissl
 clear
 load soo61
 [cvlist,cvlist2,allimdagevectors] = load_images_withsizes56(names1);
 save soo61
 end