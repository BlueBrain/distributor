function [gtot2,kernel] = analyze_gtot2(goodinds)


basedir='../image_processing/result_vectors_sagittal_nissl61/';
dirlist=dir (strcat(basedir,'/*.dat'));

length(dirlist)
startdiroffset=1;

%load('scale_ish_expression')

allgtot=[];
gtot=zeros(100,81);
gsum=0;
for count00=startdiroffset:length(dirlist)
    if (count00>300)
        break
    end
    name1=dirlist(count00).name
    
    filestr1=strcat(strcat(basedir,name1));
    
    if (~(exist(filestr1,'file')))
        continue
    end
    
    g1=dlmread(filestr1);
    gtot=gtot+g1;
    g1(:,1)=g1(:,1);%./transpose(scale_ish_expression);
    allgtot=[allgtot;g1];
    gsum=gsum+1;
    %plot(g1(:,4))
    %pause
end
gtot=gtot./gsum;
%find the scale factor
%how many mm3 is a pixel?
%I have 1300 pixels for the thickness of cortex
%DeFelipe 2011 gives 1210 um
%slices are 25 um thick
voxelVolume=0.025 * (1210/1300 *0.001).^2; %in mm3
gtot(:,2:81)=gtot(:,2:81)./voxelVolume;
allgtot(:,2:81)=allgtot(:,2:81)./voxelVolume;
kernel=lsqnonneg((allgtot(:,42:81)),allgtot(:,1));
%figure
%plot(kernel,'r')
%dsffdsfd

x=0.*gtot(:,2:21);
y=0.*gtot(:,2:21);
gtot2=gtot;
for count1=1:100
    for count2=1:20
        x(count1,count2)=count1;
        y(count1,count2)=count2;
    end
end

    if(1)
figure
surf(x,y,gtot2(:,42:61))
ylabel('Cell Area ( x10 pixels)')
xlabel('Fraction Cortical Depth')
az = 0;
el = 90;
view(az, el);
set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',[screenposition(3:4)]);
print -dpdf -painters figure1

end

load('dot_expression_matrix')
voxelVolume=0.025 * (1210/1300 *0.001).^2; %in mm3
dot_expression_matrix=dot_expression_matrix./voxelVolume;
s1=squeeze(mean(dot_expression_matrix(goodinds,:,:),1));
if(0)
figure
surf(x,y,s1(:,1:20))
ylabel('Cell Area ( x20 pixels)')
xlabel('Fraction Cortical Depth')
az = 0;
el = 90;
view(az, el);
set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',[screenposition(3:4)]);
print -dpdf -painters figure2
end

if(0)
figure
size(sum(transpose(gtot2(:,42:81))))
plot(sum(transpose(gtot2(:,42:81))))
hold on
size(sum(transpose(s1)))
plot(sum(transpose(s1)))
ylabel('Unadjusted Density (cells/mm^3)')
xlabel('Fraction Cortical Depth')
set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',[screenposition(3:4)]);
print -dpdf -painters figure3
end


%clear('s1')
%clear('dot_expression_matrix')

figure
plot(sum(s1))
hold on
plot(sum((gtot2(:,42:81))))

figure
plot(sum(s1),sum((gtot2(:,42:81))))
xlabel('dot expression')
ylabel('nissl')

end
