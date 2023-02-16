%prepare and assemble the data
clear
close all

make_gene_array

    T = readtable('medians_combined.csv');
    gene_names=T.feature;
    T.feature=[];
    ttype_labels=T.Properties.VariableNames;
    clusters=T.Variables;

    allowed_ttypes=transpose(ttype_labels);

    [~,txt,raw] = xlsread('./expression_profiles_combined.xls',1);%,'A1:CW45768');

    raw(:,1)=[];
    ish_expression=cell2mat(raw);

    [~,~,raw] = xlsread('./expression_profiles_combined.xls',2,'A1:CW45768');
    raw(:,1)=[];
    ish_expression_cv=cell2mat(raw);
    expression_cv=(sum(transpose(ish_expression_cv(:,1:99))));
    save('expression_cv','expression_cv');

    clear('txt')
    clear('raw')

    [num_genes, num_clusters]=size(clusters);
    mean_clusters=mean(transpose(clusters));
   
    goodinds=find( (~(sum(transpose(isnan(clusters)))+sum(transpose(isnan(ish_expression))))).*(mean_clusters>0) );

    %part II

    [gtot2,kernels]=analyze_gtot2(goodinds);



scaleVec=ones(num_genes,1);
ish_expression0=ish_expression;
clusters0=clusters;


load('expression_matrix')
load('dot_expression_matrix')
voxelVolume=0.025 * (1210/1300 *0.001).^2; %in mm3
dot_expression_matrix=dot_expression_matrix./voxelVolume;


%check for folded over l1
if(0)
    goodindsLower=[];
    average_gtot=sum(transpose(gtot2(:,42:81)));
    mean_gtot_l1=mean(average_gtot(1:7));
    for count2=length(goodinds):-1:1

        local_sum=sum(transpose(squeeze(dot_expression_matrix(goodinds(count2),:,:))));
        % plot(local_sum,'r')
        if (mean(local_sum(1:7))<1.2*mean_gtot_l1)

            goodindsLower=[goodindsLower,goodinds(count2)];
        end
    end
    clear('local_sum')
end

close all
%plot intensity expression
s1=squeeze(mean(expression_matrix(goodinds,:,:),1));
figure
plot(sum(transpose(s1)))
title('Mean Intensity Expression')
xlabel('Cortical Depth Bin')
ylabel('Intensity')
set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',[screenposition(3:4)]);
%print -dpdf -painters figure1



[a b]=size(clusters);
placed_types=[1:b];

%plot gene expression versus ish
figure

mean_ish_expression=mean(transpose(ish_expression));
plot(mean_clusters(goodinds),mean_ish_expression(goodinds),'k.')
p=polyfit(mean_clusters(goodinds),mean_ish_expression(goodinds),1)
%q=polyfit(mean_ish_expression(goodinds),mean_clusters(goodinds),1)
slope1=(mean_ish_expression(goodinds)/mean_clusters(goodinds));
%mdl=fitlm(mean_clusters(goodinds),mean_ish_expression(goodinds),'linear','RobustOpts','off')

hold on
plot(mean_clusters(goodinds),mean_clusters(goodinds) .*p(1) +p(2),'b')
plot([0 14],[p(2) p(2)],'r--')

xlabel('gene expression level')
ylabel('ISH intensity')
title('Comparison of Transcriptome versus ISH')
set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',[screenposition(3:4)]);
%print -dpdf -painters figure4


num_ttypes=length(ttype_labels);




goodtypes=find(sum(clusters0(goodinds,:))>0);

scaleVec=scaleVec./sum(scaleVec);
s1=squeeze(mean(expression_matrix(goodinds,:,:),1));




disp 'Intensity eval'
errlist=[];
%    scaleVec=ones(num_genes,1)./length(goodinds);
maxScale=1/num_ttypes;
ish_Scale=ones(num_genes,1);
baseline_Scale=ones(num_genes,1);

options = optimset('Display','notify','TolX',1e-12);
numReps=2;
olderr=1e27;
%load('scaleVec')
%    for count00=1:numReps

%scaleVec(goodinds)=1;
%scaleVec0=scaleVec;


num_placed_types=length(placed_types);

mean_expression_matrix=(mean(expression_matrix(goodinds,:,:),1));
mean_dot_expression_matrix=(mean(dot_expression_matrix(goodinds,:,:),1));
expression_ratio=mean(mean_expression_matrix./(mean_dot_expression_matrix+1e-9));

%plot expression ratio
figure
expression_ratio=squeeze(expression_ratio./expression_ratio(1));
plot(expression_ratio);
title('Expression Ratio')
xlabel('Cortical Depth Bin')
ylabel('Ratio')
set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',[screenposition(3:4)]);
%print -dpdf -painters figure10

mean_dot_expression_matrix=squeeze(mean_dot_expression_matrix);

%clear('clusters')
%clear('dot_expression_matrix')
%clear('expression_matrix')

load('dot_expression_matrix')
load('expression_matrix')

%make initial version of best_kernel
num_placed_types=length(placed_types);
allvec_individual= zeros(num_placed_types,100,40);
best_baseline=0;
%
%load goodindsLower


type_clusters00=([mean(transpose(clusters0(:,1:92)));mean(transpose(clusters0(:,93:270)));mean(transpose(clusters0(:,272:273)));mean(transpose(clusters0(:,274:279)));mean(transpose(clusters0(:,280:284)));mean(transpose(clusters0(:,285:287)));(transpose(clusters0(:,288)));mean(transpose(clusters0(:,289:290)))] );
type_clusters0= type_clusters00>0;
type_clusters=1./(sum(type_clusters0));
type_clusters(find(isnan(type_clusters)))=0;



OPTIONS.TolX=1e-12;
OPTIONS.TolFun=1e-12
OPTIONS.MaxIter=150000000000000;
%    options = optimset('Display','notify','TolX',1e-12)
options = optimset('Display','notify','TolX',1e-15)

if(0)
mean_expression_matrix=zeros(100,40);
scale1=mean(mean(mean(dot_expression_matrix(goodinds,:,:))));
for count1=1:100
    for count2=1:40
        mean_dot=mean(dot_expression_matrix(goodinds,count1,count2));
        mean_expression=mean(expression_matrix(goodinds,count1,count2));
 
        %        expression_matrix(goodinds,count1,count2)=expression_matrix(goodinds,count1,count2).*mean_dot./(mean_expression+1e-11)./scale1;

        mean_expression_matrix(count1,count2)=mean(expression_matrix(goodinds,count1,count2));
    end
end
end

%expr_sum=zeros(num_genes,1);
%for count1=1:length(goodinds)
%    expr_sum(goodinds(count1))=sum(sum(squeeze(expression_matrix(goodinds(count1),:,2:40))));
%end
if(0)
figure
plot(mean_clusters(goodinds),expr_sum(goodinds),'.');
xlabel('gene expression level')
ylabel('ISH intensity')
title('Comparison of Transcriptome versus bin ISH')
set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',[screenposition(3:4)]);
print -dpdf -painters figure11
end

%expression_matrix_processed=expression_matrix;
%save('expression_matrix.mat''expression_matrix')
disp 'fixing expression matrix'
%save state1
save('goodinds','goodinds')
save('clusters0','clusters0')
save('gtot2','gtot2')
%normalize_expression28

