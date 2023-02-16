%reads in image counts and turns them into csv file
[~,txt,~] = xlsread('medians.csv');
gene_names=txt(:,1);
gene_names(1)=[];
clear('txt')

expression_matrix=NaN(length(gene_names),100);
for count1=1:length(gene_names)
    filename=strcat(strcat('./vectors_collected_manual61/',gene_names(count1)),'.dat');
    if (exist(filename{1}))
        count1
        m1=dlmread(filename{1});
        expression_matrix(count1,:)=m1(:,1);
    end
end
T=num2cell([zeros(length(gene_names),1),expression_matrix]);
T(:,1)=gene_names;
filename = 'expression_profiles_combined.xls';
writecell(T,filename,'Sheet',1,'Range','A1');
filename = 'expression_profiles_combined.xlsx';
writecell(T,filename,'Sheet',1,'Range','A1');

expression_matrix=NaN(length(gene_names),100);
for count1=1:length(gene_names)
    filename=strcat(strcat('./vectors_collected_manual61/',gene_names(count1)),'.cv');
    if (exist(filename{1}))
        count1
        m1=dlmread(filename{1});
        expression_matrix(count1,:)=m1(:,1);
    end
end
T=num2cell([zeros(length(gene_names),1),expression_matrix]);
T(:,1)=gene_names;
filename = 'expression_profiles_combined.xls';
writecell(T,filename,'Sheet',2,'Range','A1');
filename = 'expression_profiles_combined.xlsx';
writecell(T,filename,'Sheet',2,'Range','A1');
