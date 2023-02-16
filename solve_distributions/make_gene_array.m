clear

%reads in image counts and turns them into csv file
if(1)
    T = readtable('medians_combined.csv');
    gene_names=T.feature;
    T.feature=[];
    ttype_labels=T.Properties.VariableNames;
    save('ttype_labels','ttype_labels')
    clusters=T.Variables;
    save('gene_names','gene_names');
else
    load('gene_names')
end


fileName='./transcriptome_process3/output.txt';
pat = '(?<gene>\w+)\s+(?<image>\w+)\s';
titles = regexp( fileread(fileName), pat, 'names');

dot_expression_matrix=nan(length(gene_names),100,40);
for count1=1:length(gene_names)
    filename=strcat(strcat('../image_processing/vectors_collected_manual61/',gene_names(count1)),'.dat');
    filename=filename{1};
    if ((exist(filename)<=0))
        found=0;
        for count2=1:length(titles)
            if (strcmp(titles(count2).gene,gene_names(count1)))
                filename=strcat(strcat('../image_processing/vectors_collected_manual61/',titles(count2).image),'.dat');
                disp ('reached')
                exist(filename)
                %disp(filename)
                %huhuhu
                if ((exist(filename)<=0))
                    disp('breaking')
                    break
                end
                found=1;                
                %disp('good')
                
                break
            end
        end
        if(found==0)
            continue
        end
        
    end
    m1=dlmread(filename);
    dot_expression_matrix(count1,:,:)=m1(:,42:81);
end
save('dot_expression_matrix','dot_expression_matrix')

%fileID = fopen('./transcriptome_process3/output.txt','r');
%A = fscanf(fileID,'%s %s')
fileName='./transcriptome_process3/output.txt';
pat = '(?<gene>\w+)\s+(?<image>\w+)\s';
titles = regexp( fileread(fileName), pat, 'names');

expression_matrix=nan(length(gene_names),100,40);
for count1=1:length(gene_names)
    filename=strcat(strcat('../image_processing/vectors_collected_manual61/',gene_names(count1)),'.dat');
    filename=filename{1};
    if ((exist(filename)<=0))
        found=0;
        for count2=1:length(titles)
            if (strcmp(titles(count2).gene,gene_names(count1)))
                filename=strcat(strcat('../image_processing/vectors_collected_manual61/',titles(count2).image),'.dat');
                disp ('reached')
                exist(filename)
                %disp(filename)
                %huhuhu
                if ((exist(filename)<=0))
                    disp('breaking')
                    break
                end
                found=1;                
                %disp('good')
                
                break
            end
        end
        if(found==0)
            continue
        end
        
    end
    m1=dlmread(filename);
    expression_matrix(count1,:,:)=m1(:,2:41);
end
save('expression_matrix','expression_matrix')

