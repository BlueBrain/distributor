%find layer information in metadata
load ttype_labels
m1=readtable("metadata.csv");

l12=7;
l23=20;
l34=35;
l45=45;
l56=70;

w1=l12;
w2=l23-l12-1;
w3=l34-l23-1;
w4=l45-l34-1;
w5=l56-l45-1;
w6=100-l56-1;


layer_array=zeros(length(ttype_labels),6);
for count1=1:length(ttype_labels)
%    count1
    target=ttype_labels{count1};
    target=strrep(target,'/','_');
    for count2=1:length(m1.cluster_label)
        local=strcat('x',strrep(m1.cluster_label{count2},' ',''));
        local=strrep(local,'/','_');
        local=strrep(local,'-','_');
        if(strcmp(target,local))
                   
            %if(length(m1.cortical_layer_label{count2}))
                disp 'found'
                m1.cortical_layer_label{count2}
 %layer_array(count1,:)
  %              pause
                %continue
   %         end
            if(strcmp(m1.cortical_layer_label{count2},'L4/5/6'))
                layer_array(count1,[4])=layer_array(count1,[4])+w4./(w4+w5+w6);
                layer_array(count1,[5])=layer_array(count1,[5])+w5./(w4+w5+w6);
                layer_array(count1,[6])=layer_array(count1,[6])+w6./(w4+w5+w6);

                continue
            end
            if(strcmp(m1.cortical_layer_label{count2},'L2/3'))
                layer_array(count1,[2])=layer_array(count1,[2])+w2./(w2+w3);
                layer_array(count1,[3])=layer_array(count1,[3])+w3./(w2+w3);

                continue
            end
            if(strcmp(m1.cortical_layer_label{count2},'L4'))
                layer_array(count1,[4])=layer_array(count1,[4])+1;
                continue
            end
            if(strcmp(m1.cortical_layer_label{count2},'L5'))
                layer_array(count1,[5])=layer_array(count1,[5])+1;
                continue
            end
            if(strcmp(m1.cortical_layer_label{count2},'L1'))
                layer_array(count1,[1])=layer_array(count1,[1])+1;
                continue
            end
            if(strcmp(m1.cortical_layer_label{count2},'L6'))
                layer_array(count1,[6])=layer_array(count1,[6])+1;
                continue
            end
            if(strcmp(m1.cortical_layer_label{count2},'L5/6'))
                layer_array(count1,[5])=layer_array(count1,[5])+w5./(w5+w6);
                layer_array(count1,[6])=layer_array(count1,[6])+w6./(w5+w6);

                continue
            end
            if(strcmp(m1.cortical_layer_label{count2},'L1/2/3/4'))
                layer_array(count1,[1])=layer_array(count1,[1])+w1./(w1+w2+w3+w4);
                layer_array(count1,[2])=layer_array(count1,[2])+w2./(w1+w2+w3+w4);
                layer_array(count1,[3])=layer_array(count1,[3])+w3./(w1+w2+w3+w4);
                layer_array(count1,[4])=layer_array(count1,[4])+w4./(w1+w2+w3+w4);

                continue
            end
            if(strcmp(m1.cortical_layer_label{count2},'All'))
                layer_array(count1,[1])=layer_array(count1,[1])+w1./(100);
                layer_array(count1,[2])=layer_array(count1,[2])+w2./(100);
                layer_array(count1,[3])=layer_array(count1,[3])+w3./(100);
                layer_array(count1,[4])=layer_array(count1,[4])+w4./(100);
                layer_array(count1,[5])=layer_array(count1,[5])+w5./(100);
                layer_array(count1,[6])=layer_array(count1,[6])+w6./(100);

                continue
            end
            if(strcmp(m1.cortical_layer_label{count2},'L1/2/3'))
                layer_array(count1,[1])=layer_array(count1,[1])+w1./(w1+w2+w3);
                layer_array(count1,[2])=layer_array(count1,[2])+w2./(w1+w2+w3);
                layer_array(count1,[3])=layer_array(count1,[3])+w3./(w1+w2+w3);

                continue
            end
            if(strcmp(m1.cortical_layer_label{count2},'L4/5'))
                layer_array(count1,[4])=layer_array(count1,[4])+w4./(w4+w5);
                layer_array(count1,[5])=layer_array(count1,[5])+w5./(w4+w5);
                continue
            end
            if(strcmp(m1.cortical_layer_label{count2},'L2/3/4'))
                layer_array(count1,[2])=layer_array(count1,[2])+w2./(w2+w3+w4);
                layer_array(count1,[3])=layer_array(count1,[3])+w3./(w2+w3+w4);
                layer_array(count1,[4])=layer_array(count1,[4])+w4./(w2+w3+w4);

                continue
            end
            if(strcmp(m1.cortical_layer_label{count2},'L2'))
                layer_array(count1,[2])=layer_array(count1,[2])+1;
                continue
            end
            if(strcmp(m1.cortical_layer_label{count2},'L6b'))
                layer_array(count1,[6])=layer_array(count1,[6])+1;
                continue
            end
            if(strcmp(m1.cortical_layer_label{count2},'L2/3/5'))
                layer_array(count1,[2])=layer_array(count1,[2])+w2./(w2+w3+w5);
                layer_array(count1,[3])=layer_array(count1,[3])+w3./(w2+w3+w5);
                layer_array(count1,[5])=layer_array(count1,[5])+w5./(w2+w3+w5);

                continue
            end
            if(strcmp(m1.cortical_layer_label{count2},'L6a'))
                layer_array(count1,[6])=layer_array(count1,[6])+1;
                continue
            end
            
            %m1.cortical_layer_label{count2}
        end
    end
end
save('layer_array','layer_array');