function [newinds] = correct_inds(uncorrected_inds,ttype_labels,clusters0,goodinds)
newinds=[];
for count1=1:length(uncorrected_inds)
    target=strcat('x',num2str(uncorrected_inds(count1)));
     target=strcat(target,'_');
     for count2=1:length(ttype_labels)
        if (sum(clusters0(goodinds,count2))>0)
            if(length(strfind(ttype_labels{count2},target)))
                newinds=[newinds;count2];
            end
        end
    end

end

