%plot excitatory types
load layer_array

sum_allvec_individual_vector5=transpose(sum(transpose(allvec_individual_vector5)));

lsum=zeros(1,100);
exc_labels={'L2','L23','L3','L4','L45','L5','L6','L6b'};
figure
count2=1;
max1=max(max(allvec_individual_vector5(exc_inds1,:)));

for count1=1:length(exc_labels)
    if (count1==1)
        inds=l2_inds1;
    end
    if (count1==2)
        inds=l23_inds1;
    end
    if (count1==3)
        inds=l3_inds1;
    end
    if (count1==4)
        inds=l4_inds1;
    end
    if (count1==5)
        inds=l45_inds1;
    end
    if (count1==6)
        inds=l5_inds1;
    end
    if (count1==7)
        inds=l6_inds1;
    end
    if (count1==8)
        inds=l6b_inds1;
    end
    if( sum(sum(allvec_individual_vector5(inds,:)))>0)
        subplot(3,4,count2);
        set(gca,'box','off')
        set(gca,'xcolor','w','ycolor','w','xtick',[],'ytick',[]);
        localmax=max(sum(allvec_individual_vector5(inds,:),1));
        max1=max(max1,localmax);

        plot_allen_layers(inds,allvec_individual_vector5,layer_array);

        set(gca,'box','off')
        lsum=lsum+squeeze(sum(allvec_individual_vector5(inds,:),1));
        title(exc_labels(count1))

      %  max1=100;
        plot([l12 l12],[0 max1],'--','color',[.7 .7 .7])
        plot([l23 l23],[0 max1],'--','color',[.7 .7 .7])
        plot([l34 l34],[0 max1],'--','color',[.7 .7 .7])
        plot([l45 l45],[0 max1],'--','color',[.7 .7 .7])
        plot([l56 l56],[0 max1],'--','color',[.7 .7 .7])
        count2=count2+1;
        
        drawnow
        refresh
    end
end
for count3=1:count2-1
    subplot(3,4,count3);
    if(max1)
        ylim([0 max1])
    end
end

for count1=1:length(exc_labels)
    subplot(3,4,count1);
    if(max1)
        ylim([0 max1])
    end
end

exc_inds2=[l2_inds1;l23_inds1;l3_inds1;l4_inds1;l45_inds1;l5_inds1;l6_inds1;l6b_inds1];
subplot(3,4,count2);
plot_allen_layers(exc_inds2,allvec_individual_vector5,layer_array)
title('Excitatory')

figure

count2=1;

max1=1e-9;
countMap=1;
allinds=[];
for count1=1:length(exc_labels)
    if (count1==1)
        inds=l2_inds1;
    end
    if (count1==2)
        inds=l23_inds1;
    end
    if (count1==3)
        inds=l3_inds1;
    end
    if (count1==4)
        inds=l4_inds1;
    end
    if (count1==5)
        inds=l45_inds1;
    end
    if (count1==6)
        inds=l5_inds1;
    end
    if (count1==7)
        inds=l6_inds1;
    end
    if (count1==8)
        inds=l6b_inds1;
    end

    newinds=[];
    for count3=1:length(inds)
        if( sum(sum(allvec_individual_vector5(inds(count3),:)))>0)
            newinds=[newinds;inds(count3)];
            allinds=[allinds;inds(count3)];
        end
    end

    if( sum(sum(allvec_individual_vector5(newinds,:)))>0)
        subplot(3,4,count2);
        localmax=max(sum(allvec_individual_vector5(newinds,:),1));
        max1=max(max1,localmax);

        %     plot_allen_layers(inds,allvec_individual_vector5,layer_array);

        plot(transpose(allvec_individual_vector5(newinds,:)));
        colororder(map(countMap:end,:))
        % set(gca,'FontSize', 2)
        % legend(ttype_labels(inds))
        countMap=countMap+length(newinds);
        m1=ttype_labels(newinds);
        for count0=1:length(m1)
            m1{count0}=strrep(m1{count0},'_',' ');
        end
        % legend(m1)

        lsum=lsum+squeeze(sum(allvec_individual_vector5(inds,:),1));
        title(exc_labels(count1))

        count2=count2+1;
        %

        max1=max(max(allvec_individual_vector5(inds,:)));
        hold on
        plot([l12 l12],[0 max1],'--','color',[.7 .7 .7])
        plot([l23 l23],[0 max1],'--','color',[.7 .7 .7])
        plot([l34 l34],[0 max1],'--','color',[.7 .7 .7])
        plot([l45 l45],[0 max1],'--','color',[.7 .7 .7])
        plot([l56 l56],[0 max1],'--','color',[.7 .7 .7])
        ylim([0 max1])
        drawnow
        refresh
    end
end

figure
plot(transpose(allvec_individual_vector5(allinds,:)), 'LineWidth', 2.5);
colororder(map(1:end,:))
m1=ttype_labels(allinds);
for count0=1:length(m1)
    m1{count0}=strrep(m1{count0},'_',' ');
end
lgd=legend(m1);
lgd.NumColumns=3
lgd.FontSize=2

