%plot sum allvec_individual_vector5(inds) with layers indicated
function [] = plot_allen_layers(inds,allvec_individual_vector5,layer_array)
l12=7;
l23=20;
l34=35;
l45=45;
l56=70;
vec_exc=sum(allvec_individual_vector5(inds,:),1);
vec_exc_parts=sum(transpose(allvec_individual_vector5(inds,:)));
vec_exc_parts=vec_exc_parts./sum(vec_exc_parts);
plot(vec_exc,'b','LineWidth',1.25)

hold on
sum_vec_exc=sum(vec_exc);
vecInds=zeros(1,100);
if(length(inds)>1)
    sum_layer_array=sum( repmat(transpose(vec_exc_parts),1,100) .* layer_array(inds,:));
else
    sum_layer_array=(layer_array(inds,:));
end

vec1=zeros(1,100);
vec1(1:l12)=sum_layer_array(1);%./length(vec1(1:l12));
vec1(l12+1:l23)=sum_layer_array(2);%./length(vec1(l12+1:l23));
vec1(l23+1:l34)=sum_layer_array(3);%./length(vec1(l23+1:l34));
vec1(l34+1:l45)=sum_layer_array(4);%./length(vec1(l34+1:l45));
vec1(l45+1:l56)=sum_layer_array(5);%./length(vec1(l45+1:l56));
vec1(l56+1:100)=sum_layer_array(6);%./length(vec1(l56+1:100));

vecInds=sum_vec_exc.*vec1./(sum(vec1)+1e-9);
plot(vecInds,'r')
hold on
max1=1.1 * max(max(vec_exc));
plot([l12 l12],[0 max1],'--','color',[.7 .7 .7])
plot([l23 l23],[0 max1],'--','color',[.7 .7 .7])
plot([l34 l34],[0 max1],'--','color',[.7 .7 .7])
plot([l45 l45],[0 max1],'--','color',[.7 .7 .7])
plot([l56 l56],[0 max1],'--','color',[.7 .7 .7])
ylim([0 max1])
end