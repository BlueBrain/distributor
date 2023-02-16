function [bar_cr,bar_som,bar_pv]= get_meyer_percentages()
%hold on
som=xlsread('meyer_som_percentage.csv');
%plot(1-som(:,2),som(:,1))
som_interp=interp1(1-som(:,2),som(:,1),(0.01:0.01:1));
som_interp(99)=som_interp(98);
som_interp(100)=som_interp(98);
pv=xlsread('meyer_pv_percentage.csv');
pv_interp=interp1(1-pv(:,2),pv(:,1),(0.01:0.01:1));
pv_interp(99)=pv_interp(98);
pv_interp(100)=pv_interp(98);
%plot(1-pv(:,2),pv(:,1))
cr=xlsread('meyer_cr_percentage.csv');
%plot(1-cr(:,2),cr(:,1))
cr_interp=interp1(1-cr(:,2),cr(:,1),(0.01:0.01:1));
cr_interp(100)=cr_interp(99);
layer=xlsread('meyer_layers_percentage.csv');
layers=round(100.*(1-layer(:,2)));
bar_cr=[mean(cr_interp(1:layers(1)));mean(cr_interp(layers(1)+1:layers(3)));mean(cr_interp(layers(3)+1:layers(4)));mean(cr_interp(layers(4)+1:layers(6)));mean(cr_interp(layers(6)+1:100))]./100;
bar_som=[mean(som_interp(1:layers(1)));mean(som_interp(layers(1)+1:layers(3)));mean(som_interp(layers(3)+1:layers(4)));mean(som_interp(layers(4)+1:layers(6)));mean(som_interp(layers(6)+1:100))]./100;
bar_pv=[mean(pv_interp(1:layers(1)));mean(pv_interp(layers(1)+1:layers(3)));mean(pv_interp(layers(3)+1:layers(4)));mean(pv_interp(layers(4)+1:layers(6)));mean(pv_interp(layers(6)+1:100))]./100;

