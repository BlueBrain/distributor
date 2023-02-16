[meyer_cr,meyer_som,meyer_pv]=get_meyer_percentages();

basedir='C:\Users\dakeller\Dropbox\Cell Densities Working Folder\Analyses\Images\';
%define constants
L0=0;
L1=10;%9
L2=19;
L3=33;
L4=43;%46
L5a=53;
L5=69;%65
L6=100;%65

map2=[0,0,0;217,83,25;237,177,32]./255;
%colormap(map2)

scale0=2e6;
interp_me_x=[L1/2,(L1+L3)./2,(L3+L4)./2,(L4+L5)./2,(L5+L6)./2];
my_inh2=scale0.*[mean(my_inh(1:L1)),mean(my_inh(L1+1:L3)),mean(my_inh(L3+1:L4)),mean(my_inh(L4+1:L5)),mean(my_inh(L5+1:L6))];
my_cr2=scale0.*[mean(my_cr(1:L1)),mean(my_cr(L1+1:L3)),mean(my_cr(L3+1:L4)),mean(my_cr(L4+1:L5)),mean(my_cr(L5+1:L6))]./my_inh2;
my_pv2=scale0.*[mean(my_pv(1:L1)),mean(my_pv(L1+1:L3)),mean(my_pv(L3+1:L4)),mean(my_pv(L4+1:L5)),mean(my_pv(L5+1:L6))]./my_inh2;
my_sst2=scale0.*[mean(my_sst(1:L1)),mean(my_sst(L1+1:L3)),mean(my_sst(L3+1:L4)),mean(my_sst(L4+1:L5)),mean(my_sst(L5+1:L6))]./my_inh2;

gonchar_sst_bottom=dlmread(strcat(basedir,'gonchar_sst_bottom.csv'));
gonchar_cr_top=dlmread(strcat(basedir,'gonchar_cr_top.csv'));
gonchar_pv_top=dlmread(strcat(basedir,'gonchar_pv_top.csv'));
gonchar_sst_top=dlmread(strcat(basedir,'gonchar_sst_top.csv'));

gonchar_sst_bottom=gonchar_sst_bottom(:,2);
gonchar_cr_top=gonchar_cr_top(:,2);
gonchar_pv_top=gonchar_pv_top(:,2);
gonchar_sst_top=gonchar_sst_top(:,2);

gonchar_pv=gonchar_pv_top;
gonchar_cr=gonchar_cr_top-gonchar_pv_top;
gonchar_sst=gonchar_sst_top-gonchar_sst_bottom;


lee_sst_high=dlmread(strcat(basedir,'lee_sst_high.csv'));        
lee_sst_low=dlmread(strcat(basedir,'lee_sst_low.csv'));
lee_5htr3_low=dlmread(strcat(basedir,'lee_5htr3_low.csv'));
lee_5htrs_high=dlmread(strcat(basedir,'lee_5htrs_high.csv'));
lee_pv_high=dlmread(strcat(basedir,'lee_pv_high.csv'));

lee_pv_high=lee_pv_high(:,1);
lee_5htrs_high=lee_5htrs_high(:,1);
lee_5htr3_low=lee_5htr3_low(:,1);
lee_sst_low=lee_sst_low(:,1);
lee_sst_high=lee_sst_high(:,1);

lee_pv=lee_pv_high;
lee_sst=lee_sst_high-lee_sst_low;
lee_5htr3=lee_5htrs_high-lee_5htr3_low;

meyer_inh=dlmread(strcat(basedir,'meyer_inh.csv'));

interp_lee_sst_x=[0,L1/2,(L1+L3)./2,(L3+L4)./2,(L4+L5)./2,(L5+L6)./2,L6];
interp_lee_sst_y=[lee_sst(1),lee_sst(1),lee_sst(2),lee_sst(3),lee_sst(4),lee_sst(5),lee_sst(5)]./100;
interp_gonchar_sst_x=[0,L1/2,(L1+L3)./2,(L3+L4)./2,(L4+L5)./2,(L5+L6)./2,L6];
interp_gonchar_sst_y=[gonchar_sst(1),gonchar_sst(1),gonchar_sst(2),gonchar_sst(3),gonchar_sst(4),gonchar_sst(5),gonchar_sst(5)]./100;

lee_sst2=(interp_lee_sst_y)   .* interp1( meyer_inh(:,2),meyer_inh(:,1),interp_lee_sst_x);
gonchar_sst2=interp1(interp_gonchar_sst_x,interp_gonchar_sst_y,interp_gonchar_sst_x)   .* interp1( meyer_inh(:,2),meyer_inh(:,1),interp_gonchar_sst_x );

figure
X = categorical({'2/3','4','5','6'});
%m_sst2=[transpose(lee_sst2(2:6)),transpose(gonchar_sst2(2:6)),transpose(my_sst2)]
%m_sst2=[transpose(my_sst2),(gonchar_sst./100),(lee_sst./100),meyer_som];
m_sst2=[transpose(my_sst2),meyer_som,(lee_sst./100)];
m_sst2=[transpose(my_sst2),meyer_som];
m_sst2(1,:)=[];
b=bar(X,m_sst2);
for k = 1:2
    b(k).FaceColor = map2(k,:);
end
title('SST')

interp_lee_pv_x=[0,L1/2,(L1+L3)./2,(L3+L4)./2,(L4+L5)./2,(L5+L6)./2,L6];
interp_lee_pv_y=[lee_pv(1),lee_pv(1),lee_pv(2),lee_pv(3),lee_pv(4),lee_pv(5),lee_pv(5)]./100;
lee_pv2=interp_lee_pv_y   .* interp1( meyer_inh(:,2),meyer_inh(:,1),interp_lee_pv_x) ...

interp_gonchar_pv_x=[0,L1/2,(L1+L3)./2,(L3+L4)./2,(L4+L5)./2,(L5+L6)./2,L6];
interp_gonchar_pv_y=[gonchar_pv(1),gonchar_pv(1),gonchar_pv(2),gonchar_pv(3),gonchar_pv(4),gonchar_pv(5),gonchar_pv(5)]./100;
gonchar_pv2=interp_gonchar_pv_y   .* interp1( meyer_inh(:,2),meyer_inh(:,1),interp_gonchar_pv_x) ;

figure
%m_pv2=[transpose(lee_pv2(2:6)),transpose(gonchar_pv2(2:6)),transpose(my_pv2)]
%m_pv2=[transpose(my_pv2),(gonchar_pv./100),(lee_pv./100)];
%m_pv2=[transpose(my_pv2),meyer_pv,(lee_pv./100)];
m_pv2=[transpose(my_pv2),meyer_pv]
m_pv2(1,:)=[];
b=bar(X,m_pv2);
for k = 1:2
    b(k).FaceColor = map2(k,:);
end
ylim([0 0.8]);
title('PV')

gonchar_cr=gonchar_cr_top-gonchar_pv_top;
interp_gonchar_cr_x=[0,L1/2,(L1+L3)./2,(L3+L4)./2,(L4+L5)./2,(L5+L6)./2,L6];
interp_gonchar_cr_y=[gonchar_cr(1),gonchar_cr(1),gonchar_cr(2),gonchar_cr(3),gonchar_cr(4),gonchar_cr(5),gonchar_cr(5)]./100;
gonchar_cr2=interp_gonchar_cr_y  .* interp1( meyer_inh(:,2),meyer_inh(:,1),interp_gonchar_cr_x);

figure
%m_cr2=[transpose(gonchar_cr2(2:6)),transpose(my_cr2)]
%m_cr2=[transpose(my_cr2),(gonchar_cr./100),meyer_cr];
m_cr2=[transpose(my_cr2),meyer_cr];
m_cr2(1,:)=[];
b=bar(X,m_cr2);
%cmap = colormap(jet);
for k = 1:2
    b(k).FaceColor = map2(k,:);
end
title('CR')