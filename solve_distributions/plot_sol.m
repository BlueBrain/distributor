%plot heat map of solution
function [] = plot_sol(mat1)

x=0.*mat1;
y=0.*mat1;
for count1=1:100
    for count2=1:40
        x(count1,count2)=count1;
        y(count1,count2)=count2;
    end
end

% figure
surf(x,y,mat1(:,:))
% ylabel('Cell Area ( x10 pixels)')
% xlabel('Fraction Cortical Depth')
az = 0;
el = 90;
view(az, el);
set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',[screenposition(3:4)]);

set(gca,'visible','off')
end
