%get thickness of region delineated by points
function [distout]=get_layer_thickness(points)
     p1=(points(1,:));
    p2=(points(2,:));
    p3=(points(3,:));
    p4=(points(4,:));

    cp1=(p1+p2)./2;
    cp2=(p3+p4)./2;

    v1=cp2-cp1;
    v1=[v1,0];
    v2=cross(v1,[0,0,1]);
    v2=v2./sqrt(sum(v2.*v2));

    top1=(p1+p4)./2;
    bottom1=(p2+p3)./2;

    allcenter1=(p1+p2+p3+p4)./4;

    vecTop=top1-allcenter1;
    vecBottom=bottom1-allcenter1;

    dist1=sum(vecTop .*v2(1:2));
    dist2=sum(vecBottom .*v2(1:2));
    distout=dist1-dist2;
if (distout<500)
  distout=1299;
end
