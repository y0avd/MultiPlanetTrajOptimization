function plotLambertArc(r1,r2,TOF,mu)
    nump = 1e3;

    [v1,v2] = lambert(r1,r2,TOF,0,mu);

    % velocity vector scaling
    K = 1e7;

    plot3(r1(1),r1(2),r1(3),'go')
    plot3([r1(1),r1(1) + K*v1(1)],...
        [r1(2),r1(2) + K*v1(2)],...
        [r1(3),r1(3) + K*v1(3)],'g')

    plot3(r2(1),r2(2),r2(3),'ro')
    plot3([r2(1),r2(1) + K*v2(1)],...
        [r2(2),r2(2) + K*v2(2)],...
        [r2(3),r2(3) + K*v2(3)],'r')
    
%     [a,eccx,eccy,inc,raan] = cart2orb(r1(1),r1(2),r1(3),v1(1),v1(2),v1(3),mu);
% 
%     th = linspace(0,360,nump);
% 
%     for i = 1:length(th)
%         [x(i),y(i),z(i)] = orb2cart(a,eccx,eccy,inc,raan,th(i),mu); %#ok<AGROW> 
%     end
% 
%     plot3(x,y,z,'b')
end

