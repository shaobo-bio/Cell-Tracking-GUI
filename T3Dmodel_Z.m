function T3Dmodel_Z(Res)

figure;
[W L]=size(Res);
L=L/3;
[sx sy sz]=sphere(30);


for i=1:L
    x=Res(:,3*i-2);
    y=Res(:,3*i-1);
    z=(10-Res(:,3*i))*4;
    
    plot3(Res(:,1),Res(:,2),Res(:,3),'w*');
    hold on;
    plot3(Res(:,W-2),Res(:,W-1),(10-Res(:,W))*4,'w*');
    box on;
    axis equal ;
    view([-20,10]);
    
    for j=1:W
        
        surf(5*sx+x(j),5*sy+y(j),5*sz+z(j));
        alpha(0.8)         %?????
        shading flat
        
    end
    hold off;
    
    
    pause(1);
end