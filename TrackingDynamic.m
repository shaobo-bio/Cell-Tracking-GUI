
handles=Pro;
YD=handles.YD;
Res=handles.Res;
Img=handles.Img;
D=handles.D;
DP=handles.DP;
Way=handles.TrackingWay;
figure;
for Index=21:-1:1
    InJ=handles.ImgIn+Way*(Index-1);
%     in=handles.in;
    
    
    imshow(Img(:,:,InJ),[]);
    hold on;
    
    P=handles.P(Index,:);
    x1=linspace(0,500);
    y1=polyval(P,x1);
    %plot(y,x,'--w');
    [L,W]=size(Res);
    PT(L+1)=plot(y1,x1,'--m');
    
    Count=1;
    for i=1:L
        Col=['r','g','b','y'];
        
        A=Res(i,1:Index*2);
        A=reshape(A,2,Index);
        A=A';
        if i > 2 && mod(i,2)
            if ~isequal(Res(i,1:2),Res(i-2,1:2))
                if ~isequal(Res(i+1,1:2),Res(i-1,1:2))
                    Count=Count+1;
                end
            end
        end
        
        M=mod(Count,4);
        PT(i)=plot(A(:,1),A(:,2),'-g','LineWidth',1);
        
        plot(A(Index,1),A(Index,2),strcat(Col(M+1),'*'),'MarkerSize',3);
        
        
    end
    hold off;
    MM(22-Index)=getframe;
end