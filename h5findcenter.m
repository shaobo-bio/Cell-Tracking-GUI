function centers=h5findcenter(img,thred)

I=zeros(size(img));
I(img>thred)=img(img>thred);


g=medfilt2(I,[2 2],'symmetric');
g1=im2uint8(g);
g1(find(g1<160))=0;
g1(find(g1>=160))=1;
imshow(g1);
bw=bwmorph(g1,'clean');
bw=bwmorph(bw,'hbreak',Inf);

 % labeled image

centers=FindCenters2(bw,50);
centers=round(centers);


% figure;
% L=bwlabel(g1);
% coloredLabels = label2rgb (L, 'hsv', 'k', 'shuffle');
% imagesc(coloredLabels);
% hold on;
% plot(Centers(:,1),Centers(:,2),'w*');

end



% [c1,~,~]=imfindcircles(L,[1 4],'Method','TwoStage','Sensitivity',0.9);
% plot(c1(:,1),c1(:,2),'r*');
% 
% A=regionprops(bw,'Area');
% A=cat(1,A.Area);
% unique(A);

% bw_perim = bwperim(bw);
% overlay1 = imoverlay(bw, bw_perim, [.3 1 .3]);
% 
% imshow(overlay1);

function Centers=FindCenters(bw,Thred)

L=bwlabel(bw);

A=regionprops(L,'Area');
A=cat(1,A.Area);

Check=max(A)>Thred;

Centers=regionprops(L,'Centroid');
Centers=cat(1,Centers.Centroid);

if Check==1 && Thred > 10
    
    Mark=find(A>Thred);
    MarkPart=ismember(L, Mark);
    bw_perim = bwperim(MarkPart);
    Fill=imfill(bw_perim,'holes');
    Fill(find(bw_perim==1))=0;
    Fill=bwmorph(Fill,'hbreak');

    Thred=Thred-16;
    CAdd=FindCenters(Fill,Thred);
%     CAdd=regionprops(Fill,'Centroid');
%     CAdd=cat(1,CAdd.Centroid);
    
    
    Centers(Mark,:)=[];
    Centers=[Centers;CAdd];
    
end

end

function Centers=FindCenters2(bw,Thred)

L=bwlabel(bw);

A=regionprops(L,'Area');
A=cat(1,A.Area);

Centers=regionprops(L,'Centroid');
Centers=cat(1,Centers.Centroid);


Mark=find(A>Thred);
MarkPart=ismember(L, Mark);

[c1,~,~]=imfindcircles(MarkPart,[1 3],'Method','TwoStage','Sensitivity',1);

Centers(Mark,:)=[];
Centers=[Centers;c1];
end