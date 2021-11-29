function [I,centers]=FindcenterSobel(im,thred1,thred2)

hy = fspecial('sobel');
hx = hy';
Iy = imfilter(double(im), hy, 'replicate');
Ix = imfilter(double(im), hx, 'replicate');
I = sqrt(Ix.^2 + Iy.^2);
thred1=thred1+1000;
I(find(I<thred1))=0;
BW=zeros(size(I));
BW(find(I>thred1))=1;
Fill=imfill(BW,'holes');
Fill(find(I>thred1))=0;

[centers,~,~]=imfindcircles(BW,[2 4],'Method','TwoStage','Sensitivity',thred2);

