function centers=FindcenterBW(I,thred)
I_eq = adapthisteq(I);
bw = im2bw(I_eq, graythresh(I_eq));
bw(find(I < thred))=0;
bw2 = imfill(bw,'holes');
bw3 = imopen(bw2, ones(2,2));
bw4 = bwareaopen(bw3, 40);
bw4_perim = bwperim(bw4);

mask_em = imextendedmax(I_eq, 30);

mask_em = imclose(mask_em, ones(5,5));
mask_em = imfill(mask_em, 'holes');
mask_em = bwareaopen(mask_em, 40);

I_eq_c = imcomplement(I_eq);
I_mod = imimposemin(I_eq_c, ~bw4 | mask_em);
L = watershed(I_mod);
L(find(L>0))=1;
Centers=regionprops(L,'Centroid');
centers=cat(1,Centers.Centroid);
end