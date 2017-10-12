function maskedImg = getMaskedImage(I, mask)
[r c d] = size(I);

maskedImg = [];
for i=1:d
    tI = I(:,:,i);
    maskedImg(:,i) = tI(mask);
end

r1 = 427;
c1 = 561;
maskedImg = reshape(maskedImg, r1, c1, d);

if(strcmp(class(I), 'uint8'))
    maskedImg = uint8(maskedImg);
end

end