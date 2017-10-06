function Igf = getRGBDGradient(rgbImg, depImg)

[Igmag(:,:,1), ~] = imgradient(rgbImg(:,:,1));
[Igmag(:,:,2), ~] = imgradient(rgbImg(:,:,2));
[Igmag(:,:,3), ~] = imgradient(rgbImg(:,:,3));

if(nargin>1)
    [Igmag(:,:,4), ~] = imgradient(depImg);
end

Igf = max(Igmag, [], 3);