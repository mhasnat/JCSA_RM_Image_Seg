function imcshow(img, op)
    if(nargin<2)
        op=1;
    end
    
    if(op==1)
        figure, imshow(label2rgb(img));
    elseif(op==2)
        imshow(label2rgb(img));
    elseif(op==3)
        img = assignRandomLabel(img);
        figure, imshow(label2rgb(img));
    elseif(op==4)
        img = assignRandomLabel(img);
        imshow(label2rgb(img));
    end
end