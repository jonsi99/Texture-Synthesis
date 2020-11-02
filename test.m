In = im2double(imread('texture/texture1.jpg'));
Out = synthesize(In,128,128,17);
imwrite(Out,'result1_17.jpg')