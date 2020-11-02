function Out = synthesize(In, outRows, outCols, windowSize)

%% INITIALIZE CONSTANTS
if mod(windowSize,2) == 0
    windowSize = windowSize + 1;
end
h = (windowSize - 1) / 2;
center = (windowSize^2 + 1) / 2;
Gauss = fspecial('gaussian', windowSize, windowSize/6);

%% GET INPUT IMAGE NEIGHBORHOODS
[inRows, inCols, channels] = size(In);
Nhoods = zeros(windowSize^2, (inRows-2*h)*(inCols-2*h), channels);
for ch = 1:channels
    Nhoods(:,:,ch) = im2col(In(:,:,ch), [windowSize windowSize]);
end

%% INITIALIZE OUTPUT IMAGE
Out = -ones(outRows, outCols, channels);
Out(1:inRows, 1:inCols, :) = In;
PadOut = -ones(outRows+2*h, outCols+2*h, channels);
PadOut(h+1:h+outRows, h+1:h+outCols, :) = Out;
Fill = Out(:,:,1) >= 0;
PadKnowns = zeros(outRows+2*h, outCols+2*h);
for i = 1:inRows
    for j = 1:inCols
        PadKnowns(i:i+2*h, j:j+2*h) = PadKnowns(i:i+2*h, j:j+2*h) + 1;
    end
end

%% ALGORITHM
for iteration = 1:(outRows*outCols - inRows*inCols)
    %% GET UNFILLED PIXEL WITH MOST KNOWN NEIGHBORS
    Knowns = PadKnowns(h+1:h+outRows, h+1:h+outCols);
    [~,ind] = max(Knowns.*~Fill, [], 'all', 'linear');
    [i,j] = ind2sub([outRows outCols], ind);
    
    %% GET WINDOW
    Window = zeros(windowSize^2, channels);
    for ch = 1:channels
        Window(:,ch) = reshape(PadOut(i:i+2*h, j:j+2*h, ch), [], 1);
    end
    valid = Window(:,1) >= 0;
    
    %% EVALUATE SIMILARITIES
    SD = zeros(nnz(valid), size(Nhoods,2));
    for ch = 1:channels
        SD = SD + (Nhoods(valid,:,ch) - Window(valid,ch)) .^ 2;
    end
    GWSD = SD .* Gauss(valid);
    
    %% PICK BEST MATCH
    [~,best] = min(sum(GWSD));
    Out(i,j,:) = Nhoods(center,best,:);
    PadOut(h+i, h+j, :) = Nhoods(center,best,:);
    Fill(i,j) = true;
    PadKnowns(i:i+2*h, j:j+2*h) = PadKnowns(i:i+2*h, j:j+2*h) + 1;
    imshow(Out)
end