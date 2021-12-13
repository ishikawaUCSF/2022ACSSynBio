function L = Neighbor26Matrix1(index)

% return length from the center point (2,2,2)
% input should be 1~27
% 

sizeIndex = size(index);
L = zeros(sizeIndex);
A = 14;
B = [5, 11, 13, 15, 17, 23];
C = [2, 4, 6, 8, 10, 12, 16, 18, 20, 22, 24, 26];
D = [1, 3, 7, 9, 19, 21, 25, 27];

for n01 = 1:length(index);
    if find(A==index(n01)) > 0
        L(n01) = 0;
    elseif find(B==index(n01)) > 0
        L(n01) = 1;
    elseif find(C==index(n01)) > 0
        L(n01) = 1.414213562373095;
    elseif find(D==index(n01)) > 0
        L(n01) = 1.732050807568877;
    else
        L(n01) = NaN;
    end
end



