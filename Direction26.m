function M = Direction26(radius)

% This function make 26 directional shapes. (170811)
% M includes 26 directional shapes and pixel number.
% Need functions 'rotation_mat3' and 'inhull'.


[xgrid, ygrid, zgrid] = meshgrid(-radius:radius);
ball = (sqrt(xgrid.^2 + ygrid.^2 + zgrid.^2) <= radius);

%figure, sphere(8);
[X1,Y1,Z1] = sphere(8);

RM = rotation_mat3([0;0;1], pi/8); %make a rotation matrix 22.5º(pi/8) around z axis
X2 = reshape(X1, [], 1);
Y2 = reshape(Y1, [], 1);
Z2 = reshape(Z1, [], 1);
XYZ1 = [X2, Y2, Z2]*RM; %rotate shpere
X3 = reshape(XYZ1(:,1), 9, 9);
Y3 = reshape(XYZ1(:,2), 9, 9);
Z3 = reshape(XYZ1(:,3), 9, 9);
%figure, surf(X3, Y3, Z3);

X4 = (X3 * radius) + radius + 1;
Y4 = (Y3 * radius) + radius + 1;
Z4 = (Z3 * radius) + radius + 1;
Center1 = [radius + 1, radius + 1, radius + 1];

direction = cell(26,1);
M = cell(26,4);

direction(1,1) = {vertcat(([reshape(X4(1:2,:),[],1), reshape(Y4(1:2,:),[],1), reshape(Z4(1:2,:),[],1)]), Center1)};
direction(2,1) = {vertcat(([reshape(X4(2:4,1:2),[],1), reshape(Y4(2:4,1:2),[],1), reshape(Z4(2:4,1:2),[],1)]), Center1)};
direction(3,1) = {vertcat(([reshape(X4(2:4,2:3),[],1), reshape(Y4(2:4,2:3),[],1), reshape(Z4(2:4,2:3),[],1)]), Center1)};
direction(4,1) = {vertcat(([reshape(X4(2:4,3:4),[],1), reshape(Y4(2:4,3:4),[],1), reshape(Z4(2:4,3:4),[],1)]), Center1)};
direction(5,1) = {vertcat(([reshape(X4(2:4,4:5),[],1), reshape(Y4(2:4,4:5),[],1), reshape(Z4(2:4,4:5),[],1)]), Center1)};
direction(6,1) = {vertcat(([reshape(X4(2:4,5:6),[],1), reshape(Y4(2:4,5:6),[],1), reshape(Z4(2:4,5:6),[],1)]), Center1)};
direction(7,1) = {vertcat(([reshape(X4(2:4,6:7),[],1), reshape(Y4(2:4,6:7),[],1), reshape(Z4(2:4,6:7),[],1)]), Center1)};
direction(8,1) = {vertcat(([reshape(X4(2:4,7:8),[],1), reshape(Y4(2:4,7:8),[],1), reshape(Z4(2:4,7:8),[],1)]), Center1)};
direction(9,1) = {vertcat(([reshape(X4(2:4,8:9),[],1), reshape(Y4(2:4,8:9),[],1), reshape(Z4(2:4,8:9),[],1)]), Center1)};
direction(10,1) = {vertcat(([reshape(X4(4:6,1:2),[],1), reshape(Y4(4:6,1:2),[],1), reshape(Z4(4:6,1:2),[],1)]), Center1)};
direction(11,1) = {vertcat(([reshape(X4(4:6,2:3),[],1), reshape(Y4(4:6,2:3),[],1), reshape(Z4(4:6,2:3),[],1)]), Center1)};
direction(12,1) = {vertcat(([reshape(X4(4:6,3:4),[],1), reshape(Y4(4:6,3:4),[],1), reshape(Z4(4:6,3:4),[],1)]), Center1)};
direction(13,1) = {vertcat(([reshape(X4(4:6,4:5),[],1), reshape(Y4(4:6,4:5),[],1), reshape(Z4(4:6,4:5),[],1)]), Center1)};
direction(14,1) = {vertcat(([reshape(X4(4:6,5:6),[],1), reshape(Y4(4:6,5:6),[],1), reshape(Z4(4:6,5:6),[],1)]), Center1)};
direction(15,1) = {vertcat(([reshape(X4(4:6,6:7),[],1), reshape(Y4(4:6,6:7),[],1), reshape(Z4(4:6,6:7),[],1)]), Center1)};
direction(16,1) = {vertcat(([reshape(X4(4:6,7:8),[],1), reshape(Y4(4:6,7:8),[],1), reshape(Z4(4:6,7:8),[],1)]), Center1)};
direction(17,1) = {vertcat(([reshape(X4(4:6,8:9),[],1), reshape(Y4(4:6,8:9),[],1), reshape(Z4(4:6,8:9),[],1)]), Center1)};
direction(18,1) = {vertcat(([reshape(X4(6:8,1:2),[],1), reshape(Y4(6:8,1:2),[],1), reshape(Z4(6:8,1:2),[],1)]), Center1)};
direction(19,1) = {vertcat(([reshape(X4(6:8,2:3),[],1), reshape(Y4(6:8,2:3),[],1), reshape(Z4(6:8,2:3),[],1)]), Center1)};
direction(20,1) = {vertcat(([reshape(X4(6:8,3:4),[],1), reshape(Y4(6:8,3:4),[],1), reshape(Z4(6:8,3:4),[],1)]), Center1)};
direction(21,1) = {vertcat(([reshape(X4(6:8,4:5),[],1), reshape(Y4(6:8,4:5),[],1), reshape(Z4(6:8,4:5),[],1)]), Center1)};
direction(22,1) = {vertcat(([reshape(X4(6:8,5:6),[],1), reshape(Y4(6:8,5:6),[],1), reshape(Z4(6:8,5:6),[],1)]), Center1)};
direction(23,1) = {vertcat(([reshape(X4(6:8,6:7),[],1), reshape(Y4(6:8,6:7),[],1), reshape(Z4(6:8,6:7),[],1)]), Center1)};
direction(24,1) = {vertcat(([reshape(X4(6:8,7:8),[],1), reshape(Y4(6:8,7:8),[],1), reshape(Z4(6:8,7:8),[],1)]), Center1)};
direction(25,1) = {vertcat(([reshape(X4(6:8,8:9),[],1), reshape(Y4(6:8,8:9),[],1), reshape(Z4(6:8,8:9),[],1)]), Center1)};
direction(26,1) = {vertcat(([reshape(X4(8:9,:),[],1), reshape(Y4(8:9,:),[],1), reshape(Z4(8:9,:),[],1)]), Center1)};

M{1,2} = [0 0 -1];
M{2,2} = [-1 0 -1];
M{3,2} = [-1 -1 -1];
M{4,2} = [0 -1 -1];
M{5,2} = [1 -1 -1];
M{6,2} = [1 0 -1];
M{7,2} = [1 1 -1];
M{8,2} = [0 1 -1];
M{9,2} = [-1 1 -1];
M{10,2} = [-1 0 0];
M{11,2} = [-1 -1 0];
M{12,2} = [0 -1 0];
M{13,2} = [1 -1 0];
M{14,2} = [1 0 0];
M{15,2} = [1 1 0];
M{16,2} = [0 1 0];
M{17,2} = [-1 1 0];
M{18,2} = [-1 0 1];
M{19,2} = [-1 -1 1];
M{20,2} = [0 -1 1];
M{21,2} = [1 -1 1];
M{22,2} = [1 0 1];
M{23,2} = [1 1 1];
M{24,2} = [0 1 1];
M{25,2} = [-1 1 1];
M{26,2} = [0 0 1];

M{1,4} = 26;
M{2,4} = 22;
M{3,4} = 23;
M{4,4} = 24;
M{5,4} = 25;
M{6,4} = 18;
M{7,4} = 19;
M{8,4} = 20;
M{9,4} = 21;
M{10,4} = 14;
M{11,4} = 15;
M{12,4} = 16;
M{13,4} = 17;
M{14,4} = 10;
M{15,4} = 11;
M{16,4} = 12;
M{17,4} = 13;
M{18,4} = 6;
M{19,4} = 7;
M{20,4} = 8;
M{21,4} = 9;
M{22,4} = 2;
M{23,4} = 3;
M{24,4} = 4;
M{25,4} = 5;
M{26,4} = 1;

sizeBall = size(ball);


for n1 = 1:26;
    MAT0 = true(sizeBall);
    PxListBall = find(MAT0);
    [A1,B1,C1] = ind2sub(sizeBall, PxListBall);
    SUB = [A1,B1,C1];
    K = convhulln(direction{n1,1});
    %figure, trisurf(K02,direction{n1,1},direction02(:,2),direction02(:,3));
    in = inhull(SUB,direction{n1,1},K, 0.25);
    A2 = SUB(in,1);
    B2 = SUB(in,2);
    C2 = SUB(in,3);
    MAT1 = false(sizeBall);
    for n2 = 1:length(A2);
        MAT1(A2(n2), B2(n2), C2(n2)) = 1;
    end
    MAT2 = ball & MAT1;
    MAT2(radius+1, radius+1, radius+1) = 0;
    M{n1,1} = MAT2;
    M{n1,3} = sum(MAT2(:));
end


    