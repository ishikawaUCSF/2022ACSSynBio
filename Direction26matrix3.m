function M = Direction26matrix3

M = zeros(26);
M(:,1) = [2;1;1;1;1;1;1;1;1;0;0;0;0;0;0;0;0;-4;-4;-4;-4;-4;-4;-4;-4;-4];
M(:,2) = [0;1;1;-1;-3;-3;-3;-1;1;1;1;0;-2;-2;-2;0;1;0;0;-2;-4;-4;-4;-2;0;-2];
M(:,3) = [0;1;1;1;-1;-3;-3;-3;-1;1;1;1;0;-2;-2;-2;0;0;0;0;-2;-4;-4;-4;-2;-2];
M(:,4) = [0;-1;1;1;1;-1;-3;-3;-3;0;1;1;1;0;-2;-2;-2;-2;0;0;0;-2;-4;-4;-4;-2];
M(:,5) = [0;-3;-1;1;1;1;-1;-3;-3;-2;0;1;1;1;0;-2;-2;-4;-2;0;0;0;-2;-4;-4;-2];
M(:,6) = [0;-3;-3;-1;1;1;1;-1;-3;-2;-2;0;1;1;1;0;-2;-4;-4;-2;0;0;0;-2;-4;-2];
M(:,7) = [0;-3;-3;-3;-1;1;1;1;-1;-2;-2;-2;0;1;1;1;0;-4;-4;-4;-2;0;0;0;-2;-2];
M(:,8) = [0;-1;-3;-3;-3;-1;1;1;1;0;-2;-2;-2;0;1;1;1;-2;-4;-4;-4;-2;0;0;0;-2];
M(:,9) = [0;1;-1;-3;-3;-3;-1;1;1;1;0;-2;-2;-2;0;1;1;0;-2;-4;-4;-4;-2;0;0;-2];
M(:,10) = [-1;1;1;-1;-3;-3;-3;-1;1;2;2;0;-2;-2;-2;0;2;1;1;-1;-3;-3;-3;-1;1;-1];
M(:,11) = [-1;1;1;1;-1;-3;-3;-3;-1;2;2;2;0;-2;-2;-2;0;1;1;1;-1;-3;-3;-3;-1;-1];
M(:,12) = [-1;-1;1;1;1;-1;-3;-3;-3;0;2;2;2;0;-2;-2;-2;-1;1;1;1;-1;-3;-3;-3;-1];
M(:,13) = [-1;-3;-1;1;1;1;-1;-3;-3;-2;0;2;2;2;0;-2;-2;-3;-1;1;1;1;-1;-3;-3;-1];
M(:,14) = [-1;-3;-3;-1;1;1;1;-1;-3;-2;-2;0;2;2;2;0;-2;-3;-3;-1;1;1;1;-1;-3;-1];
M(:,15) = [-1;-3;-3;-3;-1;1;1;1;-1;-2;-2;-2;0;2;2;2;0;-3;-3;-3;-1;1;1;1;-1;-1];
M(:,16) = [-1;-1;-3;-3;-3;-1;1;1;1;0;-2;-2;-2;0;2;2;2;-1;-3;-3;-3;-1;1;1;1;-1];
M(:,17) = [-1;1;-1;-3;-3;-3;-1;1;1;2;0;-2;-2;-2;0;2;2;1;-1;-3;-3;-3;-1;1;1;-1];
M(:,18) = [-2;0;0;-2;-4;-4;-4;-2;0;1;1;0;-2;-2;-2;0;1;1;1;-1;-3;-3;-3;-1;1;0];
M(:,19) = [-2;0;0;0;-2;-4;-4;-4;-2;1;1;1;0;-2;-2;-2;0;1;1;1;-1;-3;-3;-3;-1;0];
M(:,20) = [-2;-2;0;0;0;-2;-4;-4;-4;0;1;1;1;0;-2;-2;-2;-1;1;1;1;-1;-3;-3;-3;0];
M(:,21) = [-2;-4;-2;0;0;0;-2;-4;-4;-2;0;1;1;1;0;-2;-2;-3;-1;1;1;1;-1;-3;-3;0];
M(:,22) = [-2;-4;-4;-2;0;0;0;-2;-4;-2;-2;0;1;1;1;0;-2;-3;-3;-1;1;1;1;-1;-3;0];
M(:,23) = [-2;-4;-4;-4;-2;0;0;0;-2;-2;-2;-2;0;1;1;1;0;-3;-3;-3;-1;1;1;1;-1;0];
M(:,24) = [-2;-2;-4;-4;-4;-2;0;0;0;0;-2;-2;-2;0;1;1;1;-1;-3;-3;-3;-1;1;1;1;0];
M(:,25) = [-2;0;-2;-4;-4;-4;-2;0;0;1;0;-2;-2;-2;0;1;1;1;-1;-3;-3;-3;-1;1;1;0];
M(:,26) = [-4;-4;-4;-4;-4;-4;-4;-4;-4;0;0;0;0;0;0;0;0;1;1;1;1;1;1;1;1;2];