%% UNIVERSITY OF THESSALY,DEPARTMENT OF ELECTRICAL AND COMPUTER ENGINEERING
%%
%%            CE434 : ELECTROMAGNETIC FIELDS II
%%
%% COMPUTING PROJECT (Fall 2014)
%% Solution of Laplace equation using finite differences method 
%% in 2-dimensional areas for finding electric potential V = V(x, y) 
%%
%% Nonas Evaggelos (AEM 1321)
%% nonas@inf.uth.gr
%%


%% EXCERCISE 1

clear all; clc;

M  = 80;  % number of points
Vu = 100; % potential of the upper plate

a = 0; b = 4;
step  = (b - a)/M;
start = a + step;
stop  = b - step;
[x, y] = meshgrid(start : step : stop);

disp('Excercise 1 : Solution in rectangular area');
disp('------------------------------------------');
disp('Compare to theoritical solution : ');

% potential based on theoretical solution
Z = zeros(M-1,M-1);
for n = 1 : 2 : 19
    c = (4*Vu) / (n*pi*sinh(n*pi));
    Z = Z + c * ( sin((n*pi/b).*x) .* sinh((n*pi/b).*y) );
end

% finite differences method
N = (M-1)^2;
L  = laplacian(N, M-1);

% create RHS vector b
r = zeros(N, 1);
for i = (N-M+1) : 1 : N
    r(i) = Vu;
end
v = L\r;
V = -vec2mat(v, M-1);

% plot the results
figure(1);
surf(x, y, Z);
colorbar
title('Laplace Equation : Theoritical Solution for 10 terms');
xlabel('x axis');
ylabel('y axis');
zlabel('V(x,y)');

figure(2);
surf(x, y, V);
colorbar
title('Laplace Equation : Finite Differences Method Solution');
xlabel('x axis');
ylabel('y axis');
zlabel('V(x,y)');

% print error informations
fprintf('(1-Norm) err = %i\n', norm(Z-V, 1)/norm(Z,1)*100);
fprintf('(2-Norm) err = %i\n', norm(Z-V, 2)/norm(Z,2)*100);
fprintf('(Inf-Norm) err = %i\n', norm(Z-V, Inf)/norm(Z,Inf)*100);

disp('');
disp('Press any key to continue...');
pause;


%% EXCERCISE 2

clear all; clc;

disp('Excercise 2');
disp('-----------');

M  = 80;  % number of points
Vu = 100; % potential of the right plates

a = 0; c = 1; b = 4;
step  = (b - a)/M;
start = a + step;
stop  = b - step;
[x, y] = meshgrid(start : step : stop);

M1 = M-1;
M2 = floor((c-a)/step)-1;
M3 = M2;
M4 = M1-M2;

[L, r, N] = extract_matrix2(Vu, M1, M2, M3, M4);

v = L\r;

V = zeros(M1, M1);
for i = 1 : 1 : M2
    for j = 1 : 1 : M1
        V(i, j) = v((i-1)*M1+j);
    end
end

for i = 1 : 1 : M4
    for j = 1 : 1 : M3
        V(M2 + i, j) = v(M2*M1 + (i-1)*M3 + j); 
    end
end

figure(3);
surf(x, y, V);
colorbar;
title('Laplace Equation : Finite Differences Method Solution');
xlabel('x axis');
ylabel('y axis');
zlabel('V(x,y)');

disp('');
disp('Press any key to continue...');
pause;


%% EXCERCISE 3

clear all; clc;

disp('Excercise 3');
disp('-----------');

M  = 80;  % number of points
Vu = 100; % potential of the surround plates

a = 0; c = 1; b = 4;
step  = (b - a)/M;
start = a + step;
stop  = b - step;
[x, y] = meshgrid(start : step : stop);

M1 = M-1;
M2 = floor((c-a)/step)-1;

M3 = M2;
M4 = M1-2*M2;

[L, r, N] = extract_matrix3(Vu, M1, M2, M3, M4);

v = L\r;

V = zeros(M1, M1);

N1 = M1*M2;
N2 = M3*M4;

for i = 1 : 1 : M2
    for j = 1 : 1 : M1
        % retrieve south points
        V(i,j) = v((i-1)*M1+j);
        % retrieve north points
        V(M2+M4+i,j) = v(N1+N2+(i-1)*M1+j);
    end
end

for i = 1 : 1 : M4
    for j = 1 : 1 : M3
        % retrieve west points
        V(M2+i, j) = v(N1+(i-1)*M3+j);
        % retrieve east points
        V(M2+i, M1-M3+j) = v(2*N1+N2+(i-1)*M3+j);
    end
end

% plot the results
figure(4);
surf(x, y, V);
colorbar;
title('Laplace Equation : Finite Differences Method Solution');
xlabel('x axis');
ylabel('y axis');
zlabel('V(x,y)');

clear all; clc;
