function [A, b, N] = extract_matrix2(Vu, M1, M2, M3, M4)
    N1 = M1*M2;
    A1 = laplacian(N1, M1);
    
    N2 = M3*M4;
    A2 = laplacian(N2, M3);
    
    N = N1 + N2;
    A = zeros(N, N);
    
    % copy matrix A1 to A
    for i = 1 : 1 : N1
        for j = 1 : 1 : N1
           A(i,j) = A1(i,j); 
        end
    end
    clear A1;
    
    % copy matrix A2 to A
    for i = 1 : 1 : N2
        for j = 1 : 1 : N2
           A(i+N1, j+N1) = A2(i,j); 
        end
    end
    clear A2;
    
    % add the additional ones
    for i = 1 : 1 : M3
        A(N1+i, N1+i-M1) = 1;
        A(N1+i-M1, N1+i) = 1;
    end
    
    b1 = zeros(N1, 1);
    % the rightest point of each row reaches at potential Vu
    for i = 1 : 1 : M2-1
        b1(i*M1) = -Vu;
    end
    % exept of the last row
    for i = (M2+1) : 1 : M1 
        b1((M2-1)*M1+i) = -Vu;
    end
    b1(N1) = -2*Vu; % the rightest point reaches at potential Vu+Vu = 2Vu
    
    b2 = zeros(N2, 1);
    % the rightest point of each row reaches at potential Vu
    for i = 1 : 1 : M4-1
        b2(i*M3) = -Vu;
    end
    %the last row reaches potential Vu
    for i = N2-M3+1 : 1 : N2
        b2(i) = -Vu;
    end
    b2(N2) = -2*Vu; % the rightest point reaches at potential Vu+Vu = 2Vu
    
    b = zeros(N, 1);
    for i = 1 : 1 : N1
        b(i) = b1(i);
    end
    clear b1;
    
    for i = 1 : 1 : N2
        b(i+N1) = b2(i);
    end
    clear b2;

end