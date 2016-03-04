function [A, b, N] = extract_matrix3(Vs, M1, M2, M3, M4)
    N1 = M1*M2;
    N2 = M3*M4;
    N = 2*(N1 + N2);
    A = zeros(N, N);
    
    A1 = laplacian(N1, M1); % laplacian for "north and ""south" row sets  
    for i = 1 : 1 : N1
        for j = 1 : 1 : N1
            % copy A1-south to A (offset 0)
            A(i,j) = A1(i,j);
            % copy A1-north to A (offset N1+N2)
            A(N1+N2+i, N1+N2+j) = A1(i,j);
        end
    end
    clear A1;
    
    A2 = laplacian(N2, M3); % laplacian for "west" and "east" row sets
    for i = 1 : 1 : N2
        for j = 1 : 1 : N2
            % copy A2-west to A (offset N1)
            A(N1+i, N1+j) = A2(i,j);
            % copy A2-east to A (offset 2*N1+N2)
            A(2*N1+N2+i, 2*N1+N2+j) = A2(i,j);
        end
    end
    clear A2;
    
    % dependances between A1-south and A2-west
    for i = 1 : 1 : M2
        A(N1-M1+i, N1+i) = 1; % for the upper part
        A(N1+i, N1-M1+i) = 1; % for the lower part
    end
    
    % dependences between A2-west and A1-north
    for i = 1 : 1 : M2
        A(N1+N2-M3+i, N1+N2+i) = 1; % for the upper part
        A(N1+N2+i, N1+N2-M3+i) = 1; % for the lower part
    end
    
    % dependences between A1-north and A2-east
    for i = 1 : 1 : M2
        A(N1+N2+M1-M2+i, 2*N1+2*N2-M2+i) = 1; % for the upper part
        A(2*N1+2*N2-M2+i, N1+N2+M1-M2+i) = 1; % for the lower part
    end
    
    % dependences between A1-south and A2-east
    for i = 1 : 1 : N2
        A(N1-M2+i, 2*N1+N2+i)  = 1; % for the upper part
        A(2*N1+N2+i, N1-M2+i) = 1; % for the lower part
    end    

    b = zeros(N, 1);
    
    % RHS values corresponding to A1-south matrix
    b(1) = -2*Vs ;
    for i = 2 : 1 : M1-1
        b(i) = -Vs;
    end
    b(M1) = -2*Vs;
    for i = 2 : 1 : M2
        b((i-1)*M1+1) = -Vs;
        b(i*M1) = -Vs;
    end
    
    % RHS values corresponding to A2-west matrix
    for i = 1 : 1 : M4
        b(N1+(i-1)*M3+1) = -Vs;
    end
    
    % RHS values corresponding to A1-north matrix
    for i = 1 : 1 : M2-1
        b(N1+N2+(i-1)*M1+1) = -Vs;
        b(N1+N2+i*M1) = -Vs;
    end
    b(N1+N2+(M2-1)*M1+1) = -2*Vs;
    for i = 2 : 1 : M1-1
       b(N1+N2+(M2-1)*M1+i) = -Vs; 
    end
    b(2*N1+N2) = -2*Vs;
    
    % RHS values corresponding to A2-east matrix
    for i = 1 : 1 : M4
        b(2*N1+N2+i*M3) = -Vs;
    end
    
end

% +------------------------+
% |  A1s   Dws         Des |
% |  Dsw   A2w   Dnw       |
% |        Dwn   A1n   Den |
% |  Dse         Dne   A2e |
% +------------------------+