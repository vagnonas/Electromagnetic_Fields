function A = laplacian(N, M)

    % create matirx A
    A = zeros(N, N);
    for i = 1 : 1 : N
        for j = 1 : 1 : N
            d = abs(i-j);
            % diagonal elements
            if d == 0
                A(i, j) = -4;
            end       
            % above and below main diagonal elements
            if d == 1
                if mod(i,M) + mod(j,M) == 1
                    A(i, j) = 0;
                else
                    A(i, j) = 1;
                end
            end
            % rest of the elements
            if d == M
                A(i, j) = 1;
            end
        end
    end
        
end
