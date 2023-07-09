function C = strassen(A, B)
    [m, ~] = size(A);
    
    % Base case: If matrix size is 1x1, perform simple multiplication
    if m == 1 || mod(m,2)==0
        C = A * B;
        return;
    end
    
    % Split matrices into submatrices
    mid = floor(m/2);
    A11 = A(1:mid, 1:mid);
    A12 = A(1:mid, mid+1:end);
    A21 = A(mid+1:end, 1:mid);
    A22 = A(mid+1:end, mid+1:end);
    B11 = B(1:mid, 1:mid);
    B12 = B(1:mid, mid+1:end);
    B21 = B(mid+1:end, 1:mid);
    B22 = B(mid+1:end, mid+1:end);
    
    % Recursive steps
    size(A11),size(A22),size(B11),size(A22),pause
    P1 = strassen(A11 + A22, B11 + B22);
    P2 = strassen(A21 + A22, B11);
    P3 = strassen(A11, B12 - B22);
    P4 = strassen(A22, B21 - B11);
    P5 = strassen(A11 + A12, B22);
    P6 = strassen(A21 - A11, B11 + B12);
    P7 = strassen(A12 - A22, B21 + B22);
    
    % Combine submatrices to obtain the result
    C11 = P1 + P4 - P5 + P7;
    C12 = P3 + P5;
    C21 = P2 + P4;
    C22 = P1 - P2 + P3 + P6;
    
    % Assemble the result matrix
    C = [C11, C12; C21, C22];
end
