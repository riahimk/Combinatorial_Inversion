function [L,U] = blockSplitMatrix_ncard(matrix,nblock, ncard)
    % Get the size of the matrix
    [m, n] = size(matrix);
    if m~=n
        error('matrix should be square!\n');
    end
    
    
    L=zeros(size(matrix));
    % Split the matrix into blocks
    for i = 1:ncard
        for j = 1:ncard
            % Define the row and column ranges for the current block
            rowStart = (i-1) * nblock + 1;
            rowEnd = min(i * nblock, m);
            colStart = (j-1) * nblock + 1;
            colEnd = min(j * nblock, n);
            if i<=j
            % Extract the current block from the matrix
            U(rowStart:rowEnd, colStart:colEnd) = matrix(rowStart:rowEnd, colStart:colEnd);
            elseif i>j
            % Assign the block to the corresponding index in the output matrices
            L(rowStart:rowEnd, colStart:colEnd) = matrix(rowStart:rowEnd, colStart:colEnd);
            end
        end
    end
end