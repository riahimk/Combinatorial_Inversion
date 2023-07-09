function [L,U,iL,iU,iA]=BRSI(A,nblock,ncard)
% Description: This function computes the UL factorization for a given square matrix A of size nxn. It also computes the LU factorization for the inverse of A.
% Inputs:
%   - A: The input square matrix of size nxn.
% Outputs:
%   - L: Lower triangular matrix of the UL factorization of A.
%   - U: Upper triangular matrix of the UL factorization of A.
%   - iL: Lower triangular matrix of the LU factorization of the inverse of A.
%   - iU: Upper triangular matrix of the LU factorization of the inverse of A.
% Author: Mohamed Kamel Riahi
% Affiliation: Mathematics Department, Khalifa University
% Date: [01/02/2023]


% if nargin <1
%     A = randi(12,32,32);  % Example matrix
% end


n  = size(A,1);

% P=eye(n);
% [M,P,PA] = permute_cols(A);
% M=A;

if n<=16
    [L,U,iL,iU,iA]=RSI(A);
    return
end

% [nblock,ncard]=parametertization(n);
% ncard=2;
% nblock = ceil(n / ncard);

I  = eye(n);
Up = I;

[Ls,Us]=blockSplitMatrix_ncard(A,nblock,ncard);

for i=1:ncard
    Up=Up*Us;
    iU=INVERTE(Us,nblock,ncard);
    A_=I+iU*Ls;
    [Ls,Us]=blockSplitMatrix_ncard(A_,nblock,ncard);
end
U  = Up;
L  = I+iU*Ls;
iU = INVERTE(Up,nblock,ncard);
iL = INVERTE(L',nblock,ncard)';

iA = iL*iU;

    function [R]=INVERTE(matrix,nb,nc)
            if istril(matrix) || istriu(matrix)
                R=COMBRIT(matrix,nb,nc);
            else
                R=BlockInverseT(matrix,nb,nc);
            end
    end
    
    function [TibU] = BlockInverseT(matrixbT,nb,nc)
        TibU=eye(size(matrixbT));
        iD  =zeros(size(matrixbT));  
        for iu =1:nc
            %length([1+(iu-1)*nb:iu*nb]),pause
            Dii = matrixbT(1+(iu-1)*nb:iu*nb,1+(iu-1)*nb:iu*nb);
            %[size(Dii,1),size(matrixbT,1)]
            [~,~,~,~,iDii]=BRSI(Dii,nb/nc,nc);
            iD(1+(iu-1)*nb:iu*nb,1+(iu-1)*nb:iu*nb)=iDii;
            for ju=iu+1:nc
                TibU(1+(iu-1)*nb:iu*nb,1+(ju-1)*nb:ju*nb)= ...
                    iDii*matrixbT(1+(iu-1)*nb:iu*nb,1+(ju-1)*nb:ju*nb);
            end
        end
        ibU=COMBRIT(TibU,nb,nc);
        TibU=ibU*iD;
    end


end