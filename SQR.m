% function [S,Q,R]=SQR(A)
% 
%  if nargin<1
%     A=[3,1,-3;1,1,2;0,0,2];
%     A=rand(3);
%  end
% [n,~]=size(A);
% 
% Q=zeros(n);
% R=zeros(n);
% S=eye(n);
% D=eye(n);
% 
% Q(:,1)=A(:,1)/R(1,1);
% for j=1:n 
%     q=A(:,j);
%     for i=1:j-1
%         R(i,j) = Q(:,i)'*A(:,j);
%         q      = q - R(i,j)*Q(:,i);
%     end
%      R(j,j) = sqrt(q'*q);
%      D(j,j) = 1/R(j,j);
%      for i=1:j-1
%         S(i,j) = - D(j,j)*S(i,i:j-1)*R(i:j-1,j) ;
%      end
%          Q(:,j) = q*D(j,j);
% end
% S=D*S;
% return; 
% end
function [S, Q, R] = SQR(A)
    if nargin < 1
        A = rand(3); % Example matrix if no input is provided
    end

    [n, ~] = size(A);
    Q = zeros(n);
    R = zeros(n);
    S = eye(n);
    D = eye(n);

    Q(:, 1) = A(:, 1) / R(1, 1);
    for j = 1:n
        q = A(:, j);
        for i = 1:j-1
            R(i, j) = Q(:, i).' * A(:, j);
            q = q - R(i, j) * Q(:, i);
        end
        R(j, j) = norm(q);
        D(j, j) = 1 / R(j, j);
        S(1:j-1, j) = -D(j, j) * S(1:j-1, 1:j-1) * R(1:j-1, j);
        Q(:, j) = q * D(j, j);
    end
    S = D * S;
end
