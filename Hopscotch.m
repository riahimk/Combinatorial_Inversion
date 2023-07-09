function [z] = combinatorix(I,J,ElimMax)

% [1, -N1_2, N1_2*N2_3 - N1_3, N1_2*N2_4 - N1_4 + N1_3*N3_4 - N1_2*N2_3*N3_4]
% [0,     1,            -N2_3,                              N2_3*N3_4 - N2_4]
% [0,     0,                1,                                         -N3_4]
% [0,     0,                0,                                             1]
if nargin <1
J= 4;
I= 1;
ElimMax=power(2,16);
end

assert(J>1);
assert(I<J);

% x=[];
% for i=I:J-1
%     for j=i+1:J
%         fprintf("(%d,%d), ",i,j);
%         x=[x [i;j]];
%     end
% end
% fprintf("\n\n");

z0=[I:1:J];
Elim=length(z0)-2;
totalcmb = power(2,Elim)-1;
z=[];%[[];[]];%ones(totalcmb,1)*z0;

set=[I+1:J-1];
ii=0;
for elim=1:Elim 
    c = nchoosek(set,elim);
    m=size(c,1);
    %display(c);
    for i=1:m
        %%lia=1-double(ismember(z0,c(i,:)));
         ii=ii+1;       
         %%z(ii,:) = lia.*z0;
         comb=setxor(z0,c(i,:));
         comb(numel(z0))=0;
         z(ii,:) = comb;
        %fprintf('line %d done\n',ii);
    end
%       disp(z), pause
end  
 z=[z0;z];

% fprintf('\n');
% disp(z)
end
