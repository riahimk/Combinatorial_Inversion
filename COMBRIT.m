function [mat_invAij]=COMBRIT(A,sblock,n_card)
% Two-Block-Combinatorics-Recursive Inverse Triangular Matrix

if (nargin<3)
    if(nargin<1)
        sblock  = 40;
        n_card  = 5;
    end
    nA      = n_card * sblock;
    A       = triu(hilb(nA),1)+eye(nA);
end

isTriangular= (istril(A) || istriu(A));
if (istril(A) || istriu(A))==0
    error('The input matrix for COMBRIT must be triangular!\n');
end

if size(A,1)<=16
    mat_invAij = CRIT(A);
    %mat_invAij = inv(A);
    return;
end

%  S = zeros(size(A));

% assume we have a card for combinatorics up to n_card this will help to inverse
% triangular matrices of size (n_card)x(n_card)
tic
% card   = combinatorix(1,n_card,1);
%  n      = size(A,2);
% n_card = size(card,2);

%fprintf(' The size of the main matrix is n=%d\n and the size of the card is n_card=%d\n',n,n_card);

% Decompose the entry triangular matrix A into block matrices Aij
% sorted in a tensor format Aij(:,:,oi,oj)

Aij     = zeros(sblock,sblock,n_card,n_card);
Bij     = zeros(sblock,sblock,n_card,n_card);
invAij  = zeros(sblock,sblock,n_card,n_card);
invAij  = zeros(sblock,sblock,n_card,n_card);
invAij  = zeros(sblock,sblock,n_card,n_card);

for iu=1:n_card
    for ju=iu:n_card
        %fprintf("size (%d,%d)(%d,%d)\n",1+(i-1)*n_card,i*n_card,1+(j-1)*n_card,j*n_card);
        Aij(:,:,iu,ju) = A(1+(iu-1)*sblock:iu*sblock,1+(ju-1)*sblock:ju*sblock);
    end
end

%fprintf('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n');
%fprintf('~~ BlcokDiag Matrix inverse begin ~~\n');


for iv=1:n_card
    %invDij(:,:,i,i) = CRIT('',Aij,i,i,card,n_card);
    %invDij(:,:,iv,iv) = inv(Aij(:,:,iv,iv));
    invDij(:,:,iv,iv) = COMBRIT(Aij(:,:,iv,iv),size(Aij(:,:,iv,iv),1)/n_card,n_card);
    for jv=iv+1:n_card
        Bij(:,:,iv,jv) = invDij(:,:,iv,iv) * Aij(:,:,iv,jv);
    end
    Bij(:,:,iv,iv) = eye(sblock);
end

%fprintf('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~done!\n\n');
%fprintf('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n');
%fprintf('~~ Blcok Tensor-inverse begin ~~\n');

invBij=Inverse_card_tri_mat('B',Bij,0,0,n_card);

for jd=1:n_card
    for id=1:jd
    invAij(:,:,id,jd) = invBij(:,:,id,jd)*invDij(:,:,jd,jd);
    end
end

mat_invAij = Tensor2Matrix(invAij);


% fprintf('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~done!\n\n');

    function [M]=Tensor2Matrix(T)
        M=zeros(size(T,3,4).*size(T,1,2));
        for it=1:size(T,3)
            for jt=it:size(T,4)
                M(1+(it-1)*size(T,1):it*size(T,1),1+(jt-1)*size(T,2):jt*size(T,2)) = T(:,:,it,jt);
            end
        end
    end

    function [S]= Inverse_card_tri_mat(char,T,li,lj,nc)
        if char=='B'
            S=zeros(sblock,sblock,nc,nc);
            for ik=1:nc
                S(:,:,ik,ik)=eye(sblock,sblock);
            end
        else
            S=eye(size(T,1,2));
        end
        for jj_sc=2:nc % could be a parfor loop
            sc=sub_card(jj_sc);
            j_sc=size(sc,2);
            for kd= 0:nc-j_sc
                if char=='B'
                    S(:,:,1+kd,j_sc+kd) = val_comb_tri_matrix(char,T,li,lj,sc,kd);
                else
                    %fprintf('[kd=%d](%d,%d)\n',kd,1+kd,j_sc+kd);
                    S(1+kd,j_sc+kd) = val_comb_tri_matrix(char,T,li,lj,sc,kd);
                end
            end
        end
    end

    function [v]=val_comb_tri_matrix(char,M,li,lj,sc_k,trans)
        if char=='B'
            v=zeros(size(M,1,2));
        else
            v=0;
        end
        for ii=1:size(sc_k,1)
            pw=0;
            if char=='B'
                c1=eye(size(M,1,2));
            else
                c1=1;
            end
            for jj=1:size(sc_k,2)-1
                if(sc_k(ii,jj+1)~=0)
                    %fprintf('index are (%d,%d) \t',sc_k(i,j)+trans,sc_k(i,j+1)+trans);
                    if char=='B'
                        c1=c1*M(:,:,sc_k(ii,jj)+trans,sc_k(ii,jj+1)+trans);
                    else
                        c1=c1*M(sc_k(ii,jj)+trans,sc_k(ii,jj+1)+trans,li,lj);
                    end
                    pw=pw+1;
                end
            end
            %fprintf('\n');
            v=v+power(-1,pw)*c1;
        end
    end

    function [sc]=sub_card(n_sc)
        sc=Hopscotch(1,n_sc,1);
    end
end
