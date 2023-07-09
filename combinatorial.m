function [S]= combinatorial(T,nb,nc)
if nargin<1
    T = triu(rand(8));
    nb=2;
    nc=1;
end

S=zeros(nb,nb,nc,nc);
for i_=1:nc
    S(:,:,i_,i_)=eye(nb,nb);
end
for jj_sc=2:nc % could be a parfor loop
    sc=Hopscotch(1,jj_sc,1);
    j_sc=size(sc,2);
    for kd= 0:nc-j_sc
        S(:,:,1+kd,j_sc+kd) = INVERSE_ENTRY_VALUES_USE_CARD(T,sc,kd);
    end
end
end
