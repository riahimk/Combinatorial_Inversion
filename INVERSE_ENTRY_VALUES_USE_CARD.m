function [v]=INVERSE_ENTRY_VALUES_USE_CARD(M,sc_k,trans)
v=zeros(size(M,1,2));
for ii=1:size(sc_k,1)
    pw=0;
    c1=eye(size(M,1,2));
    for jj=1:size(sc_k,2)-1
        if(sc_k(ii,jj+1)~=0)
            %fprintf('index are (%d,%d) \t',sc_k(i,j)+trans,sc_k(i,j+1)+trans);
            c1=c1*M(:,:,sc_k(ii,jj)+trans,sc_k(ii,jj+1)+trans);
            pw=pw+1;
        end
    end
    %fprintf('\n');
    v=v+power(-1,pw)*c1;
end
end