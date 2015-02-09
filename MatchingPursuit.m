%f=a1*g1+...+ak*gk+Rk
%D is the dictionary as a matrix with columns g1, g2, ...(the atoms) 
%f must be a column vector
%The output A is a vector with entries a1, a2,...
%The output G is a matrix with columns g1, g2,...
%The output R is the first Rk so that ||Rk||^2<eps
%eps is the stopping threshold: ||Rk||^2<eps stops the recursion
%The algorithm works if f has the same size as any gk and all gk have norm
%1.
function [A, R, G]=MP(f, D, eps)
R=f; A=[]; G=[]; M=D; %we'll extract the gk from among the comumns of M. 
while (norm(R)>=eps) && (~isempty(M))
    sz=size(M); n=sz(2); %extracting the #of columns of M
    a=zeros(n); %initiallizing the vector of absolute inner products between R and the atoms of M.
    max=0; i=1;
    for k=1:n
        a(k)=sum(R.*M(:,k)); %computing each inner product
        if (max <= abs(a(k)))
            max=abs(a(k)); %selecting the maximum abs inner product
            i=k; %selecting the index of the atom which yields the max abs inner product
        end
    end
    B=A; 
    A=[B; a(i)]; %A grows with the new entry a(i) as a column vector 
    B=G;
    G=[B, M(:,i)]; %G grows with a new column
    R=R-a(i)*M(:,i);
    M(:,i)=[]; %dropping the selected column from M
end    
end
