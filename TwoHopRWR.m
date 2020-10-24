function p42=TwoHopRWR(A) %A: adjacency matrix

cc=90;%Number of nodes
w=pinv(diag(sum(A,2)))*A;

c1=0;%The number of 1-hop neighbors
c2=0;%The number of 2-hop neighbors
n1=zeros(1,cc);
A_num=zeros(cc,cc);
for i=1:cc
    A_num(i,:)=(A(i,:)~=0);
    n1(i)=sum(A_num(i,:));
    c1=c1+n1(i);
end

for i=1:cc
    for j=1:cc
    A_2=A*A;
    if A_2(i,j)~=0
        c2=c2+1;
    end
    end
end
c2=c2-cc;

c=0.85; %pagerank
pp4=zeros(cc,cc);
for k=1:cc
p40=zeros(1,cc);
p40(k)=1; 
p41=p40;
pp4(k,:)=c*((c1/(c1+c2))*p40*w+(c2/(c1+c2))*p40*w^2)+(1-c)*p40;
while max(abs(pp4(k,:)-p41))>10^(-6)
    p41=pp4(k,:);
    pp4(k,:)=c*((c1/(c1+c2))*p41*w+(c2/(c1+c2))*p41*w^2)+(1-c)*p40;
end
end
p42=zeros(1,cc);
for i=1:cc
   p42(i)=sum(pp4(:,i));
end