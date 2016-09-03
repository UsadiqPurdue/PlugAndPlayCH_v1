function [y]=circfilter2(l,x,f_size)
    
    s=(f_size-1)/2;
    
    [m n]=size(x);
    
    u=zeros(m+(2*s),n+(2*s));
    
    u(s+1:m+s,s+1:n+s)=x;
    
    u(s+1:m+s,1:s)=x(:,n-s+1:n);
    u(s+1:m+s,n+s+1:n+(2*s))=x(:,1:s);
    
    u(1:s,s+1:n+s)=x(m-s+1:m,:);
    u(m+s+1:m+(2*s),s+1:n+s)=x(1:s,:);
    
    u(1:s,1:s)=x(m-s+1:m,n-s+1:n);
    u(m+s+1:m+(2*s),n+s+1:n+(2*s))=x(1:s,1:s);
    
    u(1:s,n+s+1:n+(2*s))=x(m-s+1:m,1:s);
    u(m+s+1:m+(2*s),1:s)=x(1:s,n-s+1:n);
    
    f=filter2(l,u); 
    % f=conv2(u,l);
    
    y=f(s+1:m+s,s+1:n+s);

end
