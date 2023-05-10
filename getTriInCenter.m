function incenter=getTriInCenter(pts)
    A=pts(1,:);
    B=pts(2,:);
    C=pts(3,:);
    
    a=norm(B-C);
    b=norm(C-A);
    c=norm(A-B);
    
    incenter=(a/(a+b+c))*A+...
        (b/(a+b+c))*B+...
        (c/(a+b+c))*C;
end