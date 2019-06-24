function y = safelog(x)
    
    x(x==0) = realmin;
    y = log(x);