function g = CurvaTest(t)
    f=-2*t;
    f=2.^f;
    f=f.*sin(20*t)+0.1*sin(250*t);
    hd=(1+(t.^2));
    g=1./hd;
    f=f+g;
    f=g;
end