function [Xdot]=ecuDif(t,x)
  m=0.5;
  b=0.1;
  L=0.5;
  g=9.81;
  k=b/(m*L);
 Xdot=zeros(2,1); //#Xdot
 Xdot(1,1)=x(2);
 Xdot(2,1)=-(g/L)*sin(x(1))-(k/m)*x(2);
endfunction
