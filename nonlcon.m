function [c,ceq] = nonlcon(x)
ceq = [  x(1) * x(4) - 4*2;
         x(2) * x(4) - 5*1;
         x(3) * x(4) - 7*1;

         x(1) * x(5) - 4*3;
         x(2) * x(5) - 5*1.5;
         x(3) * x(5) - 7*1.5;

         x(1) * x(6) - 4*1;
         x(2) * x(6) - 5*0.5;
         x(3) * x(6) - 7*0.5;
         
          x(1) * x(7) - 4*1;
         x(2) * x(7) - 5*0.5;
         x(3) * x(7) - 7*0.5;
 c = [];