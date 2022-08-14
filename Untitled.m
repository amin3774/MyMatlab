mat=Modified_New_Eq_string;
n=32;
for i=1:n
    temp=char(mat(i));
    %temp=temp(1:89);
    temp=string(temp);
    mat(i)=temp;
end
for i=1:n
    mat(i)=subs(eval(mat(i))); 
end
for i=1:n
    mat(i)=strcat(mat(i),'&&');
end
for i=1:n
    temp=mat(i);
temp = regexprep(temp,'[_]','');
 mat(i)=temp;
end

 N=[6:11];
run_time=[0.4014 1.9360 9.9496 49.4758 262.291 1341.62];
plot(N,run_time,'b*')
grid on
xlim([0 12])
ylim([-2 1400])