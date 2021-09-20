which_field:=function(p,n)
        local m,max_m;
        m:=1;
        max_m:=1000;
        while (((p^m-1) mod n)<>0 ) and m<max_m do
            m:=m+1;
        od;

        if m<1000 then
            return m;
        else
            Print("no ",n,"th root of unity over GF ",p,"^m, when m<",max_m);
            return -1;
        fi;
end;

            

nth_proot_of_unity:=function(p,n)
        local m,i1,i2,i3,alpha,k;
      
        m:=which_field(p,n);

        alpha:=Z(p,m);

        for i1 in [1..(2^m-2)] do
            k:=0;

            if  alpha^(i1*n)=Z(p,m)^0 then

                for i2 in [1..(n-1)] do
                    if alpha^(i1*i2)=0 then
                        k:=1;
                    fi;
                    
                od;

                if k=0 then
                     return alpha^i1;
                fi;

            fi;
        od;
Print("no ",n,"th root of unity in GF ",p,"^m");
return 0;


end;




max_i:=6;
min_i:=1;
p:=2;
      
alpha:=[];


for i in 2*[min_i..max_i]+1 do
       #Print(i,"   ",nth_proot_of_unity(p,i),"\n");
        alpha[i]:=nth_proot_of_unity(p,i);
od;

#when n is small than 27, it is fast to calculate nth root of unity in my computer. cpu is i7 9750H.

Print(" For generators g=1+x+y+x*y \n");
for n1 in  2*[min_i..max_i]+1 do
    for n2 in 2*[min_i..max_i]+1 do
        r:=0;
        cz:=[];
      
   
        for i in [1..n1] do
            for j in [1..n2] do
                x:=alpha[n1]^i;
                y:=alpha[n2]^j;
            
                 if 1+x+y+x*y=0*Z(p) then
                 #if (x+1)*(y^2+x^2*y+1)=0*Z(p) and (x^2+x+1)*(y+1)=0*Z(p) then
                    # this is an example in "A theory of two-dimensional cyclic codes" 
                    r:=r+1;
                    Append(cz,[[i,j]]);
                 fi;
            od;
        od;

        Print("for n1=",n1," n2=",n2,"  n= ",n1*n2," k=",n1*n2-r,"\n \n");
        #Print("the exponents of the common zeros are",cz,"\n");
    od;
od;

#gamma:=alpha[3];
#beta:=alpha[5];
#x:=gamma;
#y:=beta;
#Print((gamma+1)*(beta^2+gamma^2*beta+1)=0," ",(x^2+x+1)*(y+1)=0,(gamma+1)*(beta^2+gamma^2*beta+1),(x^2+x+1)*(y+1));
#Print(0*Z(2)=0, 0*Z(2)=0*Z(2^2), 0*Z(2^3)=0*Z(2^5));
