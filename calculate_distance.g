
x:=Indeterminate(GF(2),"x");

#H1 from h1x, H2^T from h2xt





#hx:=1+x+x^2+x^3+x^4+x^5+x^6+x^7+x^8+x^9+x^10;
#n0:=11;



#functions:
############################################


        
#cyclic code check matrix from hx with x^n+1
        
cyclic_H:=function(hx,n)
        local hx_vector_gf2, H, k,r,i1,i2,i3,h;
        
        H:=[];

  

        hx_vector_gf2:=CoefficientsOfUnivariatePolynomial(hx);
       
        k:=Length(hx_vector_gf2)-1;
        r:=n-k;

        for i1 in [1..r] do
            H[i1]:=[];
            if (n-k-i1)<>0 then
                for i2 in [1..(n-k-i1)] do Add(H[i1], 0*Z(2)); od;
            fi;

            h:=x^(i1-1)*hx;
            for i3 in [1..(k+i1)] do Add(H[i1], CoefficientsOfUnivariatePolynomial( h)[k+i1+1-i3]); od;
           

            #Print("for i1=",i1,"\n");
            #Print(G[i1],"\n");
                
        od;
        #Print(DimensionsMat(H));
    
        #if n<>DimensionsMat(H)[2] then
        #    Print("error: cyclic_H: \n  n=",n,"  column number=", DimensionsMat(H)[2], "\n");
        #fi;

        #Display(H);
        return H;

#H looks like:
        
#   000...h_k...h_0
#   .
#   .
#   .
#   h_k...h_0...000

#notice: CoefficientsOfUnivariatePolynomial:  1+2x+3x^3-> 1,2,0,3

        
end;

#not finished yet;
h_to_g:=function(hx,n)
        

end;

cyclic_G:=function(hx,n)
        local gx,fx,gx_vector_gf2, G,g, k,r,i1,i2;
        fx:=x^(n)+1;
        G:=[];

        gx:=Quotient(fx,hx);

        #Print("n=",n,"\n","gx=");
        #Print(gx,"\n");

        gx_vector_gf2:=CoefficientsOfUnivariatePolynomial(gx);
       
        r:=Length(gx_vector_gf2)-1;
        k:=n-r;

        for i1 in [1..k] do

            g:=x^(i1-1)*gx;
            G[i1]:=ShallowCopy(CoefficientsOfUnivariatePolynomial( g));

            #Print("for i1=",i1,"\n");
            #Print(G[i1],"\n");

            for i2 in [1..(k-i1)] do Add(G[i1], 0*Z(2)); od;
                
        od;

#check:
        if cyclic_H(hx,n)*TransposedMat(G)<>NullMat( r, k,Z(2) ) then
            Print("error:cyclic_G \n");
        fi;

        #Display(G);

        return G;
end;

#get Hz from H1, H2^T and n1

HP_Hz:=function(h1x,h2xt,n1)
        local di,n2,i,H1,H1t,H2,H2t,E1t,E2t,Hz_left,Hz_right;

#H1t: H1_Transposed
#E1t: E1_tilt

        H1:=cyclic_H(h1x,n1);
        H2t:=cyclic_H(h2xt,n1);
        H1t:=TransposedMat(H1);
        H2:=TransposedMat(H2t);

        di:=DimensionsMat( H2);
        n2:=di[2];
        
        E1t:=IdentityMat(n1,Z(2));
        E2t:=IdentityMat(n2,Z(2));
        Hz_left:=KroneckerProduct(H2t,E1t);
        Hz_right:=KroneckerProduct(E2t,H1t);
        
        if (Length(Hz_left)=Length(Hz_right) )  then
            for i in [1..Length(Hz_left)] do
                Append(Hz_left[i],Hz_right[i]);
            od;
            #Print("\n hz for n1=",n1," is \n");
            #Display(Hz_left);
            return Hz_left;

        else
            Print("error :function HP_Hz \n");
            return 0;
        fi;
end;

#get Hz from H1, H2 and n1

HP_Hz2:=function(h1x,h2x,n1)
        local di,n2,i,H1,H1t,H2,H2t,E1t,E2t,Hz_left,Hz_right;

#H1t: H1_Transposed
#E1t: E1_tilt

        H1:=cyclic_H(h1x,n1);
        H2:=cyclic_H(h2x,n1);
        H1t:=TransposedMat(H1);
        H2t:=TransposedMat(H2);

        di:=DimensionsMat( H2);
        n2:=di[2];
        
        E1t:=IdentityMat(n1,Z(2));
        E2t:=IdentityMat(n2,Z(2));
        Hz_left:=KroneckerProduct(H2t,E1t);
        Hz_right:=KroneckerProduct(E2t,H1t);
        
        if (Length(Hz_left)=Length(Hz_right) )  then
            for i in [1..Length(Hz_left)] do
                Append(Hz_left[i],Hz_right[i]);
            od;

            #Print("\n hz for n1=",n1," is \n");
            #Display(Hz_left);

            return Hz_left;

        else
            Print("error :function HP_Hz2 \n");
            return 0;
        fi;
end;



HP_Hx:=function(h1x,h2xt,n1)
        local r1,r2,di,n2,i,H1,H1t,H2,H2t,E1,E2,Hx_left,Hx_right;


        H1:=cyclic_H(h1x,n1);
        H2t:=cyclic_H(h2xt,n1);
        H1t:=TransposedMat(H1);
        H2:=TransposedMat(H2t);
        di:=DimensionsMat( H2);
        n2:=di[2];
        r1:=Length(H1);
        r2:=Length(H2);
        
        E1:=IdentityMat(r1,Z(2));
        E2:=IdentityMat(r2,Z(2));
        Hx_left:=KroneckerProduct(E2,H1);
        Hx_right:=KroneckerProduct(H2,E1);
        
        if (Length(Hx_left)=Length(Hx_right) )  then
            for i in [1..Length(Hx_left)] do
                Append(Hx_left[i],Hx_right[i]);
            od;
            #Print("\n hx for n1=",n1," is \n");
            #Display(Hx_left);
            return Hx_left;

        else
            Print("error :function HP_Hx \n");
            return 0;
        fi;
end;

HP_Hx2:=function(h1x,h2x,n1)
        local r1,r2,di,n2,i,H1,H1t,H2,H2t,E1,E2,Hx_left,Hx_right;


        H1:=cyclic_H(h1x,n1);
        H2:=cyclic_H(h2x,n1);
        H1t:=TransposedMat(H1);
        H2t:=TransposedMat(H2);
        di:=DimensionsMat( H2);
        n2:=di[2];
        r1:=Length(H1);
        r2:=Length(H2);
        
        E1:=IdentityMat(r1,Z(2));
        E2:=IdentityMat(r2,Z(2));
        Hx_left:=KroneckerProduct(E2,H1);
        Hx_right:=KroneckerProduct(H2,E1);
        
        if (Length(Hx_left)=Length(Hx_right) )  then
            for i in [1..Length(Hx_left)] do
                Append(Hx_left[i],Hx_right[i]);
            od;
            #Print("\n hx for n1=",n1," is \n");
            #Display(Hx_left);
            return Hx_left;

        else
            Print("error :function HP_Hx2 \n");
            return 0;
        fi;
end;


rand_dist:=function(G,perm_num)
        local k,n,r,i,i2,perm_mat,perm_G,gaussian,d,wt,perm_v,range_v,rand_num;
     
        
        k:=Length(G);
        n:=Length(G[1]);

        r:=RankMat(G);


        
  
        d:=n;
        
      
        perm_G:=ShallowCopy(G);
        for i in [1..perm_num] do

            range_v:=[1..n];
            perm_v:=[];

            for i2 in [1..n] do
                rand_num:=Random([1..(n-i2+1)]);
                    Add(perm_v,range_v[rand_num]);
                    Remove(range_v,rand_num);
            od;

                perm_mat:=PermutationMat(  PermList(perm_v), n, GF(2) );
                perm_G:=perm_G*perm_mat;
                gaussian:=TriangulizedMat( perm_G );
                #Print("\n",i,"\n");
                #Display(gaussian);
                #Print("\n");
         


            wt:=WeightVecFFE( gaussian[r]);

           # for i2 in [1..k] do
           #     i3:=k+1-i2;
           #     if WeightVecFFE( gaussian[i3])<>0 then  wt:=WeightVecFFE( gaussian[i3]);break;fi;
           #od;
            
            if wt<d then d:=wt;  fi;
        od;
        return d;
end;

        

int_to_gf2:= function(int_vector)
        
              local l,gf2_vector,i;
              l:=Length(int_vector);
              gf2_vector:=[];

              for i in[1..l] do
                             if int_vector[i]=1 then
                                gf2_vector[i]:=Z(2)^0;
                            elif int_vector[i]=0 then
                                 gf2_vector[i]:=0*Z(2);
                            else
                                Print("error, this is not a binary  vector");
                            fi;

               od;
               return gf2_vector;

end;

cyclic_dist:=function(hx,n,perm_n)
        local d,G;
        G:=cyclic_G(hx,n);
        d:=rand_dist(G,perm_n);
        return d;

end;


CSS_dist:=function(Hx,Hz,perm_num)
        local G,d;
        G:=ShallowCopy(Hx);
        Append(G,Hz);

        #Display(G);

        d:=rand_dist(G,perm_num);
        return d;
end;


calculate_dist:=function(G,k,n)
        local bi_vec,gf2_vec,temp,i1,i2,d,wt,codeword;
         d:=n;
         bi_vec:=[];
         gf2_vec:=[];

         for i1 in [1..(2^(k)-1)] do
            bi_vec[1]:=i1 mod 2;
            temp:=(i1-bi_vec[1])/2;
    

            for i2 in [2..k] do
    
                bi_vec[i2]:=temp mod 2;
                temp:= (temp-bi_vec[i2])/2;
            od;

            gf2_vec:=ShallowCopy(int_to_gf2(bi_vec));
            codeword:=gf2_vec*G;
            wt:=WeightVecFFE( codeword);
            if wt<d then d:=wt;  fi;
          od;
          return d;
end;
########################################

#where max_n=n0*max_m
#use H1,H2^2:

QHP_code:=function(h1x,h2xt,n0,max_m,perm_num)
    local i,n1,n2t,H1,H2,H1t,H2t,di_1,di_2,n2,r1,r2,s1,s2,k1,k2,n,k,Hx,Hz,d_upper_bound,d1,d2,d1_t,d2_t,G1,G2,G1t,G2t;
    for i in [1..max_m] do
        n1:=i*n0;
        n2t:=n1;
#they are the row_numbers for H1, H2^T

        H1:=cyclic_H(h1x,n1);
        H2t:=cyclic_H(h2xt,n1);
        H1t:=TransposedMat(H1);
        H2:=TransposedMat(H2t);
        G1:=cyclic_G(h1x,n1);
        G2t:=cyclic_G(h2xt,n1);
        G1t:=TransposedMat(G1);
        G2:=TransposedMat(G2t);
        di_1:=DimensionsMat( H1);
        di_2:=DimensionsMat( H2);
        n1:=di_1[2];
        n2:=di_2[2];
        r1:=di_1[1];
        r2:=di_2[1];

        d1:=rand_dist(G1,perm_num);
        d2:=rand_dist(G2,perm_num);
        d1_t:=rand_dist(G1t,perm_num);
        d2_t:=rand_dist(G2t,perm_num);

#notice: r1,r2 are not the ranks

        s1:=n1-r1;
        k1:=n1-RankMat(H1);
        s2:=n2-r2;
        k2:=n1-RankMat(H2);
    
        n:=r2*n1+r1*n2;
        k:=2*k1*k2-k1*s2-k2*s1;

        Hx:=HP_Hx(h1x,h2xt,n1);
        Hz:=HP_Hz(h1x,h2xt,n1);

        #Display(Hx*TransposedMat(Hz));
        #Print(DimensionsMat(Hx*TransposedMat(Hz)));
        d_upper_bound:=CSS_dist(Hx,Hz,perm_num);;
  
        Print("H1=",h1x, "H2^T=",h2xt, "   mod  x^(",n0,"*",i,")+1 : ", "d1=", d1," d2=",d2," d1_t=" ,d1_t, " d2_t=", d2_t,  "\n","n=",n,", k=",k," ,d<=",d_upper_bound,"\n");

    #Print("for n1=",n1," G1=\n");
    #Display(G1);
  
    od;
end;

#from H1, H2:
QHP_code2:=function(h1x,h2x,n0,max_m,perm_num)
    local i,n1,n2t,H1,H2,H1t,H2t,di_1,di_2,n2,r1,r2,s1,s2,k1,k2,n,k,Hx,Hz,d_upper_bound;
    for i in [1..max_m] do
        n1:=i*n0;
        n2t:=n1;

        H1:=cyclic_H(h1x,n1);
        H2:=cyclic_H(h2x,n1);
        H1t:=TransposedMat(H1);
        H2t:=TransposedMat(H2);
        di_1:=DimensionsMat( H1);
        di_2:=DimensionsMat( H2);
        n1:=di_1[2];
        n2:=di_2[2];
        r1:=di_1[1];
        r2:=di_2[1];

        s1:=n1-r1;
        k1:=n1-RankMat(H1);
        s2:=n2-r2;
        k2:=n1-RankMat(H2);
    
        n:=r2*n1+r1*n2;
        k:=2*k1*k2-k1*s2-k2*s1;

        Hx:=HP_Hx(h1x,h2xt,n1);
        Hz:=HP_Hz(h1x,h2xt,n1);

        #Display(Hx*TransposedMat(Hz));
        #Print(DimensionsMat(Hx*TransposedMat(Hz)));
        d_upper_bound:=CSS_dist(Hx,Hz,perm_num);;
  
        Print("for H1 from h1x=",h1x, "and H2^T from h2xt=",h2xt, "  with x^(",n0,"*",i,")+1 we have:\n","n=",n,", k=",k," ,d<=",d_upper_bound,"\n");

    #Print("for n1=",n1," G1=\n");
    #Display(G1);
  
    od;
end;

#main function:
####################################################################################################



h1_1:=1+x+x^3;
n0_1:=7;
h2t_1:=1+x;


perm_num:=30;
max_m:=5;
Print("perm_num=",perm_num,"  \n \n");
QHP_code(h1_1,h2t_1,n0_1,max_m,perm_num);


#for i in [1..max_m] do
#    d1:=cyclic_dist(h1_1,n0_1*i,perm_num);

#    d2_t:=cyclic_dist(h2t_1,n0_1*i,perm_num);
#
 #   Print("For n=",n0_1*i,"  d1=",d1,"  d2_t=",d2_t,"\n");

#od;






