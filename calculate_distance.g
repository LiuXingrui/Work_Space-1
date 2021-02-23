
x:=Indeterminate(GF(2),"x");

#H1 from h1x, H2^T from h2xt

h1x:=1+x+x^3;
n0:=7;
h2xt:=1+x;



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
            Append(H[i1],CoefficientsOfUnivariatePolynomial( h));

            #Print("for i1=",i1,"\n");
            #Print(G[i1],"\n");
            if i1<>1 then
                for i3 in [1..(i1-1)] do Add(H[i1], 0*Z(2)); od;
            fi;

           
                
        od;
        Print(DimensionsMat(H));
    
        #if n<>DimensionsMat(H)[2] then
        #    Print("error: cyclic_H: \n  n=",n,"  column number=", DimensionsMat(H)[2], "\n");
        #fi;

        return H;

#H looks like:
        
#   000...h_k...h_0
#   .
#   .
#   .
#   h_k...h_0...000

#notice: CoefficientsOfUnivariatePolynomial:  1+2x+3x^3-> 3,0,2,1

        
end;

#don't use this function now, because it will give the same form H as the last function, just replace hx by gx,
#but we need 
#   g_0...g_r000...
#   0g_0...g_r00...
#   .
#   .
#   000...g_0...g_r
#02/23


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

        #Display(G);

        return G;
end;

#get Hz from H1, H2^T and n1

HP_Hz:=function(h1x,h2xt,n1)
        local di,n2,i,H1,H1t,H2,H2t,E1t,E2t,Hz_left,Hz_right;

#H1t: H1_Transposed, E1t: E1_tilt

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
                Append(Hz_left,Hz_right);
            od;
            return Hz_left;

        else
            Print("error :function HP_Hz \n");
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
                Append(Hx_left,Hx_right);
            od;
            return Hx_left;

        else
            Print("error :function HP_Hx \n");
            return 0;
        fi;
end;


rand_dist:=function(G)
        local k,n,perm_num,i,perm,perm_mat,perm_G,gaussian,d,wt,rand_a,rand_b;
        k:=Length(G);
        n:=Length(G[1]);
        d:=n;
        perm_num:=20;
        perm_G:=ShallowCopy(G);
        for i in [1..perm_num] do
            rand_a:=Random([1..n]);
            rand_b:=Random([1..n]);

            #Print(rand_a," ",rand_b,"\n");
        
            if rand_a=rand_b then  gaussian:=TriangulizedMat( perm_G);
            else
                perm:=(rand_a,rand_b);
                perm_mat:=PermutationMat( perm, n, GF(2) );
                perm_G:=perm_G*perm_mat;
                gaussian:=TriangulizedMat( perm_G );
                #Display(gaussian);
                #Print("\n");
            fi;
        
            wt:=WeightVecFFE( gaussian[k]);
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


CSS_dist:=function(Hx,Hz)
        local G,i1,i2;
        G:=[];
        return G;
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
#main function:

   
for i in [1..5] do
    n1:=i*n0;
    n2t:=n1;
#they are the row_numbers for H1, H2^T


     H1:=cyclic_H(h1x,n1);
     H2t:=cyclic_H(h2xt,n1);
     H1t:=TransposedMat(H1);
     H2:=TransposedMat(H2t);
     di_1:=DimensionsMat( H1);
     di_2:=DimensionsMat( H2);
     n1:=di_1[2];
     n2:=di_2[2];
     r1:=di_1[1];
     r2:=di_2[1];

#notice: r1,r2 are not the ranks


    s1:=n1-r1;
    k1:=n1-RankMat(H1);
    s2:=n2-r2;
    k2:=n1-RankMat(H2);
    
    n:=r2*n1+r1*n2;
    k:=2*k1*k2-k1*s2-k2*s1;

    Hx:=HP_Hx(h1x,h2xt,n1);
    Hz:=HP_Hz(h1x,h2xt,n1);
    d_upper_bound:=CSS_dist(Hx,Hz);;
  
    Print("for H1 from h1x=",h1x, "and H2^T from h2xt=",h2xt, "  with x^(",n1,"*",i,")+1 we have:\n","n=",n,", k=",k," ,d<=",d_upper_bound,"\n");


    #Print("for n1=",n1," G1=\n");
    #Display(G1);
  


od;




