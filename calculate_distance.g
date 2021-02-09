x:=Indeterminate(GF(2),"x");

#hx:=1+x+x^2+x^4;
#n0:=7;
#very strange results: d=3*i, for x^(n0*i)+1?

hx:=1+x+x^2+x^3+x^4;
n0:=10;

hx_vector:=CoefficientsOfUnivariatePolynomial( hx);
k1:=Length(hx_vector)-1;

#functions:
############################################

#construct gx=(x^n+1)/hx
        
#construct generator matrix G1 for H1.

construct_G_from_hx:=function(hx,n)
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

      
        return G;
end;


rand_dist:=function(G,k,n)
        local perm_num,i,perm,perm_mat,perm_G,gaussian,d,wt,rand_a,rand_b;
        d:=n;
        perm_num:=10;
        for i in [1..perm_num] do
            rand_a:=Random([1..n]);
            rand_b:=Random([1..n]);

            #Print(rand_a," ",rand_b,"\n");
        
            if rand_a=rand_b then  gaussian:=TriangulizedMat( G);
            else
                perm:=(rand_a,rand_b);
                perm_mat:=PermutationMat( perm, n, GF(2) );
                perm_G:=G*perm_mat;
                gaussian:=TriangulizedMat( perm_G );
                #Display(gaussian);
                #Print("\n");
            fi;
        
            wt:=WeightVecFFE( gaussian[k]);
            if wt<d then d:=wt;  fi;
        od;
        return d;
end;
        
########################################
#main function:

   
for i in [1..100] do
        
#construct G1
        
    n1:=n0*i;
    n:=(n1-k1)^2+n1^2;
    k:=k1^2;
    
    G1:=construct_G_from_hx(hx,n1);
    d_upper_bound:=rand_dist(G1,k1,n1);

    Print("for hx=",hx, " with x^(",n0,"*",i,")+1 we have:\n","n=",n,", k=",k," ,d<=",d_upper_bound,"\n");


    #Print("for n1=",n1," G1=\n");
    #Display(G1);
    Print("\n");

#for H1 cyclic, H2=H1^T, d=d1

od;




