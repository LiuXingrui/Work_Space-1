#o:=0*Z(2);
#id:=Z(2)^0;


hx:=[1,1,0,1];

# a test for h(x)=1+x^2+x^3 and x^7+1
n1:=7;
r1:=Length(hx);
k1:=n1-r1;

# H2=H1^T

n2:=r1;
k2:=0;
r2:=n1;

n:=(n1-k1)^2+n1^2;
k:=k1^2;
r:=n-k;

zero_vector:=[];
i:=1;
while i<=n do
    zero_vector[i]:=0;
    i:=i+1;
od;


h:=[];
H1:=[];


for i in [1..k1] do
    h[i]:=0;
od;

Append(h,hx);

H1[1]:=ShallowCopy(h);
i2:=1;
temp:=[];

for i1 in [2..r1] do
    for i2 in [1..n1] do
        index:=((i1+i2-2) mod n1)+1;
        temp[i2]:=h[index];
        i2:=i2+1;
    od;
    
    H1[i1]:=ShallowCopy(temp);
   
od;

###########
H2:=TransposedMat( H1);


E1:=IdentityMat(r1);
E2:=IdentityMat(r2);
E1_tilt:=IdentityMat(n1);
E2_tilt:=IdentityMat(n2);

Gx_left:=KroneckerProduct( E2,H1 );
Gx_right:=KroneckerProduct( H2,E1 );
Gz_left:=KroneckerProduct(  H1,E1_tilt);
Gz_right:=KroneckerProduct(E2_tilt,H2);

#Print(DimensionsMat(Gx_left));
#Print(DimensionsMat(Gx_right));
#Print(Length(zero_vector));
#Print(DimensionsMat(Gz_left));
#Print(DimensionsMat(Gz_right));


#get stabilizers
S:=[];
for i in [1..r1*r2] do

    S[i]:=[];
    Append(S[i],Gx_left[i]);
    Append(S[i],Gx_right[i]);
    Append(S[i],zero_vector);

od;


for i in [(r1*r2+1)..(r1*r2+r1*n1)] do

    S[i]:=[];
    Append(S[i],zero_vector);
    Append(S[i],Gz_left[i-r1*r2]);
    Append(S[i],Gz_right[i-r1*r2]);
   
od;

#Print(S[1]{[129]});


d:=99999;

i1:=1;


po:=[];
#enumerate all pauli operators to calculate the distance;
#error, Range: <last> must be an integer less than 2^60 (not a integer (>= 2^60)) in for i1 in [ 0 .. 2 ^ (2 * n) - 1 ] do
#so i have to use while loop

#while i1 <=2^(2*n)-1 do
#test:
while i1<=2^15 do;
    is_NS:=false;
    is_S:=false;
    commute_or_not:=0;

#    break_or_not:=false;

    po[1]:=i1 mod 2;
    temp:=(i1-po[1])/2;
    

    for i2 in [2..2*n] do
    
        po[i2]:=temp mod 2;
        temp:= (temp-po[i2])/2;
        i2:=i2+1;

    od;

# if po is in S?

    for i3  in [1..(r1*r2+r1*n1)] do
   
        if po=S[i3] then
           is_S:=true;
           is_NS:=true;

#          break_or_not:=true;
#          break;


        fi;
    od;

#   if break_or_not=true then
#        i1:=i1+1;
#        break;
#   fi;

    if is_S=false then
        for i4 in [1..r] do

 #           Print(po{[1..n]});
 #           Print(S[i4]{[(n+1)..2*n]});

            commute_or_not:=commute_or_not+((po{[1..n]}*S[i4]{[(n+1)..2*n]}+po{[(n+1)..2*n]}*S[i4]{[1..n]}) mod 2);
        od;
    fi;


    if commute_or_not=0 then
        is_NS:=true;
    fi;

    if is_NS and (not is_S) then
        wt:=0;
        for i5 in [1..n] do
            if (not (po[i5]+po[n+i5]=0) ) then
                wt:=wt+1;
            fi;

        od;
        Print(wt,"  ",i1,"\n");
        
        d:=Minimum(d,wt);
    fi;
    i1:=i1+1;
   

od;


Print("for h(x)=1+x^2+x^3 and x^7+1 d=",d);


    

