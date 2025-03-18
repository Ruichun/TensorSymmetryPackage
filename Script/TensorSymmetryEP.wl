(* ::Package:: *)

(*Rui-Chun Xiao, Anhui University, xiaoruichun@ahu.edu.cn*)
(*This module is to establish the relationships between response tenors and the magnetic order (Neel vector).*)
AHEdirection[SymmOperation_,Basis_, Magup_] := Module[{OsMatCry, OsVectCry, isSwitch, BasisT,
OsMat, Index, index, Tensor2, Tensor2New, mysolveRi, Tensor4, Tensor4New, Tensor4tmp, AA, Tensor2rot, Ri,
LengthOfOsMat, m,p,q,l,i,j,k, \[Sigma]1,\[Sigma]3, Neel},

{Nsym,n}=Dimensions[SymmOperation];

OsMatCry=ConstantArray[0,{Nsym,3,3}];
OsVectCry=ConstantArray[0,{Nsym,3}];
(**To construct the R and t**)
For [i=1,i<=Nsym,i++,
OsMatCry[[i,;;,;;]]=SymmOperation[[i,2]];
OsVectCry[[i,;;]]=SymmOperation[[i,3]];
];

(**To determine whether the aom postion is change or not**)
{Nmagup,n}=Dimensions[Magup];
isSwitch=ConstantArray[-1,Nsym];
For[i=1,i<=Nsym,i++,
 k=1; (*Check the first atom is OK*)
 newposition=OsMatCry[[i,;;,;;]] . Magup[[k]]+OsVectCry[[i,;;]];
  For[j=1,j<=Nmagup,j++,
     If[ Mod[newposition[[1]]-Magup[[j,1]],1]==0 && Mod[newposition[[2]]-Magup[[j,2]],1]==0 && Mod[newposition[[3]]-Magup[[j,3]],1]==0,
     isSwitch[[i]]=1;
]]];

OsMat=ConstantArray[0,{Nsym,3,3}];
BasisT=Transpose[Basis];
For[i=1,i<=Nsym,i++,
   OsMat[[i,;;,;;]]=BasisT . OsMatCry[[i,;;,;;]] . Inverse[BasisT];
];
LengthOfOsMat = Length[OsMat];  
(*Define the 1st order matrix*)
 Index={"X","Y","Z"}; index={"x","y","z"};
Tensor2=ConstantArray[0,{3,3}];
For[i=1,i<=3,i++,
 For[j=1,j<=3,j++,
   Tensor2[[i,j]]=Subscript[T, Index[[i]] <> index[[j]]];
]];
(*Solve the 1st order matrix*)
Tensor2New=Tensor2;
For [r=1,r<=LengthOfOsMat,r++,
 Ri=OsMat[[r]];
 Tensor2rot=isSwitch[[r]]*Ri . Tensor2New . Inverse[Ri];
 mysolveRi=Reduce[Tensor2rot==Tensor2New,Flatten[Tensor2]];
 Tensor2New=Tensor2New/.ToRules[mysolveRi];
];

(*******************************)
(*Define the 3nd order term tensor*)
Tensor4=ConstantArray[0,{3,3,3,3}];
For[i=1,i<=3,i++,
 For[j=1,j<=3,j++, For[k=1,k<=3,k++, For[m=1,m<=3,m++,
    AA=Sort[{j,k,m}];
    Tensor4[[i,j,k,m]]=Subscript[T, Index[[i]]<> index[[AA[[1]]]]<> index[[AA[[2]]]]<> index[[AA[[3]]]]];
]]]];
Tensor4New=Tensor4;
For [r=1,r<=LengthOfOsMat,r++,  
Ri=OsMat[[r]];
(* to rotate the matrix*)
Tensor4tmp=ConstantArray[0,{3,3,3,3}];
 For[i=1,i<=3,i++,
  For[j=1,j<=3,j++,
   For[k=1,k<=3,k++,
    For[l=1,l<=3,l++,
     For[m=1,m<=3,m++,
      For[n=1,n<=3,n++,
       For[p=1,p<=3,p++,
        For[q=1,q<=3,q++,
Tensor4tmp[[i,j,k,l]]=Tensor4tmp[[i,j,k,l]]+isSwitch[[r]]*Ri[[i,m]]*Ri[[j,n]]*Ri[[k,p]]*Ri[[l,q]]*Tensor4New[[ m,n,p,q ]];
]]]]]]]];
mysolveRi=Reduce[Tensor4New==Tensor4tmp,Flatten[Tensor4]];
Tensor4New=Tensor4New/.ToRules[mysolveRi];
];

(***print the AHE conductivity vector***)
Clear [i,j,k,l,m,n];
\[Sigma]1={0, 0, 0}; \[Sigma]3={0, 0, 0}; 
Neel={Subscript[n, x],Subscript[n, y],Subscript[n, z]};
(***The 1st order Taylor term**)
For[i=1,i<=3,i++,
  For[j=1,j<=3,j++,
   \[Sigma]1[[i]]=\[Sigma]1[[i]]+Tensor2New[[i,j]]*Neel[[j]];
]];
(***The 3rd order Taylor term**)
For[i=1,i<=3,i++,
 For[j=1,j<=3,j++,
  For[k=1,k<=3,k++,
   For[m=1,m<=3,m++,
    \[Sigma]3[[i]]=\[Sigma]3[[i]]+Tensor4New[[i,j,k,m]]*Neel[[j]]*Neel[[k]]*Neel[[m]];
]]]];

Print["Independent 1st order terms of AHE vector: ",Variables[Tensor2New],", number: ",Length[Variables[Tensor2New]]];
Print["Independent 3rd order terms of AHE vector: ",Variables[Tensor4New],", number: ",Length[Variables[Tensor4New]]];
Print["\!\(\*SubscriptBox[\(\[Sigma]\), \(H\)]\)=", \[Sigma]1//MatrixForm,"+",\[Sigma]3//MatrixForm];
sigma=\[Sigma]1+\[Sigma]3;
sigma
](*The end*)


SHEdirection[SymmOperation_,Basis_, Magup_] := Module[{OsMatCry, OsVectCry, isSwitch, BasisT,
OsMat, Index, index, mysolveRi, Tensor4, Tensor4New, Tensor4tmp, \[Sigma]1,\[Sigma]3, Neel, sigma,
a, b, c, d, e, f, m, N, p, q, l, i, j, k, r, Nsym, NewTensorCache, iS,
Tensor6, Tensor6New,Tensor6tmp,LengthOfOsMat},

{Nsym,n}=Dimensions[SymmOperation];

OsMatCry=ConstantArray[0,{Nsym,3,3}];
OsVectCry=ConstantArray[0,{Nsym,3}];
(**To construct the R and t**)
For [i=1,i<=Nsym,i++,
OsMatCry[[i,;;,;;]]=SymmOperation[[i,2]];
OsVectCry[[i,;;]]=SymmOperation[[i,3]];
];

(**To determine whether the aom postion is change or not**)
{Nmagup,n}=Dimensions[Magup];
isSwitch=ConstantArray[-1,Nsym];
For[i=1,i<=Nsym,i++,
 k=1; (*Chech the First atom is OK*)
 newposition=OsMatCry[[i,;;,;;]] . Magup[[k]]+OsVectCry[[i,;;]];
  For[j=1,j<=Nmagup,j++,
     If[ Mod[newposition[[1]]-Magup[[j,1]],1]==0 && Mod[newposition[[2]]-Magup[[j,2]],1]==0&&Mod[newposition[[3]]-Magup[[j,3]],1]==0,
     isSwitch[[i]]=1;
]]];
(*Print[isSwitch];*)
OsMat=ConstantArray[0,{Nsym,3,3}];
BasisT=Transpose[Basis];
For[i=1,i<=Nsym,i++,
   OsMat[[i,;;,;;]]=BasisT . OsMatCry[[i,;;,;;]] . Inverse[BasisT]
];
LengthOfOsMat = Length[OsMat];  

(**********First Taylor expand*******)
Index={"X","Y","Z"}; index={"x","y","z"};
(*Define the First-order term tensor*)
Tensor4=ConstantArray[0,{3,3,3,3}];
For[i=1,i<=3,i++,
 For[j=1,j<=3,j++,
  For[k=1,k<=3,k++,
    For[m=1,m<=3,m++,
    Tensor4[[i,j,k,m]]=Subscript[T, Index[[i]]<> index[[j]]<> index[[k]]<> index[[m]] ];
]]]];
Tensor4New=Tensor4;
(*******************************)
For [r=1,r<=LengthOfOsMat,r++,  
Ri=OsMat[[r]];
  RotInv=Transpose[Ri];
(* to rotate the matrix*)
Tensor4tmp=ConstantArray[0,{3,3,3,3}];
 For[i=1,i<=3,i++,
  For[j=1,j<=3,j++,
   For[k=1,k<=3,k++,
    For[l=1,l<=3,l++,
     For[m=1,m<=3,m++,
      For[n=1,n<=3,n++,
       For[p=1,p<=3,p++,
        For[q=1,q<=3,q++,
   Tensor4tmp[[i,j,k,l]]=Tensor4tmp[[i,j,k,l]]+isSwitch[[r]]*Ri[[i,m]]*Ri[[j,n]]*Ri[[k,p]]*Ri[[l,q]]*Tensor4New[[m,n,p,q]];
]]]]]]]];
mysolveRi=Reduce[Tensor4New==Tensor4tmp,Flatten[Tensor4]];
Tensor4New=Tensor4New/.ToRules[mysolveRi];
];
(**********End first order Taylor expansion on Neel vector *******)

(***************Third order Taylor expansion on Neel vector*************)
(*Define the 6th order term tensor*)
Tensor6=ConstantArray[0,{3,3,3,3,3,3}]; Tensor6New=ConstantArray[0,{3,3,3,3,3,3}];
For[i=1,i<=3,i++, For[j=1,j<=3,j++, For[k=1,k<=3,k++,
  For[l=1,l<=3,l++,  For[m=1,m<=3,m++, For[n=1,n<=3,n++,
    AA=Sort[{l,m,n}];
    Tensor6[[i,j,k,l,m,n]]=Subscript[T, Index[[i]]<> index[[j]]<> index[[k]]<> index[[AA[[1]]]]
    <> index[[AA[[2]]]]<> index[[AA[[3]]]]];
]]]]]];
Tensor6New=Tensor6;
(****Solve the tensor******)
{time, result} = AbsoluteTiming[
For [r=1,r<=LengthOfOsMat,r++,  
 Print["Calculating the ", r, "-th opreation on SHE, ", "total opreations: ", LengthOfOsMat];
 Ri=OsMat[[r]]; iS=isSwitch[[r]];
 Tensor6tmp=ConstantArray[0,{3,3,3,3,3,3}];
 (********************)
NewTensorCache[i_, j_, k_, l_, m_, n_]:=
  Sum[Sum[Sum[ Sum[Sum[Sum[ 
     iS*Ri[[i,a]]*Ri[[j,b]]*Ri[[k,c]]*Ri[[l,d]]*Ri[[m,e]]*Ri[[n,f]]*Tensor6New[[a,b,c,d,e,f]], 
     {a,1,3}], {b,1,3}], {c,1,3}], {d,1,3}], {e,1,3}], {f,1,3}];
 Tensor6tmp=Table[ NewTensorCache[i,j,k,l,m,n], {i,1,3}, {j,1,3}, {k,1,3}, {l,1,3}, {m,1,3}, {n,1,3}];
 (************************)
 mysolveRi=Reduce[Tensor6New==Tensor6tmp,Flatten[Tensor6]];
 Tensor6New=Tensor6New/.ToRules[mysolveRi];
 ];
];
Print["Total time of 3rd Taylor expansion of SHE: ", time, " seconds."];
(***************End 3rd order Taylor expansion on Neel vector*************)

(***print the \[Sigma] vector***)
Clear [i,j,k,l,m,a,b,c,d,e,f,n];
\[Sigma]1=ConstantArray[0,{3,3,3}]; \[Sigma]3=ConstantArray[0,{3,3,3}]; 
Neel={Subscript[n, x],Subscript[n, y],Subscript[n, z]};
(***The 1st order term**)
For[i=1,i<=3,i++,  For[j=1,j<=3,j++,  
  For[k=1,k<=3,k++, For[m=1,m<=3,m++,
   \[Sigma]1[[i,j,k]]=\[Sigma]1[[i,j,k]]+Tensor4New[[i,j,k,m]]*Neel[[m]];
]]]];

For[i=1,i<=3,i++,  For[j=1,j<=3,j++,  For[k=1,k<=3,k++, 
   For[l=1,l<=3,l++, For[m=1,m<=3,m++, For[N=1,N<=3,N++,
   \[Sigma]3[[i,j,k]]=\[Sigma]3[[i,j,k]]+Tensor6New[[i,j,k,l,m,N]]*Neel[[l]]*Neel[[m]]*Neel[[N]];
]]]]]];

Print["Independent 1st order Taylor elements of SHE tensor: ",Variables[Tensor4New],", number: ",Length[Variables[Tensor4New]]];
Print["\!\(\*SuperscriptBox[\(\[Sigma]\), \(X \((1)\)\)]\)= ",MatrixForm[\[Sigma]1[[1]]],", ",
"\!\(\*SuperscriptBox[\(\[Sigma]\), \(Y \((1)\)\)]\)= ",MatrixForm[\[Sigma]1[[2]]], ", ",
"\!\(\*SuperscriptBox[\(\[Sigma]\), \(Z \((1)\)\)]\)= ",MatrixForm[\[Sigma]1[[3]]], ". "
];

Print["Independent 3rd order Taylor elements of SHE tensor: ",Variables[Tensor6New],", number: ",Length[Variables[Tensor6New]]];
Print["\!\(\*SuperscriptBox[\(\[Sigma]\), \(X \((3)\)\)]\)= ",MatrixForm[\[Sigma]3[[1]]], ", "]; 
Print["\!\(\*SuperscriptBox[\(\[Sigma]\), \(Y \((3)\)\)]\)= ",MatrixForm[\[Sigma]3[[2]]], ", "]; 
Print["\!\(\*SuperscriptBox[\(\[Sigma]\), \(Z \((3)\)\)]\)= ",MatrixForm[\[Sigma]3[[3]]], ". " ]; 
sigma=\[Sigma]1+\[Sigma]3;
](*The end*)


SHGdirection[SymmOperation_,Basis_, Magup_] := Module[{OsMatCry, OsVectCry, isSwitch, BasisT,
OsMat, Index, index, mysolveRi, Tensor4, Tensor4New, Tensor4tmp, \[Chi]1,\[Chi]3, Neel, sigma,
a, b, c, d, e, f, m, N, p, q, l, i, j, k, r, Nsym, AA, BB, SHG1st, SHG3rd, NewTensorCache, iS,
Tensor6, Tensor6New,Tensor6tmp,LengthOfOsMat},

{Nsym,n}=Dimensions[SymmOperation];

OsMatCry=ConstantArray[0,{Nsym,3,3}];
OsVectCry=ConstantArray[0,{Nsym,3}];
(**To construct the R and t**)
For [i=1,i<=Nsym,i++,
OsMatCry[[i,;;,;;]]=SymmOperation[[i,2]];
OsVectCry[[i,;;]]=SymmOperation[[i,3]];
];

(**To determine whether the aom postion is change or not**)
{Nmagup,n}=Dimensions[Magup];
isSwitch=ConstantArray[-1,Nsym];
For[i=1,i<=Nsym,i++,
 k=1; (*Chech the First atom is OK*)
 newposition=OsMatCry[[i,;;,;;]] . Magup[[k]]+OsVectCry[[i,;;]];
  For[j=1,j<=Nmagup,j++,
     If[ Mod[newposition[[1]]-Magup[[j,1]],1]==0 && Mod[newposition[[2]]-Magup[[j,2]],1]==0&&Mod[newposition[[3]]-Magup[[j,3]],1]==0,
     isSwitch[[i]]=1;
]]];

OsMat=ConstantArray[0,{Nsym,3,3}];
BasisT=Transpose[Basis];
For[i=1,i<=Nsym,i++,
   OsMat[[i,;;,;;]]=BasisT . OsMatCry[[i,;;,;;]] . Inverse[BasisT]
];
LengthOfOsMat = Length[OsMat];  

(*******************************)
Index={"X","Y","Z"}; index={"x","y","z"};
(*Define the First-order term tensor*)
Tensor4=ConstantArray[0,{3,3,3,3}];
For[i=1,i<=3,i++,
 For[j=1,j<=3,j++,
  For[k=1,k<=3,k++,
      AA=Sort[{j,k}];
    For[m=1,m<=3,m++,
    Tensor4[[i,j,k,m]]=Subscript[T, Index[[i]]<> index[[AA[[1]]]]<> index[[AA[[2]]]]<> index[[m]] ];
]]]];
Tensor4New=Tensor4;
(*******************************)
For [r=1,r<=LengthOfOsMat,r++,  
Ri=OsMat[[r]];
  RotInv=Transpose[Ri];
(* to rotate the matrix*)
Tensor4tmp=ConstantArray[0,{3,3,3,3}];
 For[i=1,i<=3,i++,
  For[j=1,j<=3,j++,
   For[k=1,k<=3,k++,
    For[l=1,l<=3,l++,
     For[m=1,m<=3,m++,
      For[n=1,n<=3,n++,
       For[p=1,p<=3,p++,
        For[q=1,q<=3,q++,
Tensor4tmp[[i,j,k,l]]=Tensor4tmp[[i,j,k,l]]+isSwitch[[r]]*Det[Ri]*Ri[[i,m]]*Ri[[j,n]]*Ri[[k,p]]*Ri[[l,q]]*Tensor4New[[m,n,p,q]];
]]]]]]]];
mysolveRi=Reduce[Tensor4New==Tensor4tmp,Flatten[Tensor4]];
Tensor4New=Tensor4New/.ToRules[mysolveRi];
];

(***************3rd order Taylor explaned*************)
(*Define the 6th order tensor*)
Tensor6=ConstantArray[0,{3,3,3,3,3,3}]; Tensor6New=ConstantArray[0,{3,3,3,3,3,3}];
For[i=1,i<=3,i++, For[j=1,j<=3,j++, For[k=1,k<=3,k++,
   For[l=1,l<=3,l++,  For[m=1,m<=3,m++, For[n=1,n<=3,n++,
    AA=Sort[{j,k}];  BB=Sort[{l,m,n}];
    Tensor6[[i,j,k,l,m,n]]=Subscript[T, Index[[i]]<> index[[AA[[1]]]]<> index[[AA[[2]]]]<> index[[BB[[1]]]]
    <> index[[BB[[2]]]]<> index[[BB[[3]]]]];
]]] ]]];

Tensor6New=Tensor6;
(****Solve the Tensor******)
{time, result} = AbsoluteTiming[
For [r=1,r<=LengthOfOsMat,r++,  
 Print["Calculating the ", r, "-th opreation on SHG, ", "total opreations: ", LengthOfOsMat];
 Ri=OsMat[[r]]; iS=isSwitch[[r]]; DetR=Det[Ri];
 Tensor6tmp=ConstantArray[0,{3,3,3,3,3,3}];
 (********************)
NewTensorCache[i_, j_, k_, l_, m_, n_]:=
  Sum[Sum[Sum[ Sum[Sum[Sum[ 
     iS*DetR*Ri[[i,a]]*Ri[[j,b]]*Ri[[k,c]]*Ri[[l,d]]*Ri[[m,e]]*Ri[[n,f]]*Tensor6New[[a,b,c,d,e,f]], 
     {a,1,3}], {b,1,3}], {c,1,3}], {d,1,3}], {e,1,3}], {f,1,3}];
 Tensor6tmp=Table[ NewTensorCache[i,j,k,l,m,n], {i,1,3}, {j,1,3}, {k,1,3}, {l,1,3}, {m,1,3}, {n,1,3}];
 (************************)
 mysolveRi=Reduce[Tensor6New==Tensor6tmp,Flatten[Tensor6]];
 Tensor6New=Tensor6New/.ToRules[mysolveRi];
 ];
];
Print["Total time: ", time, " seconds."];
(***************End 3rd order Taylor explaned*************)

(***print the \[Sigma] vector***)
Clear [i,j,k,l,m,n,a,b,c,d,N];
\[Chi]1=ConstantArray[0,{3,3,3}]; \[Chi]3=ConstantArray[0,{3,3,3}]; 
Neel={Subscript[n, x],Subscript[n, y],Subscript[n, z]};
(***The 1st order term**)
For[i=1,i<=3,i++,  For[j=1,j<=3,j++,  
  For[k=1,k<=3,k++, For[m=1,m<=3,m++,
   \[Chi]1[[i,j,k]]=\[Chi]1[[i,j,k]]+Tensor4New[[i,j,k,m]]*Neel[[m]];
]]]];

SHG1st=ConstantArray[0,{3,6}]; 
alpha={1,2,3,2,1,1};
beta={1,2,3,3,3,2};
For[i=1,i<=3,i++, For[j=1,j<=6,j++,
   k=alpha[[j]]; l=beta[[j]];
   SHG1st[[i,j]]=\[Chi]1[[i,k,l]];
]];
Print["Independent 1st order terms of Neel vector on SHG: ",Variables[Tensor4New],", number: ",Length[Variables[Tensor4New]]];
Print["SHG1st =",MatrixForm[SHG1st]];

(************3rd**************)
For[i=1,i<=3,i++,  For[j=1,j<=3,j++,  For[k=1,k<=3,k++, 
   For[l=1,l<=3,l++, For[m=1,m<=3,m++, For[N=1,N<=3,N++,
   \[Chi]3[[i,j,k]]=\[Chi]3[[i,j,k]]+Tensor6New[[i,j,k,l,m,N]]*Neel[[l]]*Neel[[m]]*Neel[[N]];
]]]]]];
SHG3rd=ConstantArray[0,{3,6}]; 
For[i=1,i<=3,i++, For[j=1,j<=6,j++,
   k=alpha[[j]]; l=beta[[j]];
   SHG3rd[[i,j]]=\[Chi]3[[i,k,l]];
]];
Print["Independent 3rd order terms of Neel vector on SHG: ",Variables[Tensor6New],", number: ",Length[Variables[Tensor6New]]];
Print["SHG3rd =",MatrixForm[SHG3rd]];
](*The end*)
