(* ::Package:: *)

(*Rui-Chun Xiao, Anhui University, xiaoruichun@ahu.edu.cn*)
(* This module is the program to calculate the response tensors constrained by the spin groups. *)
SpinGroupSHE[SG_] := Module[{index, Index, SHEtensor, i,j,k,l,m,n,r, 
Tensor3NewOdd, Tensor3NewEven, TensorTmp, Us, Rc,mysolve },

(*To define the Spin Hall tensor*)
index={"x","y","z"}; Index={"\[Alpha]","\[Beta]","\[Gamma]"};
SHEtensor=ConstantArray[0,{3,3,3}];
For[i=1,i<=3,i++,
 For[j=1,j<=3,j++,
  For[k=1,k<=3,k++,
       SHEtensor[[i,j,k]]=Overscript[Subscript["\[Sigma]", index[[j]]<>index[[k]]],Index[[i]] ];
 ]]];
(******To find the nonzero elements******)

Tensor3NewOdd=SHEtensor;
For [r=1,r<= Dimensions[SG][[1]],r++,  
  Us=SG[[r,1]]; Rc=SG[[r,2]];
  TensorTmp=ConstantArray[0,{3,3,3}];  
   For[i=1,i<=3,i++,  For[j=1,j<=3,j++, For[k=1,k<=3,k++,
     For[l=1,l<=3,l++, For[m=1,m<=3,m++, For[n=1,n<=3,n++,
            TensorTmp[[i,j,k]]=TensorTmp[[i,j,k]]+Us[[i,l]]*Rc[[j,m]]*Rc[[k,n]]*Tensor3NewOdd[[l,m,n]]; (*Det[Us]*)
         ]]]]]];
   mysolve=Reduce[Tensor3NewOdd==TensorTmp,Flatten[SHEtensor]];
   Tensor3NewOdd=Tensor3NewOdd/. ToRules[mysolve];
 ]; 

(*Output the odd spin Hall matrix*)
Print["Independent odd SHE tensor elements without SOC effect: ",Variables[Tensor3NewOdd],". Number: ",Length[Variables[Tensor3NewOdd]], "."];
Print["\!\(\*SuperscriptBox[\(\[Sigma]\), \(\[Alpha]\)]\)= ",MatrixForm[Tensor3NewOdd[[1]]], ", ",
"\!\(\*SuperscriptBox[\(\[Sigma]\), \(\[Beta]\)]\)= ",MatrixForm[Tensor3NewOdd[[2]]], ", ",
"\!\(\*SuperscriptBox[\(\[Sigma]\), \(\[Gamma]\)]\)= ",MatrixForm[Tensor3NewOdd[[3]]], "."
];

(***********Even part**********************)
Tensor3NewEven=SHEtensor;
For [r=1,r<= Dimensions[SG][[1]],r++,  
Us=SG[[r,1]]; Rc=SG[[r,2]];
  TensorTmp=ConstantArray[0,{3,3,3}];  
  For[i=1,i<=3,i++,   For[j=1,j<=3,j++, For[k=1,k<=3,k++,
     For[l=1,l<=3,l++,  For[m=1,m<=3,m++,  For[n=1,n<=3,n++,
            TensorTmp[[i,j,k]]=TensorTmp[[i,j,k]]+Det[Us]*Us[[i,l]]*Rc[[j,m]]*Rc[[k,n]]*Tensor3NewEven[[l,m,n]]; 
         ]]]]]];
   mysolve=Reduce[Tensor3NewEven==TensorTmp,Flatten[SHEtensor]];
   Tensor3NewEven=Tensor3NewEven/. ToRules[mysolve];
 ]; 

(*Output the even spin Hall matrix*)
Print["Independent even SHE tensor elements without SOC effect: ",Variables[Tensor3NewEven], ". Number: ",Length[Variables[Tensor3NewEven]], "."];
Print["\!\(\*SuperscriptBox[\(\[Sigma]\), \(\[Alpha]\)]\)= ",MatrixForm[Tensor3NewEven[[1]]], ", ",
"\!\(\*SuperscriptBox[\(\[Sigma]\), \(\[Beta]\)]\)= ",MatrixForm[Tensor3NewEven[[2]]], ", ",
"\!\(\*SuperscriptBox[\(\[Sigma]\), \(\[Gamma]\)]\)= ",MatrixForm[Tensor3NewEven[[3]]], "."
];
](*The end*)


SpinGroupSPGE[SG_] := Module[{index, Index, IndexSpin, SPGEtensor, i,j,k,l,m,n,r,a,b,c,d,s,
Tensor4New, TensorTmp, Us, Rc,mysolve, SPGE36Odd\:ff0cSPGE36Even }, 

(*To define the SPGE tensor*)
IndexSpin={"\[Alpha]","\[Beta]","\[Gamma]"}; Index={"X","Y","Z"}; index={"x","y","z"};
SPGEtensor=ConstantArray[0,{3,3,3,3}];
For[i=1,i<=3,i++, For[j=1,j<=3,j++, 
  For[k=1,k<=3,k++, For[l=1,l<=3,l++,
  maxValue=Max[k,l];
  minValue=Min[k,l];
  SPGEtensor[[i,j,k,l]]=Overscript[Subscript["\[Sigma]",Index[[j]]<> index[[minValue]]<>index[[maxValue]]],IndexSpin[[i]]];
]]]];
(*Print["Independent SPGE tensor elements : ",Variables[SPGEtensor],", Number: ",Length[Variables[SPGEtensor]]];*)
(******To find the nonzero odd SPGE elements******)
Tensor4New=SPGEtensor;
For [r=1,r<= Dimensions[SG][[1]],r++,  
  Us=SG[[r,1]]; Rc=SG[[r,2]];
  TensorTmp=ConstantArray[0,{3,3,3,3}];  
  For[i=1,i<=3,i++,  For[j=1,j<=3,j++, For[k=1,k<=3,k++, For[m=1,m<=3,m++,  
       For[a=1,a<=3,a++, For[b=1,b<=3,b++, For[c=1,c<=3,c++, For[d=1,d<=3,d++,
          TensorTmp[[i,j,k,m]]=TensorTmp[[i,j,k,m]]+Us[[i,a]]*Rc[[j,b]]*Rc[[k,c]]*Rc[[m,d]]*Tensor4New[[a,b,c,d]]; 
         ]]]]]]]];
   mysolve=Reduce[Tensor4New==TensorTmp,Flatten[SPGEtensor]];
   Tensor4New=Tensor4New/. ToRules[mysolve];
 ]; 

(*Print the SHG tensro in the 3*6 form*)
SPGE36Odd=ConstantArray[0,{3,3,6}];
alpha={1,2,3,2,1,1};
beta={1,2,3,3,3,2};

For[s=1,s<=3,s++,
  For[i=1,i<=3,i++,
  For[j=1,j<=6,j++,
    k=alpha[[j]]; l=beta[[j]];
    SPGE36Odd[[s,i,j]]=Tensor4New[[s,i,k,l]];
]]];

Print["Independent odd SPGE tensor elements :",Variables[Tensor4New],". Number: ",Length[Variables[Tensor4New]], "."];
Print["\!\(\*SuperscriptBox[\(\[Sigma]\), \(Sx \)]\)=  ",MatrixForm[SPGE36Odd[[1,;;,;;]]]];
Print["\!\(\*SuperscriptBox[\(\[Sigma]\), \(Sy \)]\)= ",MatrixForm[SPGE36Odd[[2,;;,;;]]]];
Print["\!\(\*SuperscriptBox[\(\[Sigma]\), \(Sz \)]\)= ",MatrixForm[SPGE36Odd[[3,;;,;;]]]];


(******To find the nonzero even SPGE elements******)
Tensor4New=SPGEtensor;
For [r=1,r<= Dimensions[SG][[1]],r++,  
Us=SG[[r,1]]; Rc=SG[[r,2]];
  TensorTmp=ConstantArray[0,{3,3,3,3}];  
  For[i=1,i<=3,i++,  For[j=1,j<=3,j++, For[k=1,k<=3,k++, For[m=1,m<=3,m++,  
       For[a=1,a<=3,a++, For[b=1,b<=3,b++, For[c=1,c<=3,c++, For[d=1,d<=3,d++,
          TensorTmp[[i,j,k,m]]=TensorTmp[[i,j,k,m]]+Det[Us]*Us[[i,a]]*Rc[[j,b]]*Rc[[k,c]]*Rc[[m,d]]*Tensor4New[[a,b,c,d]];  (*Det[Us]*)
         ]]]]]]]];
   mysolve=Reduce[Tensor4New==TensorTmp,Flatten[SPGEtensor]];
   Tensor4New=Tensor4New/. ToRules[mysolve];
 ]; 

(*Print the SHG tensro in the 3*6 form*)
SPGE36Even=ConstantArray[0,{3,3,6}];
alpha={1,2,3,2,1,1};
beta={1,2,3,3,3,2};

For[s=1,s<=3,s++,
  For[i=1,i<=3,i++,
  For[j=1,j<=6,j++,
    k=alpha[[j]]; l=beta[[j]];
    SPGE36Even[[s,i,j]]=Tensor4New[[s,i,k,l]];
]]];

Print["Independent even SPGE tensor elements :",Variables[Tensor4New],". Number: ",Length[Variables[Tensor4New]], "."];
Print["\!\(\*SuperscriptBox[\(\[Sigma]\), \(Sx \)]\)=  ",MatrixForm[SPGE36Even[[1,;;,;;]]]];
Print["\!\(\*SuperscriptBox[\(\[Sigma]\), \(Sy \)]\)= ",MatrixForm[SPGE36Even[[2,;;,;;]]]];
Print["\!\(\*SuperscriptBox[\(\[Sigma]\), \(Sz \)]\)= ",MatrixForm[SPGE36Even[[3,;;,;;]]]];
](*The end*)


SpinGroupSHG[SG_] := Module[{index, Index, SHGtensor, i,j,k,l,m,n,r, 
Tensor3, Tensor3Even, Tensor3Odd, TensorTmpEven, TensorTmpOdd, Us, Rc,mysolveEven, mysolveOdd,
SHG36Even, SHG36Odd},

(*To define the SHG matrix*)
 Tensor3 = ConstantArray[0, {3, 3, 3}];  
 For[i=1,i<=3,i++,  For[j=1,j<=3,j++, For[k=1,k<=3,k++,
    maxValue=Max[j,k];
    minValue=Min[j,k];
    Tensor3[[i,j,k]]=Subscript[\[Chi], ToString[i]<> ToString[minValue] <> ToString[maxValue]];
  ]]];
(******To find the nonzero elements******)

Tensor3Even=Tensor3; Tensor3Odd=Tensor3;
For [r=1,r<= Dimensions[SG][[1]],r++,  
Us=SG[[r,1]]; Rc=SG[[r,2]];
  TensorTmpEven=ConstantArray[0,{3,3,3}]; 
  TensorTmpOdd=ConstantArray[0,{3,3,3}];  
  For[i=1,i<=3,i++,  For[j=1,j<=3,j++,  For[k=1,k<=3,k++,
     For[l=1,l<=3,l++, For[m=1,m<=3,m++,  For[n=1,n<=3,n++,
            TensorTmpEven[[i,j,k]]=TensorTmpEven[[i,j,k]]+Rc[[i,l]]*Rc[[j,m]]*Rc[[k,n]]*Tensor3Even[[l,m,n]];
            TensorTmpOdd[[i,j,k]]=TensorTmpOdd[[i,j,k]]+Det[Us]*Rc[[i,l]]*Rc[[j,m]]*Rc[[k,n]]*Tensor3Odd[[l,m,n]];
         ]]]]]];
   mysolveEven=Reduce[Tensor3Even==TensorTmpEven,Flatten[Tensor3]];
   Tensor3Even=Tensor3Even/. ToRules[mysolveEven];
   mysolveOdd=Reduce[Tensor3Odd==TensorTmpOdd,Flatten[Tensor3]];
   Tensor3Odd=Tensor3Odd/. ToRules[mysolveOdd];
 ]; 

SHG36Even=ConstantArray[0,{3,6}]; SHG36Odd=ConstantArray[0,{3,6}];
alpha={1,2,3,2,1,1};
beta={1,2,3,3,3,2};
For[i=1,i<=3,i++,
 For[j=1,j<=6,j++,
   k=alpha[[j]]; l=beta[[j]];
   SHG36Even[[i,j]]=Tensor3Even[[i,k,l]];
   SHG36Odd[[i,j]]=Tensor3Odd[[i,k,l]]
]];

(*Output the SHG matrix*)
Print["Independent even SHG tensor elements without SOC effect: ",Variables[Tensor3Even],". Number: ",Length[Variables[Tensor3Even]], "."];
Print[MatrixForm[SHG36Even]];
Print["Independent odd SHG tensor elements without SOC effect: ",Variables[Tensor3Odd],". Number: ",Length[Variables[Tensor3Odd]], "."];
Print[MatrixForm[SHG36Odd]];
](*The end*)
