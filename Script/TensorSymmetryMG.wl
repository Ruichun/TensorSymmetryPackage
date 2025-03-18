(* ::Package:: *)

(*Rui-Chun Xiao, Anhui University, xiaoruichun@ahu.edu.cn*)
(* This module is the program to calculate the response tensors constrained by the magnetic groups. *)
ConductivityTensor[MPGop_,Basis_] := Module[{TensorSym, TensorAsym, TensorSymTemp,TensorAsymTemp,
TensorSymOrg, TensorAsymOrg, RiCry, RiDescart, solutionSym, solutionAsym, IsUnitary},
  
 (*Define the 2rd tensor*)
 TensorSymOrg=ConstantArray[0, {3, 3}];  TensorAsymOrg=ConstantArray[0, {3, 3}]; 
  TensorSymOrg=({
 {Subscript[\[Sigma], 11], Subscript[\[Sigma], 12], Subscript[\[Sigma], 13]},
 {Subscript[\[Sigma], 12], Subscript[\[Sigma], 22], Subscript[\[Sigma], 23]},
 {Subscript[\[Sigma], 13], Subscript[\[Sigma], 23], Subscript[\[Sigma], 33]}
});  
  TensorAsymOrg=({
 {0, Subscript[A, 12], -Subscript[A, 31]},
 {-Subscript[A, 12], 0, Subscript[A, 23]},
 {Subscript[A, 31], -Subscript[A, 23], 0}
});  
  TensorSym=TensorSymOrg; TensorAsym=TensorAsymOrg;
  
  For[r = 1, r <= Dimensions[MPGop][[1]], r++,  
     RiCry = MPGop[[r,1]];  
     RiDescart=Basis . RiCry . Inverse[Basis];
     If[ MPGop[[r,2]]==0,  IsUnitary=1 ];
     If[ MPGop[[r,2]]==1,  IsUnitary=-1 ];
     
     TensorSymTemp = ConstantArray[0, {3, 3}];  
     TensorAsymTemp = ConstantArray[0, {3, 3}]; 
     
     TensorSymTemp=RiDescart . TensorSym . Inverse[RiDescart];
     TensorAsymTemp=IsUnitary*RiDescart . TensorAsym . Inverse[RiDescart];
     
    solutionSym = Reduce[TensorSymTemp ==TensorSym , Flatten[TensorSymOrg]];
    solutionAsym = Reduce[TensorAsymTemp == TensorAsym, {Subscript[A, 12],Subscript[A, 31],Subscript[A, 23]}];  
      
    TensorSym  = TensorSym  /. ToRules[solutionSym];
    TensorAsym = TensorAsym /. ToRules[solutionAsym];  
    ];  
    
Print["Independent tensor elements of symmetry conductivity tensor: ",Variables[TensorSym],". Number: ",Length[Variables[TensorSym]],"."]; 
Print[MatrixForm[TensorSym]];

Print["Independent tensor elements of antisymmetry conductivity tensor: ",Variables[TensorAsym],". Number: ",Length[Variables[TensorAsym]],"."]; 
Print[MatrixForm[TensorAsym]];

] (*The end*)


SHGtensor[MPGop_,Basis_] := Module[{Tensor3, Tensor3NewEven, Tensor3tmpEven,Tensor3NewOdd, Tensor3tmpOdd,
RiCry, RiDescart, solutionsEven,  solutionsOdd, SHG36Even, SHG36Odd, TensorList},
  
 (*Define the 3rd tensor*)
 Tensor3 = ConstantArray[0, {3, 3, 3}];  
 For[i=1,i<=3,i++,
  For[j=1,j<=3,j++,
   For[k=1,k<=3,k++,
    maxValue=Max[j,k];
    minValue=Min[j,k];
    Tensor3[[i,j,k]]=Subscript[\[Chi], ToString[i]<> ToString[minValue] <> ToString[maxValue]];
   ]]];
  
    Tensor3NewEven = Tensor3;  Tensor3NewOdd = Tensor3;
    For[r = 1, r <= Dimensions[MPGop][[1]], r++,  
        RiCry = MPGop[[r,1]];  
        RiDescart=Basis . RiCry . Inverse[Basis];
        If[ MPGop[[r,2]]==0,   IsUnitary=1 ];
        If[ MPGop[[r,2]]==1,  IsUnitary=-1 ];
            
        Tensor3tmpEven = ConstantArray[0, {3, 3, 3}];  
        Tensor3tmpOdd = ConstantArray[0, {3, 3, 3}]; 
        For[i = 1, i <= 3, i++,  
            For[j = 1, j <= 3, j++,  
                For[k = 1, k <= 3, k++,  
                    For[l = 1, l <= 3, l++,  
                        For[m = 1, m <= 3, m++,  
                            For[n = 1, n <= 3, n++,  
                                Tensor3tmpEven[[i, j, k]] += RiDescart[[i, l]] * RiDescart[[j, m]] * RiDescart[[k, n]] * Tensor3NewEven[[l, m, n]];
                                Tensor3tmpOdd[[i, j, k]] += IsUnitary*RiDescart[[i, l]] * RiDescart[[j, m]] * RiDescart[[k, n]] * Tensor3NewOdd[[l, m, n]];    
                             ] ] ] ] ] ];    
        solutionsEven = Reduce[Flatten[Tensor3tmpEven - Tensor3NewEven] == 0, Flatten[Tensor3]];
        solutionsOdd  = Reduce[Flatten[Tensor3tmpOdd - Tensor3NewOdd] == 0, Flatten[Tensor3]];   
      
        Tensor3NewEven = Tensor3NewEven /.ToRules[solutionsEven];
        Tensor3NewOdd = Tensor3NewOdd /.ToRules[solutionsOdd];  
    ];  
    
SHG36Even=ConstantArray[0,{3,6}]; SHG36Odd=ConstantArray[0,{3,6}];
alpha={1,2,3,2,1,1};
beta={1,2,3,3,3,2};
For[i=1,i<=3,i++,
For[j=1,j<=6,j++,
   k=alpha[[j]]; l=beta[[j]];
   SHG36Even[[i,j]]=Tensor3NewEven[[i,k,l]];
   SHG36Odd[[i,j]]=Tensor3NewOdd[[i,k,l]]
]];

Print["Independent elements of even SHG tensor: ",Variables[Tensor3NewEven],". Number: ",Length[Variables[Tensor3NewEven]],"."]; 
Print[MatrixForm[SHG36Even]];

Print["Independent tensor elements of odd SHG tensor: ",Variables[Tensor3NewOdd],". Number: ",Length[Variables[Tensor3NewOdd]],"."]; 
Print[MatrixForm[SHG36Odd]];

] (*The end*)


SPGEtensor[MPGop_,Basis_] := Module[{SPGEtensor, Tensor4NewOdd, Tensor4NewEven, Tensor4tmpEven, Tensor4tmpOdd,
RiCry, Ri, mysolveEven,  mysolveOdd, SPGE36Even, SPGE36Odd, TensorList},
(*SPGE: Spin Spin photogalvanic effect*)
IndexSpin={"Sx","Sy","Sz"}; Index={"X","Y","Z"}; index={"x","y","z"};
SPGEtensor=ConstantArray[0,{3,3,3,3}];
For[i=1,i<=3,i++,
For[j=1,j<=3,j++,
 For[k=1,k<=3,k++,
 For[l=1,l<=3,l++,
  maxValue=Max[k,l];
  minValue=Min[k,l];
  SPGEtensor[[i,j,k,l]]=Overscript[Subscript["\[Sigma]",Index[[j]]<> index[[minValue]]<>index[[maxValue]]],IndexSpin[[i]]];
]]]];

(******To find the nonzero elements******)
Tensor4NewOdd=SPGEtensor; Tensor4NewEven=SPGEtensor;
For [r=1,r<= Dimensions[MPGop][[1]],r++, 
    RiCry = MPGop[[r,1]];  
    Ri=Basis . RiCry . Inverse[Basis];
    If[ MPGop[[r,2]]==0,   IsUnitary=1 ];
    If[ MPGop[[r,2]]==1,  IsUnitary=-1 ];
        
(* to rotate the matrix*)
Tensor4tmpEven=ConstantArray[0,{3,3,3,3}];
Tensor4tmpOdd=ConstantArray[0,{3,3,3,3}];
For[i=1,i<=3,i++, For[j=1,j<=3,j++, For[k=1,k<=3,k++, For[l=1,l<=3,l++, 
  For[m=1,m<=3,m++,  For[n=1,n<=3,n++, For[p=1,p<=3,p++, For[q=1,q<=3,q++,
     Tensor4tmpEven[[i,j,k,l]]=Tensor4tmpEven[[i,j,k,l]]+Det[Ri]*Ri[[i,m]]*Ri[[j,n]]*Ri[[k,p]]*Ri[[l,q]]*Tensor4NewEven[[m,n,p,q]];
     Tensor4tmpOdd[[i,j,k,l]]=Tensor4tmpOdd[[i,j,k,l]]+IsUnitary*Det[Ri]*Ri[[i,m]]*Ri[[j,n]]*Ri[[k,p]]*Ri[[l,q]]*Tensor4NewOdd[[m,n,p,q]];
]]]]]]]];
  mysolveEven=Reduce[Tensor4NewEven==Tensor4tmpEven,Flatten[SPGEtensor]];
  Tensor4NewEven=Tensor4NewEven/. ToRules[mysolveEven];
  mysolveOdd=Reduce[Tensor4NewOdd==Tensor4tmpOdd,Flatten[SPGEtensor]];
  Tensor4NewOdd=Tensor4NewOdd/. ToRules[mysolveOdd];
];

(*Print the SHG tensor in the 3*6 form*)
SPGE36Even=ConstantArray[0,{3,3,6}];
SPGE36Odd=ConstantArray[0,{3,3,6}];
alpha={1,2,3,2,1,1};
beta={1,2,3,3,3,2};

For[s=1,s<=3,s++,
 For[i=1,i<=3,i++,
  For[j=1,j<=6,j++,
    k=alpha[[j]]; l=beta[[j]];
    SPGE36Even[[s,i,j]]=Tensor4NewEven[[s,i,k,l]];
    SPGE36Odd[[s,i,j]]=Tensor4NewOdd[[s,i,k,l]];
]]];

Print["Independent tensor elements for SPGE even part: ",Variables[Tensor4NewEven],". Number: ",Length[Variables[Tensor4NewEven]],"."];
Print["\!\(\*SuperscriptBox[\(\[Sigma]\), \(Sx \((even)\)\)]\)= ",MatrixForm[SPGE36Even[[1,;;,;;]]]];
Print["\!\(\*SuperscriptBox[\(\[Sigma]\), \(Sy \((even)\)\)]\)= ",MatrixForm[SPGE36Even[[2,;;,;;]]]];
Print["\!\(\*SuperscriptBox[\(\[Sigma]\), \(Sz \((even)\)\)]\)= ",MatrixForm[SPGE36Even[[3,;;,;;]]]];
Print["Independent tensor elements for SPGE odd part: ",Variables[Tensor4NewOdd],". Number: ",Length[Variables[Tensor4NewOdd]],"."];
Print["\!\(\*SuperscriptBox[\(\[Sigma]\), \(Sx \((odd)\)\)]\)= ", MatrixForm[SPGE36Odd[[1,;;,;;]]]];
Print["\!\(\*SuperscriptBox[\(\[Sigma]\), \(Sy \((odd)\)\)]\)= ", MatrixForm[SPGE36Odd[[2,;;,;;]]]];
Print["\!\(\*SuperscriptBox[\(\[Sigma]\), \(Sz \((odd)\)\)]\)= ", MatrixForm[SPGE36Odd[[3,;;,;;]]]];

] (*The end*)


SHEtensor[MPGop_,Basis_] := Module[{SHEtensor, Tensor3NewEven, Tensor3NewOdd, TensorTmpEven, TensorTmpOdd,
RiCry, Ri, mysolveEven, mysolveOdd, SPGE36Even, SPGE36Odd, TensorList},
(*SHE: spin Hall effect*)
index={"x","y","z"}; Index={"X","Y","Z"};
SHEtensor=ConstantArray[0,{3,3,3}];
For[i=1,i<=3,i++,
 For[j=1,j<=3,j++,
  For[k=1,k<=3,k++,
  (*SHEtensor[[i,j,k]]=Overscript[Subscript["\[Sigma]", index[[j]]<>index[[k]]],Index[[i]] ]*)
  (*SHEtensor[[i,j,k]]=Subsuperscript[\[Sigma], index[[j]]<>index[[k]], Index[[i]]];*)
  SHEtensor[[i,j,k]]=Subscript["\[Sigma]", Index[[i]]<>index[[j]]<>index[[k]]];
 ]]];

(******To find the nonzero elements******)
Tensor3NewEven=SHEtensor; Tensor3NewOdd=SHEtensor;
For [r=1,r<= Dimensions[MPGop][[1]],r++,  
    RiCry = MPGop[[r,1]];  
    Ri=Basis . RiCry . Inverse[Basis];
    If[ MPGop[[r,2]]==0,   IsUnitary=1 ];
    If[ MPGop[[r,2]]==1,  IsUnitary=-1 ];

  TensorTmpEven=ConstantArray[0,{3,3,3}];  TensorTmpOdd=ConstantArray[0,{3,3,3}];
  For[i=1,i<=3,i++, For[j=1,j<=3,j++, For[k=1,k<=3,k++,
     For[l=1,l<=3,l++, For[m=1,m<=3,m++, For[n=1,n<=3,n++,
      TensorTmpEven[[i,j,k]]=TensorTmpEven[[i,j,k]]+Det[Ri]*Ri[[i,l]]*Ri[[j,m]]*Ri[[k,n]]*Tensor3NewEven[[l,m,n]];
      TensorTmpOdd[[i,j,k]]=TensorTmpOdd[[i,j,k]]+IsUnitary*Det[Ri]*Ri[[i,l]]*Ri[[j,m]]*Ri[[k,n]]*Tensor3NewOdd[[l,m,n]];
]]]]]];
   mysolveEven=Reduce[Tensor3NewEven==TensorTmpEven,Flatten[SHEtensor]];
   mysolveOdd=Reduce[Tensor3NewOdd==TensorTmpOdd,Flatten[SHEtensor]];   
   Tensor3NewEven=Tensor3NewEven/. ToRules[mysolveEven];
   Tensor3NewOdd=Tensor3NewOdd/. ToRules[mysolveOdd];
 ]; 
 
Print["Independent even SHE tensor elements: ",Variables[Tensor3NewEven],". Number: ",Length[Variables[Tensor3NewEven]],"."];
Print["\!\(\*SuperscriptBox[\(\[Sigma]\), \(X\)]\)= ",MatrixForm[Tensor3NewEven[[1]]], ", ",
"\!\(\*SuperscriptBox[\(\[Sigma]\), \(Y\)]\)= ",MatrixForm[Tensor3NewEven[[2]]], ", ",
"\!\(\*SuperscriptBox[\(\[Sigma]\), \(Z\)]\)= ",MatrixForm[Tensor3NewEven[[3]]]];

Print["Independent odd SHE tensor elements: ",Variables[Tensor3NewOdd],". Number: ",Length[Variables[Tensor3NewOdd]],"."];
Print["\!\(\*SuperscriptBox[\(\[Sigma]\), \(X\)]\)= ",MatrixForm[Tensor3NewOdd[[1]]],", ",
"\!\(\*SuperscriptBox[\(\[Sigma]\), \(Y\)]\)= ",MatrixForm[Tensor3NewOdd[[2]]], ", ",
"\!\(\*SuperscriptBox[\(\[Sigma]\), \(Z\)]\)= ",MatrixForm[Tensor3NewOdd[[3]]]
];

](*The end*)


AxialTensor[MPGop_,Basis_] := Module[{Tensor2, TensorEven, TensorOdd, TensorEvenTemp, TensorOddTemp,
RiCry, RiDescart, solutionEven, solutionOdd},
(*Axial tensor: CPGE, Edelstein effect or magnetoelectric effect tensors*)
 (*Define the 2rd tensor*)
 Tensor2=ConstantArray[0, {3, 3}]; 
 For[i=1,i<=3,i++,
   For[j=1,j<=3,j++,
   Tensor2[[i,j]]=Subscript[\[Beta], ToString[i]<> ToString[j] ]
]];

  TensorEven=Tensor2; TensorOdd=Tensor2;
  
  For[r = 1, r <= Dimensions[MPGop][[1]], r++,  
     RiCry = MPGop[[r,1]];  
     RiDescart=Basis . RiCry . Inverse[Basis];
     If[ MPGop[[r,2]]==0,   IsUnitary=1 ];
     If[ MPGop[[r,2]]==1,  IsUnitary=-1 ];
            
     TensorEvenTemp = ConstantArray[0, {3, 3}];  
     TensorOddTemp = ConstantArray[0, {3, 3}]; 
        
     TensorEvenTemp=Det[RiDescart]*RiDescart . TensorEven . Inverse[RiDescart];
     TensorOddTemp=IsUnitary*Det[RiDescart]*RiDescart . TensorOdd . Inverse[RiDescart];
     
     solutionEven = Reduce[TensorEvenTemp ==TensorEven , Flatten[Tensor2]];
     solutionOdd = Reduce[TensorOddTemp == TensorOdd, Flatten[Tensor2]];  
      
    TensorEven  = TensorEven  /. ToRules[solutionEven];
    TensorOdd = TensorOdd /. ToRules[solutionOdd];  
    ];  
    
Print["Independent tensor elements of even axial tensor: ",Variables[TensorEven],". Number: ",Length[Variables[TensorEven]],"."]; 
Print[MatrixForm[TensorEven]];

Print["Independent tensor elements of odd aixal tensor: ",Variables[TensorOdd],". Number: ",Length[Variables[TensorOdd]],"."]; 
Print[MatrixForm[TensorOdd]];

] (*The end*)
