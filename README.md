# TensorSymmetryPackage

See Ref. [[**TensorSymmetry: a package to get the symmetry-adapted tensors to disentangle the spin-orbit coupling effect and establish the relationship with magnetic order**](https://arxiv.org/abs/2411.10147) for more information.  

The software package consists of four kinds of modules:

```
(1)  SGData.wl, MPGop.wl: These modules primarily contain all the symmetry operations for space groups and magnetic point groups.
(2) TensorSymmetryMG.wl: This module is the program for calculating the response tensors constrained by the magnetic groups.
(3) TensorSymmetrySG.wl: This module is to study the response tensors under the spin groups.
(4) TensorSymmetryEP.wl: This module is to establish the relationships between response tensors and the magnetic orders.
```

Here, we introduce the usage of the TensorSymmetry briefly. 

##A.  *Response tensors constrained by magnetic groups.* ##

(1)   Open the user interface file *ToUseMagneticGroup.nb*. This file should be put in the same directory as the script folder.

(2)   Enter the magnetic point group number, and use the “MPGmassage” function to get the symmetry operation matrix (MPGop) and basis vectors (Basis). One can use FINDSYM (https://stokes.byu.edu/iso/findsym.php) or a related website to determine the magnetic point group and corresponding magnetic point group of the magnetic material.

```
MPGno={14,3};  

{MPGop,Basis}=MPGmassage[MPGno];  
```

(3)  Input the desired tensor function and press enter to calculate it.  

```
ConductivityTensor[MPGop,  Basis]; (*Calculate the normal  conductivity and anomalous Hall tensor*)  
SHGtensor[MPGop,  Basis]; (*Calculate the second  harmonic generation tensor*)   
SPGEtensor[MPGop,  Basis]; (*Calculate the spin  photovoltaic effect tensor*)  SHEtensor[MPGop,  Basis];  (*Calculate the spin Hall  tensor*)  
AxialTensor[MPGop,  Basis];  (*Calculate the axial tensor,  such as CPGE, Edelstein effect, magnetoelectric effect*)  
```

We have also provided the table dictionary on dielectric tensors, AHE tensors, SHG tensors, SHE tensors, and SPGE tensors for all magnetic point groups. These can be found on the website https://ruichun.github.io/TensorSymmetry/.

##*B.*  *Response tensors constrained by the spin group.*

​    Since there is no comprehensive operation database on spin groups, we use third-party tools like FINDSPINGROUP (https://findspingroup.com/)  to determine the spin group operations for specific magnetically ordered materials. Then, the user needs to manually input the generator operation elements   $\{U\|R\}$. All the possible matrices of *U* and *R* are already given in our program. 

​    Then, with these matrices, we construct the generator operations for a specific spin group. For example  the "ToUseSpinGroup.nb":  

```
SGop1={S2,C4z}; (*The first of generator operation of  spin group operation in {U||R} form*)  
SGop2={C1,mz}; (*The second of generator operation of spin  group operation*)  
SGop3={S2,mx};  
SGop4={C1,mxy};  
SGop5={C∞,C1};  
SGop6={mx,C1};  
RuO2SG={SGop1, SGop2, SGop3, SGop4, SGop5, SGop6};  
SpinGroupSHE[RuO2SG]    
SpinGroupSHG[RuO2SG]    
SpinGroupSPGE[RuO2SG]    
SpinGroupSHE[RuO2SG]   
```

##*C.*  *Response tensor with the magnetic order. 

We should input its space group number and the positions of the magnetic atoms with spin-up and spin-down states. For example the "ToUseExtrinsicParameters.nb":

```
{spgop, Basis}=SpaceGroupData[136]; (* 136 is the space group number of RuO2*)  
Magup={{0,0,0}}; (*Magnetic atom position with spin up  state, in crystal coordinate*)  
Magdn={{1/2,1/2,1/2}};   (*Magnetic atom position with spin  down state, in crystal coordinate*)  
AHEdirection[spgop, Basis, Magup];  (*AHE tensor with the magnetic order *) 
SHEdirection[spgop, Basis, Magup];  (*SHE tensor with the magnetic order *) 
SHGdirection[spgop, Basis, Magup];  (* SHG tensor with the  magnetic order *)  
```

​       Due to the large dimensionality of the response tensors, the calculation will require several tens of seconds to complete. Additionally, our program analysis has the capability to analyze response tensors with magmatic moment order for collinear ferromagnets, and users simply need to configure all magnetic atoms in the "Magup" matrix while leaving “Magdn” as an empty set. 