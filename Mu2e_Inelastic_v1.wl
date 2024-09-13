(* ::Package:: *)

minus=-1;

JNorm[j_]:=2j+1;

QNorm[j_]:=Sqrt[2 j+1];


Triangular[x_,y_,z_]:=Abs[x-y]<=z<=x+y


RangeM[j_,m_]:=-Abs[j]<=m<=Abs[j]


SixJConditionTriad[j1_,j2_,j3_]:=IntegerQ[j1+j2+j3]&&Triangular[j1,j2,j3]

SixJCondition[j1_,j2_,j3_,J1_,J2_,J3_]:=SixJConditionTriad[j1,j2,j3]&&SixJConditionTriad[j1,J2,J3]&&SixJConditionTriad[J1,j2,J3]&&SixJConditionTriad[J1,J2,j3]


SixJ[{j1_,j2_,j3_},{J1_,J2_,J3_}]:=Which[SixJCondition[j1,j2,j3,J1,J2,J3],SixJSymbol[{j1,j2,j3},{J1,J2,J3}],True,0]


ThreeJCondition[{j1_,m1_},{j2_,m2_},{j3_,m3_}]:=Triangular[j1,j2,j3]&&RangeM[j1,m1]&&RangeM[j2,m2]&&RangeM[j3,m3]&&m1+m2+m3==0

ThreeJ[{j1_,m1_},{j2_,m2_},{j3_,m3_}]:=Which[ThreeJCondition[{j1,m1},{j2,m2},{j3,m3}],ThreeJSymbol[{j1,m1},{j2,m2},{j3,m3}],True,0]


cutsum=50;

SummandNineJ[j1_,j2_,j12_,j3_,j4_,j34_,j13_,j24_,j_,g_]:=minus^(2*g) JNorm[g] SixJ[{j1,j2,j12},{j34,j,g}]*SixJ[{j3,j4,j34},{j2,g,j24}]*SixJ[{j13,j24,j},{g,j1,j3}]


NineJSymbol[{j1_,j2_,j12_},{j3_,j4_,j34_},{j13_,j24_,j_}]:=Sum[SummandNineJ[j1,j2,j12,j3,j4,j34,j13,j24,j,g],{g,0,cutsum+1/2,1/2}]


BesselFactor1[y_,{np_,lp_},{n_,l_},lcap_]:=(2^lcap)/(JNorm[lcap]!!) y^(lcap/2) Exp[-y] Sqrt[(np-1)! (n-1)!]

BesselFactor2[y_,{np_,lp_},{n_,l_},lcap_]:=Sqrt[Gamma[np+lp+1/2]*Gamma[n+l+1/2]]

Summand1[y_,{np_,lp_,mp_},{n_,l_,m_},lcap_]:=minus^(m+mp)/(m! mp! (n-1-m)! (np-1-mp)!)

Summand2[y_,{np_,lp_,mp_},{n_,l_,m_},lcap_]:=Gamma[(l+lp+lcap+2*m+2*mp+3)/2]/(Gamma[l+m+3/2]*Gamma[lp+mp+3/2])

Summand3[y_,{np_,lp_,mp_},{n_,l_,m_},lcap_]:=Hypergeometric1F1[(lcap-lp-l-2*mp-2*m)/2,lcap+3/2,y]

BesselFactor3[y_,{np_,lp_},{n_,l_},lcap_]:=Sum[(Summand1[y,{np,lp,mp},{n,l,m},lcap]*Summand2[y,{np,lp,mp},{n,l,m},lcap]*Summand3[y,{np,lp,mp},{n,l,m},lcap]),{m,0,n-1},{mp,0,np-1}]


BesselElement[y_,{np_,lp_},{n_,l_},lcap_]:=BesselFactor1[y,{np,lp},{n,l},lcap]*BesselFactor2[y,{np,lp},{n,l},lcap]*BesselFactor3[y,{np,lp},{n,l},lcap]


BesselFactor1A[y_,{np_,lp_},{n_,l_},lcap_]:=(2^(lcap-1))/(JNorm[lcap]!!) y^((lcap-1)/2) Exp[-y] Sqrt[(np-1)! (n-1)!]

BesselFactor2[y_,{np_,lp_},{n_,l_},lcap_]:=Sqrt[Gamma[np+lp+1/2]*Gamma[n+l+1/2]]

(**)


Summand1[y_,{np_,lp_,mp_},{n_,l_,m_},lcap_]:=(-1)^(m+mp)/(m! mp! (n-1-m)! (np-1-mp)!)

Summand2A[y_,{np_,lp_,mp_},{n_,l_,m_},lcap_]:=Gamma[(l+lp+lcap+2*m+2*mp+2)/2]/(Gamma[l+m+3/2]*Gamma[lp+mp+3/2])

Summand3A[y_,{np_,lp_,mp_},{n_,l_,m_},lcap_]:=-(l+lp+lcap+2*m+2*mp+2)/2*Hypergeometric1F1[(lcap-lp-l-2*mp-2*m-1)/2,lcap+3/2,y]+2*m*Hypergeometric1F1[(lcap-lp-l-2*mp-2*m+1)/2,lcap+3/2,y]


BesselFactor3A[y_,{np_,lp_},{n_,l_},lcap_]:=Sum[(Summand1[y,{np,lp,mp},{n,l,m},lcap]*Summand2A[y,{np,lp,mp},{n,l,m},lcap]*Summand3A[y,{np,lp,mp},{n,l,m},lcap]),{m,0,n-1},{mp,0,np-1}]

(**)

BesselElementMinus[y_,{np_,lp_},{n_,l_},lcap_]:=BesselFactor1A[y,{np,lp},{n,l},lcap]*BesselFactor2[y,{np,lp},{n,l},lcap]*BesselFactor3A[y,{np,lp},{n,l},lcap]


(**)

Summand4A[y_,{np_,lp_,mp_},{n_,l_,m_},lcap_]:=-(l+lp+lcap+2*m+2*mp+2)/2*Hypergeometric1F1[(lcap-lp-l-2*mp-2*m-1)/2,lcap+3/2,y]+(2*l+2*m+1)*Hypergeometric1F1[(lcap-lp-l-2*mp-2*m+1)/2,lcap+3/2,y]

BesselFactor4A[y_,{np_,lp_},{n_,l_},lcap_]:=Sum[(Summand1[y,{np,lp,mp},{n,l,m},lcap]*Summand2A[y,{np,lp,mp},{n,l,m},lcap]*Summand4A[y,{np,lp,mp},{n,l,m},lcap]),{m,0,n-1},{mp,0,np-1}]

BesselElementPlus[y_,{np_,lp_},{n_,l_},lcap_]:=BesselFactor1A[y,{np,lp},{n,l},lcap]*BesselFactor2[y,{np,lp},{n,l},lcap]*BesselFactor4A[y,{np,lp},{n,l},lcap]


Lnumber[NPrincipal_,j_]:=Which[EvenQ[NPrincipal-(j+1/2)],(j+1/2),True,j-1/2]

Nodal[NPrincipal_,j_]:=(NPrincipal-Lnumber[NPrincipal,j])/2+1


ParityState[NPrincipal_,j_]:=minus^(Lnumber[NPrincipal,j])


ParityNormal[jcap_]:=minus^(jcap)

ParityConsNormal[NPrincipalp_,jp_,NPrincipal_,j_,jcap_]:=ParityState[NPrincipal,j]*ParityState[NPrincipalp,jp]*ParityNormal[jcap]


PhysicalConditionsNormal[NPrincipalp_,jp_,NPrincipal_,j_,jcap_]:=ParityConsNormal[NPrincipal,j,NPrincipalp,jp,jcap]==1&&Nodal[NPrincipal,j]>0&&Nodal[NPrincipalp,jp]>0&&Abs[j-jp]<=jcap<=j+jp


ParityAbnormal[jcap_]:=minus^(jcap+1)

ParityConsAbnormal[NPrincipalp_,jp_,NPrincipal_,j_,jcap_]:=ParityState[NPrincipal,j]*ParityState[NPrincipalp,jp]*ParityAbnormal[jcap]


PhysicalConditionsAbnormal[NPrincipalp_,jp_,NPrincipal_,j_,jcap_]:=ParityConsAbnormal[NPrincipal,j,NPrincipalp,jp,jcap]==1&&Nodal[NPrincipal,j]>0&&Nodal[NPrincipalp,jp]>0&&Abs[j-jp]<=jcap<=j+jp


(* ::Text:: *)
(*Operators: normal parity*)
(*MJ: *)


mjelement[y_,{np_,lp_,jp_},{n_,l_,j_},jcap_]:=minus^(1/2+j+jcap) Sqrt[(JNorm[j]*JNorm[jp]*JNorm[l]*JNorm[lp]*JNorm[jcap])/(4 Pi)]*ThreeJ[{lp,0},{jcap,0},{l,0}]*SixJ[{lp,jp,1/2},{j,l,jcap}]*BesselElement[y,{np,lp},{n,l},jcap]


MJE[y_,{ncapp_,jp_},{ncap_,j_},jcap_]:=mjelement[y,{Nodal[ncapp,jp],Lnumber[ncapp,jp],jp},{Nodal[ncap,j],Lnumber[ncap,j],j},jcap]


MJ[y_,{ncapp_,jp_},{ncap_,j_},jcap_]:=If[PhysicalConditionsNormal[ncapp,jp,ncap,j,jcap],MJE[y,{ncapp,jp},{ncap,j},jcap],0]


(* ::Text:: *)
(*SigmaJ:*)


MJLSigma[y_,{np_,lp_,jp_},{n_,l_,j_},jcap_,lcap_]:=minus^lp Sqrt[(JNorm[j]*JNorm[jp]*JNorm[l]*JNorm[lp]*JNorm[jcap]*JNorm[lcap])/(4 Pi)]*Sqrt[6]*ThreeJ[{lp,0},{lcap,0},{l,0}]*NineJSymbol[{lp,l,lcap},{1/2,1/2,1},{jp,j,jcap}]*BesselElement[y,{np,lp},{n,l},lcap]


SigmaJE[y_,{ncapp_,jp_},{ncap_,j_},jcap_]:=MJLSigma[y,{Nodal[ncapp,jp],Lnumber[ncapp,jp],jp},{Nodal[ncap,j],Lnumber[ncap,j],j},jcap,jcap]


SigmaJ[y_,{ncapp_,jp_},{ncap_,j_},jcap_]:=If[PhysicalConditionsNormal[ncapp,jp,ncap,j,jcap],SigmaJE[y,{ncapp,jp},{ncap,j},jcap],0]


(* ::Text:: *)
(*DeltaPJ:*)


MJLDivQoverall[y_,{np_,lp_,jp_},{n_,l_,j_},jcap_,lcap_]:=minus^(lcap+j+1/2)*QNorm[lp]*QNorm[jp]*QNorm[j]*QNorm[jcap]*QNorm[lcap]*SixJ[{lp,jp,1/2},{j,l,jcap}]/Sqrt[4 Pi]

MJLDivQsummand1[y_,{np_,lp_,jp_},{n_,l_,j_},jcap_,lcap_]:=-Sqrt[l+1] QNorm[l+1] SixJ[{lcap,1,jcap},{l,lp,l+1}] ThreeJ[{lp,0},{lcap,0},{l+1,0}] BesselElementMinus[y,{np,lp},{n,l},lcap]

MJLDivQsummand2[y_,{np_,lp_,jp_},{n_,l_,j_},jcap_,lcap_]:=Sqrt[l] QNorm[l-1] SixJ[{lcap,1,jcap},{l,lp,l-1}] ThreeJ[{lp,0},{lcap,0},{l-1,0}] BesselElementPlus[y,{np,lp},{n,l},lcap]

MJLDivQ[y_,{np_,lp_,jp_},{n_,l_,j_},jcap_,lcap_]:=MJLDivQoverall[y,{np,lp,jp},{n,l,j},jcap,lcap]*(MJLDivQsummand1[y,{np,lp,jp},{n,l,j},jcap,lcap]+MJLDivQsummand2[y,{np,lp,jp},{n,l,j},jcap,lcap])

DeltaPrime[y_,{np_,lp_,jp_},{n_,l_,j_},jcap_]:=-Sqrt[jcap]/QNorm[jcap]*MJLDivQ[y,{np,lp,jp},{n,l,j},jcap,jcap+1]+Sqrt[jcap+1]/QNorm[jcap]*MJLDivQ[y,{np,lp,jp},{n,l,j},jcap,jcap-1]


DeltaPrimeJE[y_,{ncapp_,jp_},{ncap_,j_},jcap_]:=DeltaPrime[y,{Nodal[ncapp,jp],Lnumber[ncapp,jp],jp},{Nodal[ncap,j],Lnumber[ncap,j],j},jcap]


DeltaPJ[y_,{ncapp_,jp_},{ncap_,j_},jcap_]:=If[PhysicalConditionsNormal[ncapp,jp,ncap,j,jcap],DeltaPrimeJE[y,{ncapp,jp},{ncap,j},jcap],0]


(* ::Text:: *)
(*DeltaPPJ and DeltaTPPJ:*)


DeltaPrimePrime[y_,{np_,lp_,jp_},{n_,l_,j_},jcap_]:=Sqrt[jcap+1]/QNorm[jcap]*MJLDivQ[y,{np,lp,jp},{n,l,j},jcap,jcap+1]+Sqrt[jcap]/QNorm[jcap]*MJLDivQ[y,{np,lp,jp},{n,l,j},jcap,jcap-1]


DeltaPrimePrimeJE[y_,{ncapp_,jp_},{ncap_,j_},jcap_]:=DeltaPrimePrime[y,{Nodal[ncapp,jp],Lnumber[ncapp,jp],jp},{Nodal[ncap,j],Lnumber[ncap,j],j},jcap]


DeltaPPJ[y_,{ncapp_,jp_},{ncap_,j_},jcap_]:=If[PhysicalConditionsNormal[ncapp,jp,ncap,j,jcap],DeltaPrimePrimeJE[y,{ncapp,jp},{ncap,j},jcap],0]


DeltaTPPJ[y_,{ncapp_,jp_},{ncap_,j_},jcap_]:=DeltaPPJ[y,{ncapp,jp},{ncap,j},jcap]-MJ[y,{ncapp,jp},{ncap,j},jcap]/2


(* ::Text:: *)
(*PhiPPJ:*)


PhiPPoverall[y_,{np_,lp_,jp_},{n_,l_,j_},jcap_]:=minus^(lp+1)*6 *QNorm[lp]*QNorm[jp]*QNorm[j]/Sqrt[4 Pi]

PhiPPsummand1[y_,{np_,lp_,jp_},{n_,l_,j_},jcap_]:=QNorm[l+1]*Sqrt[l+1]*ThreeJ[{lp,0},{jcap+1,0},{l+1,0}]*BesselElementMinus[y,{np,lp},{n,l},jcap+1]*Sum[minus^(jcap-lcap+1)*JNorm[lcap]*SixJ[{jcap+1,1,lcap},{1,jcap,1}]*SixJ[{jcap+1,1,lcap},{l,lp,l+1}]*NineJSymbol[{lp,l,lcap},{1/2,1/2,1},{jp,j,jcap}],{lcap,jcap,jcap+1}]

PhiPPsummand2[y_,{np_,lp_,jp_},{n_,l_,j_},jcap_]:=If[l==0,0,QNorm[l-1]* Sqrt[l]*ThreeJ[{lp,0},{jcap+1,0},{l-1,0}]*BesselElementPlus[y,{np,lp},{n,l},jcap+1]*Sum[minus^(jcap-lcap)*JNorm[lcap]*SixJ[{jcap+1,1,lcap},{1,jcap,1}]*SixJ[{jcap+1,1,lcap},{l,lp,l-1}]*NineJSymbol[{lp,l,lcap},{1/2,1/2,1},{jp,j,jcap}],{lcap,jcap,jcap+1}]]

PhiPPsummand3[y_,{np_,lp_,jp_},{n_,l_,j_},jcap_]:=If[jcap==0,0,QNorm[l+1]*Sqrt[l+1]*ThreeJ[{lp,0},{jcap-1,0},{l+1,0}]*BesselElementMinus[y,{np,lp},{n,l},jcap-1]*Sum[minus^(jcap-lcap+1)*JNorm[lcap]*SixJ[{jcap-1,1,lcap},{1,jcap,1}]*SixJ[{jcap-1,1,lcap},{l,lp,l+1}]*NineJSymbol[{lp,l,lcap},{1/2,1/2,1},{jp,j,jcap}],{lcap,jcap-1,jcap}]]

PhiPPsummand4[y_,{np_,lp_,jp_},{n_,l_,j_},jcap_]:=If[jcap==0,0,If[l==0,0,QNorm[l-1]* Sqrt[l]*ThreeJ[{lp,0},{jcap-1,0},{l-1,0}]*BesselElementPlus[y,{np,lp},{n,l},jcap-1]*Sum[minus^(jcap-lcap)*JNorm[lcap]*SixJ[{jcap-1,1,lcap},{1,jcap,1}]*SixJ[{jcap-1,1,lcap},{l,lp,l-1}]*NineJSymbol[{lp,l,lcap},{1/2,1/2,1},{jp,j,jcap}],{lcap,jcap-1,jcap}]]]

PhiPPJX[y_,{np_,lp_,jp_},{n_,l_,j_},jcap_]:=PhiPPoverall[y,{np,lp,jp},{n,l,j},jcap]*(QNorm[jcap+1]*Sqrt[jcap+1]*(PhiPPsummand1[y,{np,lp,jp},{n,l,j},jcap]+PhiPPsummand2[y,{np,lp,jp},{n,l,j},jcap])+If[jcap==0,0,QNorm[jcap-1]*Sqrt[jcap]*(PhiPPsummand3[y,{np,lp,jp},{n,l,j},jcap]+PhiPPsummand4[y,{np,lp,jp},{n,l,j},jcap])])


PhiPPJE[y_,{ncapp_,jp_},{ncap_,j_},jcap_]:=PhiPPJX[y,{Nodal[ncapp,jp],Lnumber[ncapp,jp],jp},{Nodal[ncap,j],Lnumber[ncap,j],j},jcap]


PhiPPJ[y_,{ncapp_,jp_},{ncap_,j_},jcap_]:=If[PhysicalConditionsNormal[ncapp,jp,ncap,j,jcap],PhiPPJE[y,{ncapp,jp},{ncap,j},jcap],0]


(* ::Text:: *)
(*PhiPJ and PhiTPJ:*)


PhiPJX[y_,{np_,lp_,jp_},{n_,l_,j_},jcap_]:=PhiPPoverall[y,{np,lp,jp},{n,l,j},jcap]*(-QNorm[jcap+1]*Sqrt[jcap]*(PhiPPsummand1[y,{np,lp,jp},{n,l,j},jcap]+PhiPPsummand2[y,{np,lp,jp},{n,l,j},jcap])+If[jcap==0,0,QNorm[jcap-1]*Sqrt[jcap+1]*(PhiPPsummand3[y,{np,lp,jp},{n,l,j},jcap]+PhiPPsummand4[y,{np,lp,jp},{n,l,j},jcap])])


PhiPJE[y_,{ncapp_,jp_},{ncap_,j_},jcap_]:=PhiPJX[y,{Nodal[ncapp,jp],Lnumber[ncapp,jp],jp},{Nodal[ncap,j],Lnumber[ncap,j],j},jcap]


PhiPJ[y_,{ncapp_,jp_},{ncap_,j_},jcap_]:=If[PhysicalConditionsNormal[ncapp,jp,ncap,j,jcap],PhiPJE[y,{ncapp,jp},{ncap,j},jcap],0]


PhiTPJ[y_,{ncapp_,jp_},{ncap_,j_},jcap_]:=PhiPJ[y,{ncapp,jp},{ncap,j},jcap]+SigmaJ[y,{ncapp,jp},{ncap,j},jcap]/2


(* ::Text:: *)
(*Operators: abnormal parity*)
(*DeltaJ:*)


Deltaop[y_,{np_,lp_,jp_},{n_,l_,j_},jcap_]:=MJLDivQ[y,{np,lp,jp},{n,l,j},jcap,jcap]


DeltaJE[y_,{ncapp_,jp_},{ncap_,j_},jcap_]:=Deltaop[y,{Nodal[ncapp,jp],Lnumber[ncapp,jp],jp},{Nodal[ncap,j],Lnumber[ncap,j],j},jcap]


DeltaJ[y_,{ncapp_,jp_},{ncap_,j_},jcap_]:=If[PhysicalConditionsAbnormal[ncapp,jp,ncap,j,jcap],DeltaJE[y,{ncapp,jp},{ncap,j},jcap],0]


(* ::Text:: *)
(*SigmaPJ*)


SigmaPrime[y_,{np_,lp_,jp_},{n_,l_,j_},jcap_]:=-Sqrt[jcap]/QNorm[jcap]*MJLSigma[y,{np,lp,jp},{n,l,j},jcap,jcap+1]+Sqrt[jcap+1]/QNorm[jcap]*MJLSigma[y,{np,lp,jp},{n,l,j},jcap,jcap-1]


SigmaPrimeJE[y_,{ncapp_,jp_},{ncap_,j_},jcap_]:=SigmaPrime[y,{Nodal[ncapp,jp],Lnumber[ncapp,jp],jp},{Nodal[ncap,j],Lnumber[ncap,j],j},jcap]


SigmaPJ[y_,{ncapp_,jp_},{ncap_,j_},jcap_]:=If[PhysicalConditionsAbnormal[ncapp,jp,ncap,j,jcap],SigmaPrimeJE[y,{ncapp,jp},{ncap,j},jcap],0]


(* ::Text:: *)
(*SigmaPPJ*)


SigmaSecond[y_,{np_,lp_,jp_},{n_,l_,j_},jcap_]:=(Sqrt[jcap+1] MJLSigma[y,{np,lp,jp},{n,l,j},jcap,jcap+1]+Sqrt[jcap] MJLSigma[y,{np,lp,jp},{n,l,j},jcap,jcap-1])/QNorm[jcap]


SigmaSecondJE[y_,{ncapp_,jp_},{ncap_,j_},jcap_]:=SigmaSecond[y,{Nodal[ncapp,jp],Lnumber[ncapp,jp],jp},{Nodal[ncap,j],Lnumber[ncap,j],j},jcap]


SigmaPPJ[y_,{ncapp_,jp_},{ncap_,j_},jcap_]:=If[PhysicalConditionsAbnormal[ncapp,jp,ncap,j,jcap],SigmaSecondJE[y,{ncapp,jp},{ncap,j},jcap],0]


(* ::Text:: *)
(*OmegaJ and OmegaTJ:*)


MJLDivSigoverall[y_,{np_,lp_,jp_},{n_,l_,j_},jcap_]:=minus^(lp)*QNorm[lp]*QNorm[jp]*QNorm[j]*QNorm[2 j-l]*QNorm[jcap]*SixJ[{lp,jp,1/2},{j,2j-l,jcap}]*ThreeJ[{lp,0},{jcap,0},{2j-l,0}]/Sqrt[4 Pi]

MJLDivSigsummand1[y_,{np_,lp_,jp_},{n_,l_,j_},jcap_]:=-KroneckerDelta[j,l+1/2]*BesselElementMinus[y,{np,lp},{n,l},jcap]

MJLDivSigsummand2[y_,{np_,lp_,jp_},{n_,l_,j_},jcap_]:=KroneckerDelta[j,l-1/2]*BesselElementPlus[y,{np,lp},{n,l},jcap]

MJLDivSig[y_,{np_,lp_,jp_},{n_,l_,j_},jcap_]:=MJLDivSigoverall[y,{np,lp,jp},{n,l,j},jcap]*(MJLDivSigsummand1[y,{np,lp,jp},{n,l,j},jcap]+MJLDivSigsummand2[y,{np,lp,jp},{n,l,j},jcap])


OmegaJE[y_,{ncapp_,jp_},{ncap_,j_},jcap_]:=MJLDivSig[y,{Nodal[ncapp,jp],Lnumber[ncapp,jp],jp},{Nodal[ncap,j],Lnumber[ncap,j],j},jcap]


OmegaJ[y_,{ncapp_,jp_},{ncap_,j_},jcap_]:=If[PhysicalConditionsAbnormal[ncapp,jp,ncap,j,jcap],OmegaJE[y,{ncapp,jp},{ncap,j},jcap],0]


OmegaTJ[y_,{ncapp_,jp_},{ncap_,j_},jcap_]:=OmegaJ[y,{ncapp,jp},{ncap,j},jcap]+SigmaPPJ[y,{ncapp,jp},{ncap,j},jcap]/2


(* ::Text:: *)
(*PhiJ and PhiTJ:*)


Phioverall[y_,{np_,lp_,jp_},{n_,l_,j_},jcap_]:=minus^(lp+1)*6 *QNorm[lp]*QNorm[jp]*QNorm[j]*JNorm[jcap]/Sqrt[4 Pi]

Phisummand1[y_,{np_,lp_,jp_},{n_,l_,j_},jcap_]:=QNorm[l+1]*Sqrt[l+1]*ThreeJ[{lp,0},{jcap,0},{l+1,0}]*BesselElementMinus[y,{np,lp},{n,l},jcap]*Sum[minus^(jcap-lcap+1)*JNorm[lcap]*SixJ[{jcap,1,lcap},{1,jcap,1}]*SixJ[{jcap,1,lcap},{l,lp,l+1}]*NineJSymbol[{lp,l,lcap},{1/2,1/2,1},{jp,j,jcap}],{lcap,Max[0,jcap-1],jcap+1}]

Phisummand2[y_,{np_,lp_,jp_},{n_,l_,j_},jcap_]:=If[l==0,0,QNorm[l-1]* Sqrt[l]*ThreeJ[{lp,0},{jcap,0},{l-1,0}]*BesselElementPlus[y,{np,lp},{n,l},jcap]*Sum[minus^(jcap-lcap)*JNorm[lcap]*SixJ[{jcap,1,lcap},{1,jcap,1}]*SixJ[{jcap,1,lcap},{l,lp,l-1}]*NineJSymbol[{lp,l,lcap},{1/2,1/2,1},{jp,j,jcap}],{lcap,Max[0,jcap-1],jcap+1}]]

PhiJX[y_,{np_,lp_,jp_},{n_,l_,j_},jcap_]:=Phioverall[y,{np,lp,jp},{n,l,j},jcap]*(Phisummand1[y,{np,lp,jp},{n,l,j},jcap]+Phisummand2[y,{np,lp,jp},{n,l,j},jcap])


PhiJE[y_,{ncapp_,jp_},{ncap_,j_},jcap_]:=PhiJX[y,{Nodal[ncapp,jp],Lnumber[ncapp,jp],jp},{Nodal[ncap,j],Lnumber[ncap,j],j},jcap]


PhiJ[y_,{ncapp_,jp_},{ncap_,j_},jcap_]:=If[PhysicalConditionsAbnormal[ncapp,jp,ncap,j,jcap],PhiJE[y,{ncapp,jp},{ncap,j},jcap],0]


PhiTJ[y_,{ncapp_,jp_},{ncap_,j_},jcap_]:=PhiJ[y,{ncapp,jp},{ncap,j},jcap]-SigmaPJ[y,{ncapp,jp},{ncap,j},jcap]/2


DMLoad[ISM_]:=Module[{},
Clear[DM];Do[DMin=ReadList[DMstring[ISM,iFS],Number,RecordLists->True];
LcMax=Length[DMin];
jj=0;
Do[
If[Length[DMin[[Lc]]]==6,(jj=jj+1;NumDMs[iFS]=jj;kk=0;JF[iFS]=DMin[[Lc,1]]/2;TF[iFS]=DMin[[Lc,2]]/2;JI[iFS]=DMin[[Lc,3]]/2;TI[iFS]=DMin[[Lc,4]]/2;J0[iFS,jj]=DMin[[Lc,5]]/2;T0[iFS,jj]=DMin[[Lc,6]]/2;
If[JF[iFS]==JSF[[iFS]],True,Abort[]];
If[TF[iFS]==TMTS[[1]],True,Abort[]];
If[JI[iFS]==JSI[[iFS]],True,Abort[]];
If[TI[iFS]==TMTS[[1]],True,Abort[]]),True]
If[Length[DMin[[Lc]]]==5,(kk=kk+1;LDM[iFS,jj]=kk;
NB[iFS,jj,kk]=DMin[[Lc,1]];JB[iFS,jj,kk]=DMin[[Lc,2]]/2;NK[iFS,jj,kk]=DMin[[Lc,3]];JK[iFS,jj,kk]=DMin[[Lc,4]]/2;DMval[iFS,jj,kk]=DMin[[Lc,5]]),True],
{Lc,1,LcMax}];
Do[DM[iFS,jj]=Table[{{NB[iFS,jj,kk],JB[iFS,jj,kk]},{NK[iFS,jj,kk],JK[iFS,jj,kk]},DMval[iFS,jj,kk]},{kk,1,LDM[iFS,jj]}],{jj,1,NumDMs[iFS]}];,
{iFS,1,NumFS}]]


ISOT[TT_,MT_,T0_]:=KroneckerDelta[T0,0]Sqrt[2]/QNorm[TT]+KroneckerDelta[T0,1]  If[MT==0,0,MT Sqrt[6 /((2TT+1)(TT+1)TT)]]


(* ::Text:: *)
(*Elastic responses*)


FM[iFS_,jj_,y_]:=ISOT[TMTS[[1]],TMTS[[2]],T0[iFS,jj]]Sum[DM[iFS,jj][[kk,3]] MJ[y,DM[iFS,jj][[kk,1]],DM[iFS,jj][[kk,2]],J0[iFS,jj]],{kk,1,LDM[iFS,jj]}]


FPhiPP[ii_,jj_,y_]:=ISOT[TMTS[[1]],TMTS[[2]],T0[ii,jj]]Sum[DM[ii,jj][[kk,3]] PhiPPJ[y,DM[ii,jj][[kk,1]],DM[ii,jj][[kk,2]],J0[ii,jj]],{kk,1,LDM[ii,jj]}]


FPhiTP[ii_,jj_,y_]:=ISOT[TMTS[[1]],TMTS[[2]],T0[ii,jj]]Sum[DM[ii,jj][[kk,3]] PhiTPJ[y,DM[ii,jj][[kk,1]],DM[ii,jj][[kk,2]],J0[ii,jj]],{kk,1,LDM[ii,jj]}]


FSigmaP[ii_,jj_,y_]:=ISOT[TMTS[[1]],TMTS[[2]],T0[ii,jj]]Sum[DM[ii,jj][[kk,3]] SigmaPJ[y,DM[ii,jj][[kk,1]],DM[ii,jj][[kk,2]],J0[ii,jj]],{kk,1,LDM[ii,jj]}]


FSigmaPP[ii_,jj_,y_]:=ISOT[TMTS[[1]],TMTS[[2]],T0[ii,jj]]Sum[DM[ii,jj][[kk,3]] SigmaPPJ[y,DM[ii,jj][[kk,1]],DM[ii,jj][[kk,2]],J0[ii,jj]],{kk,1,LDM[ii,jj]}]


FDelta[ii_,jj_,y_]:=ISOT[TMTS[[1]],TMTS[[2]],T0[ii,jj]]Sum[DM[ii,jj][[kk,3]] DeltaJ[y,DM[ii,jj][[kk,1]],DM[ii,jj][[kk,2]],J0[ii,jj]],{kk,1,LDM[ii,jj]}]


(* ::Text:: *)
(*Elastic CLFV response functions: Direct*)


RM[i_,j_,cs_]:= cs[1,i]Conjugate[cs[1,j]] + cs[11,i] Conjugate[cs[11,j]]


RPhiPP[i_,j_,cs_]:=  cs[3,i] Conjugate[cs[3,j]] + (cs[12,i]- cs[15,i]) Conjugate[cs[12,j]- cs[15,j]]


RPhiTP[i_,j_,cs_]:= cs[12,i] Conjugate[cs[12,j]] + cs[13,i] Conjugate[cs[13,j]]


RSigmaP[i_,j_,cs_]:=cs[4,i] Conjugate[cs[4,j]]+  cs[9,i] Conjugate[cs[9,j]] 


RSigmaPP[i_,j_,cs_]:=cs[10,i] Conjugate[cs[10,j]] +(cs[4,i]-cs[6,i])Conjugate[ cs[4,j]-cs[6,j]] 


RDelta[i_,j_,cs_]:= cs[5,i] Conjugate[cs[5,j]] +cs[8,i] Conjugate[cs[8,j]]


(* ::Text:: *)
(*Elastic CLFV response functions : Interference*)


RPhiPPM[i_,j_,cs_]:=Re[ cs[3,i]Conjugate[cs[1,j]] - (cs[12,i]-cs[15,i]) Conjugate[cs[11,j]]]


RDeltaSigmaP[i_,j_,cs_]:= Re[cs[5,i] Conjugate[cs[4,j]]+ cs[8,i]Conjugate[cs[9,j]]]


(* ::Text:: *)
(*Elastic nuclear responses: Direct*)


WM[iFS_,qm_,y_,cs_]:=(4 Pi)/(2JI[iFS]+1)Sum[If[J0[iFS,jj]==J0[iFS,jjp]&&(RM[T0[iFS,jj],T0[iFS,jjp],cs]!=0.0),FM[iFS,jj,y] FM[iFS,jjp,y]RM[T0[iFS,jj],T0[iFS,jjp],cs],0],{jjp,1,NumDMs[iFS]},{jj,1,NumDMs[iFS]}]


WPhiPP[iFS_,qm_,y_,cs_]:=qm^2 *(4 Pi)/(2JI[iFS]+1)Sum[If[J0[iFS,jj]==J0[iFS,jjp]&&(RPhiPP[T0[iFS,jj],T0[iFS,jjp],cs]!=0.0),FPhiPP[iFS,jj,y] FPhiPP[iFS,jjp,y]RPhiPP[T0[iFS,jj],T0[iFS,jjp],cs],0],{jjp,1,NumDMs[iFS]},{jj,1,NumDMs[iFS]}]


WPhiTP[iFS_,qm_,y_,cs_]:=qm^2 *(4 Pi)/(2JI[iFS]+1)Sum[If[J0[iFS,jj]==J0[iFS,jjp]&&(RPhiTP[T0[iFS,jj],T0[iFS,jjp],cs]!=0.0),FPhiTP[iFS,jj,y] FPhiTP[iFS,jjp,y]RPhiTP[T0[iFS,jj],T0[iFS,jjp],cs],0],{jjp,1,NumDMs[iFS]},{jj,1,NumDMs[iFS]}]


WSigmaP[iFS_,qm_,y_,cs_]:=(4 Pi)/(2JI[iFS]+1)Sum[If[J0[iFS,jj]==J0[iFS,jjp]&&(RSigmaP[T0[iFS,jj],T0[iFS,jjp],cs]!=0.0),FSigmaP[iFS,jj,y] FSigmaP[iFS,jjp,y]RSigmaP[T0[iFS,jj],T0[iFS,jjp],cs],0],{jjp,1,NumDMs[iFS]},{jj,1,NumDMs[iFS]}]


WSigmaPP[iFS_,qm_,y_,cs_]:=(4 Pi)/(2JI[iFS]+1)Sum[If[J0[iFS,jj]==J0[iFS,jjp]&&(RSigmaPP[T0[iFS,jj],T0[iFS,jjp],cs]!=0.0),FSigmaPP[iFS,jj,y] FSigmaPP[iFS,jjp,y]RSigmaPP[T0[iFS,jj],T0[iFS,jjp],cs],0],{jjp,1,NumDMs[iFS]},{jj,1,NumDMs[iFS]}]


WDelta[iFS_,qm_,y_,cs_]:=qm^2 *(4 Pi)/(2JI[iFS]+1)Sum[If[J0[iFS,jj]==J0[iFS,jjp]&&(RDelta[T0[iFS,jj],T0[iFS,jjp],cs]!=0.0),FDelta[iFS,jj,y] FDelta[iFS,jjp,y]RDelta[T0[iFS,jj],T0[iFS,jjp],cs],0],{jjp,1,NumDMs[iFS]},{jj,1,NumDMs[iFS]}]


(* ::Text:: *)
(*Elastic nuclear responses: Interference*)


WPhiPPM[iFS_,qm_,y_,cs_]:=-2 qm *4 Pi/(2JI[iFS]+1)Sum[If[J0[iFS,jj]==J0[iFS,jjp]&&(RPhiPPM[T0[iFS,jj],T0[iFS,jjp],cs]!=0.0),FPhiPP[iFS,jj,y] FM[iFS,jjp,y]RPhiPPM[T0[iFS,jj],T0[iFS,jjp],cs],0],{jjp,1,NumDMs[iFS]},{jj,1,NumDMs[iFS]}]


WDeltaSigmaP[iFS_,qm_,y_,cs_]:=-2 qm *(4 Pi)/(2JI[iFS]+1)Sum[If[J0[iFS,jj]==J0[iFS,jjp]&&(RDeltaSigmaP[T0[iFS,jj],T0[iFS,jjp],cs]!=0.0),FDelta[iFS,jj,y] FSigmaP[iFS,jjp,y]RDeltaSigmaP[T0[iFS,jj],T0[iFS,jjp],cs],0],{jjp,1,NumDMs[iFS]},{jj,1,NumDMs[iFS]}]


(* ::Text:: *)
(*Inelastic responses*)


FSigma[ii_,jj_,y_]:=ISOT[TMTS[[1]],TMTS[[2]],T0[ii,jj]]Sum[DM[ii,jj][[kk,3]] SigmaJ[y,DM[ii,jj][[kk,1]],DM[ii,jj][[kk,2]],J0[ii,jj]],{kk,1,LDM[ii,jj]}]


FDeltaP[ii_,jj_,y_]:=ISOT[TMTS[[1]],TMTS[[2]],T0[ii,jj]]Sum[DM[ii,jj][[kk,3]] DeltaPJ[y,DM[ii,jj][[kk,1]],DM[ii,jj][[kk,2]],J0[ii,jj]],{kk,1,LDM[ii,jj]}]


FDeltaTPP[ii_,jj_,y_]:=ISOT[TMTS[[1]],TMTS[[2]],T0[ii,jj]]Sum[DM[ii,jj][[kk,3]] DeltaTPPJ[y,DM[ii,jj][[kk,1]],DM[ii,jj][[kk,2]],J0[ii,jj]],{kk,1,LDM[ii,jj]}]


FOmegaT[ii_,jj_,y_]:=ISOT[TMTS[[1]],TMTS[[2]],T0[ii,jj]]Sum[DM[ii,jj][[kk,3]] OmegaTJ[y,DM[ii,jj][[kk,1]],DM[ii,jj][[kk,2]],J0[ii,jj]],{kk,1,LDM[ii,jj]}]


FPhiT[ii_,jj_,y_]:=ISOT[TMTS[[1]],TMTS[[2]],T0[ii,jj]]Sum[DM[ii,jj][[kk,3]] PhiTJ[y,DM[ii,jj][[kk,1]],DM[ii,jj][[kk,2]],J0[ii,jj]],{kk,1,LDM[ii,jj]}]


(* ::Text:: *)
(*Inelastic CLFV response functions: Direct*)


ROmegaT[i_,j_,cs_]:= cs[7,i] Conjugate[cs[7,j]]+ cs[14,i]Conjugate[cs[14,j]]


RDeltaTPP[i_,j_,cs_]:= cs[2,i] Conjugate[cs[2,j]]+ (cs[8,i]-cs[16,i])*Conjugate[cs[8,j]-cs[16,j]]


(* ::Text:: *)
(*(Note we do not write RDeltaP, RSigma, and RPhiT because they are redundant with RDelta, RSigmaP, and RPhi, respectively.)*)


(* ::Text:: *)
(*Inelastic CLFV response functions: Interference*)


RDeltaTPPM[i_,j_,cs_]:=Im[ cs[2,i] Conjugate[cs[1,j]]- (cs[8,i]-cs[16,i])*Conjugate[cs[11,j]]]


RDeltaTPPPhiPP[i_,j_,cs_]:= Im[cs[2,i] Conjugate[cs[3,j]]+ (cs[8,i]-cs[16,i])*Conjugate[cs[12,j]-cs[15,j]]]


ROmegaTSigmaPP[i_,j_,cs_]:= Im[cs[7,i]*Conjugate[cs[10,j]]-cs[14,i]* Conjugate[cs[4,j]-cs[6,j] ]]


RDeltaPhiT[i_,j_,cs_]:=Im[ cs[5,i] Conjugate[cs[13,j]]+ cs[8,i]*Conjugate[cs[12,j]]]


RPhiTPSigma[i_,j_,cs_]:=Im[cs[13,i] Conjugate[cs[4,j]]+ cs[12,i]*Conjugate[cs[9,j]]]


(* ::Text:: *)
(*(Note we do not write RDeltaPSigma, RPhiTSigmaP, and RDeltaPPhiTP because they are redundant with RDeltaSigmaP, RPhiTPSigma, and RDeltaPhiT, respectively .)*)


(* ::Text:: *)
(*Inelastic nuclear responses: Direct*)


WOmegaT[iFS_,qm_,y_,cs_]:=qm^2 *(4 Pi)/(2JI[iFS]+1)Sum[If[J0[iFS,jj]==J0[iFS,jjp]&&(ROmegaT[T0[iFS,jj],T0[iFS,jjp],cs]!=0.0),FOmegaT[iFS,jj,y] FOmegaT[iFS,jjp,y]ROmegaT[T0[iFS,jj],T0[iFS,jjp],cs],0],{jjp,1,NumDMs[iFS]},{jj,1,NumDMs[iFS]}]


WDeltaTPP[iFS_,qm_,y_,cs_]:=qm^2 *(4 Pi)/(2JI[iFS]+1)Sum[If[J0[iFS,jj]==J0[iFS,jjp]&&(RDeltaTPP[T0[iFS,jj],T0[iFS,jjp],cs]!=0.0),FDeltaTPP[iFS,jj,y] FDeltaTPP[iFS,jjp,y]RDeltaTPP[T0[iFS,jj],T0[iFS,jjp],cs],0],{jjp,1,NumDMs[iFS]},{jj,1,NumDMs[iFS]}]


WDeltaP[iFS_,qm_,y_,cs_]:=qm^2 *(4 Pi)/(2JI[iFS]+1)Sum[If[J0[iFS,jj]==J0[iFS,jjp]&&(RDelta[T0[iFS,jj],T0[iFS,jjp],cs]!=0.0),FDeltaP[iFS,jj,y] FDeltaP[iFS,jjp,y]RDelta[T0[iFS,jj],T0[iFS,jjp],cs],0],{jjp,1,NumDMs[iFS]},{jj,1,NumDMs[iFS]}]


WSigma[iFS_,qm_,y_,cs_]:=(4 Pi)/(2JI[iFS]+1)Sum[If[J0[iFS,jj]==J0[iFS,jjp]&&(RSigmaP[T0[iFS,jj],T0[iFS,jjp],cs]!=0.0),FSigma[iFS,jj,y] FSigma[iFS,jjp,y]RSigmaP[T0[iFS,jj],T0[iFS,jjp],cs],0],{jjp,1,NumDMs[iFS]},{jj,1,NumDMs[iFS]}]


WPhiT[iFS_,qm_,y_,cs_]:=qm^2 *(4 Pi)/(2JI[iFS]+1)Sum[If[J0[iFS,jj]==J0[iFS,jjp]&&(RPhiTP[T0[iFS,jj],T0[iFS,jjp],cs]!=0.0),FPhiT[iFS,jj,y] FPhiT[iFS,jjp,y]RPhiTP[T0[iFS,jj],T0[iFS,jjp],cs],0],{jjp,1,NumDMs[iFS]},{jj,1,NumDMs[iFS]}]


(* ::Text:: *)
(*Inelastic nuclear responses: Interference*)


WDeltaTPPM[iFS_,qm_,y_,cs_]:=2 qm *4 Pi/(2JI[iFS]+1)Sum[If[J0[iFS,jj]==J0[iFS,jjp]&&(RDeltaTPPM[T0[iFS,jj],T0[iFS,jjp],cs]!=0.0),FDeltaTPP[iFS,jj,y] FM[iFS,jjp,y]RDeltaTPPM[T0[iFS,jj],T0[iFS,jjp],cs],0],{jjp,1,NumDMs[iFS]},{jj,1,NumDMs[iFS]}]


WDeltaTPPPhiPP[iFS_,qm_,y_,cs_]:=2 qm^2 *4 Pi/(2JI[iFS]+1)Sum[If[J0[iFS,jj]==J0[iFS,jjp]&&(RDeltaTPPPhiPP[T0[iFS,jj],T0[iFS,jjp],cs]!=0.0), FDeltaTPP[iFS,jj,y]FPhiPP[iFS,jjp,y]RDeltaTPPPhiPP[T0[iFS,jj],T0[iFS,jjp],cs],0],{jjp,1,NumDMs[iFS]},{jj,1,NumDMs[iFS]}]


WDeltaPSigma[iFS_,qm_,y_,cs_]:=2 qm *(4 Pi)/(2JI[iFS]+1)Sum[If[J0[iFS,jj]==J0[iFS,jjp]&&(RDeltaSigmaP[T0[iFS,jj],T0[iFS,jjp],cs]!=0.0),FDeltaP[iFS,jj,y] FSigma[iFS,jjp,y]RDeltaSigmaP[T0[iFS,jj],T0[iFS,jjp],cs],0],{jjp,1,NumDMs[iFS]},{jj,1,NumDMs[iFS]}]


WOmegaTSigmaPP[iFS_,qm_,y_,cs_]:=-2 qm *(4 Pi)/(2JI[iFS]+1)Sum[If[J0[iFS,jj]==J0[iFS,jjp]&&(ROmegaTSigmaPP[T0[iFS,jj],T0[iFS,jjp],cs]!=0.0), FOmegaT[iFS,jj,y]FSigmaPP[iFS,jjp,y]ROmegaTSigmaPP[T0[iFS,jj],T0[iFS,jjp],cs],0],{jjp,1,NumDMs[iFS]},{jj,1,NumDMs[iFS]}]


WDeltaPhiT[iFS_,qm_,y_,cs_]:=2 qm^2 *4 Pi/(2JI[iFS]+1)Sum[If[J0[iFS,jj]==J0[iFS,jjp]&&(RDeltaPhiT[T0[iFS,jj],T0[iFS,jjp],cs]!=0.0),FDelta[iFS,jj,y] FPhiT[iFS,jjp,y]RDeltaPhiT[T0[iFS,jj],T0[iFS,jjp],cs],0],{jjp,1,NumDMs[iFS]},{jj,1,NumDMs[iFS]}]


WDeltaPPhiTP[iFS_,qm_,y_,cs_]:=2 qm^2 *4 Pi/(2JI[iFS]+1)Sum[If[J0[iFS,jj]==J0[iFS,jjp]&&(RDeltaPhiT[T0[iFS,jj],T0[iFS,jjp],cs]!=0.0),FDeltaP[iFS,jj,y] FPhiTP[iFS,jjp,y]RDeltaPhiT[T0[iFS,jj],T0[iFS,jjp],cs],0],{jjp,1,NumDMs[iFS]},{jj,1,NumDMs[iFS]}]


WPhiTSigmaP[iFS_,qm_,y_,cs_]:=2 qm *(4 Pi)/(2JI[iFS]+1)Sum[If[J0[iFS,jj]==J0[iFS,jjp]&&(RPhiTPSigma[T0[iFS,jj],T0[iFS,jjp],cs]!=0.0), FPhiT[iFS,jj,y]FSigmaP[iFS,jjp,y]RPhiTPSigma[T0[iFS,jj],T0[iFS,jjp],cs],0],{jjp,1,NumDMs[iFS]},{jj,1,NumDMs[iFS]}]


WPhiTPSigma[iFS_,qm_,y_,cs_]:=-2 qm *(4 Pi)/(2JI[iFS]+1)Sum[If[J0[iFS,jj]==J0[iFS,jjp]&&(RPhiTPSigma[T0[iFS,jj],T0[iFS,jjp],cs]!=0.0), FPhiTP[iFS,jj,y]FSigma[iFS,jjp,y]RPhiTPSigma[T0[iFS,jj],T0[iFS,jjp],cs],0],{jjp,1,NumDMs[iFS]},{jj,1,NumDMs[iFS]}]


(* ::Text:: *)
(*Alternate expressions for DeltaTPP in terms of M, after applying current conservation*)


WDeltaTPPCC[iFS_,qm_,y_,cs_]:=1/qm^2*(\[CapitalDelta]ENuc[[iFS]]/mN)^2 *(4 Pi)/(2JI[iFS]+1)Sum[If[J0[iFS,jj]==J0[iFS,jjp]&&(RDeltaTPP[T0[iFS,jj],T0[iFS,jjp],cs]!=0.0),FM[iFS,jj,y] FM[iFS,jjp,y]RDeltaTPP[T0[iFS,jj],T0[iFS,jjp],cs],0],{jjp,1,NumDMs[iFS]},{jj,1,NumDMs[iFS]}]


WDeltaTPPPhiPPCC[iFS_,qm_,y_,cs_]:=-2 *(\[CapitalDelta]ENuc[[iFS]]/mN) *4 Pi/(2JI[iFS]+1)Sum[If[J0[iFS,jj]==J0[iFS,jjp]&&(RDeltaTPPPhiPP[T0[iFS,jj],T0[iFS,jjp],cs]!=0.0), FM[iFS,jj,y]FPhiPP[iFS,jjp,y]RDeltaTPPPhiPP[T0[iFS,jj],T0[iFS,jjp],cs],0],{jjp,1,NumDMs[iFS]},{jj,1,NumDMs[iFS]}]


WDeltaTPPMCC[iFS_,qm_,y_,cs_]:=2/qm *(\[CapitalDelta]ENuc[[iFS]]/mN) *4 Pi/(2JI[iFS]+1)Sum[If[J0[iFS,jj]==J0[iFS,jjp]&&(RDeltaTPPM[T0[iFS,jj],T0[iFS,jjp],cs]!=0.0),FM[iFS,jj,y] FM[iFS,jjp,y]RDeltaTPPM[T0[iFS,jj],T0[iFS,jjp],cs],0],{jjp,1,NumDMs[iFS]},{jj,1,NumDMs[iFS]}]


Heff[iFS_,qm_,y_,cs_]:=Module[{Wtot},
Wtot=0.0;

Wtot+=WM[iFS,qm,y,cs]+WSigmaPP[iFS,qm,y,cs]+WSigma[iFS,qm,y,cs]+WSigmaP[iFS,qm,y,cs];

Wtot+=WDeltaTPPCC[iFS,qm,y,cs]+WOmegaT[iFS,qm,y,cs]+WPhiPP[iFS,qm,y,cs]+WDeltaTPPPhiPPCC[iFS,qm,y,cs]+WDelta[iFS,qm,y,cs]+WDeltaP[iFS,qm,y,cs]+WPhiTP[iFS,qm,y,cs]+WPhiT[iFS,qm,y,cs]+WDeltaPhiT[iFS,qm,y,cs]+WDeltaPPhiTP[iFS,qm,y,cs];

Wtot+=WDeltaTPPMCC[iFS,qm,y,cs]+WPhiPPM[iFS,qm,y,cs]+WOmegaTSigmaPP[iFS,qm,y,cs]+WDeltaSigmaP[iFS,qm,y,cs]+WDeltaPSigma[iFS,qm,y,cs]+WPhiTPSigma[iFS,qm,y,cs]+WPhiTSigmaP[iFS,qm,y,cs];

Return[Wtot];
];


ShiftSignal[Mu2eResponse_,\[CapitalDelta]E_,BR_]:=Module[{nShift,Mu2eShift},
nShift=Ceiling[\[CapitalDelta]E/0.05];
Mu2eShift={};
Do[
If[is<Length[Mu2eResponse]-nShift,
AppendTo[Mu2eShift,{Mu2eResponse[[is,1]],PlusMinus[BR[[1]]*Mu2eResponse[[is+nShift,2,1]],Sqrt[BR[[2]]^2*Mu2eResponse[[is+nShift,2,1]]^2+BR[[1]]^2*Mu2eResponse[[is+nShift,2,2]]^2]]}];
,AppendTo[Mu2eShift,{Mu2eResponse[[is,1]],PlusMinus[0,0]}]];
,{is,1,Length[Mu2eResponse]}];
Return[Mu2eShift];
]
