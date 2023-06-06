(* ::Package:: *)

$Assumptions=p>6&&e>0;
SetDirectory[NotebookDirectory[]];
orbitdata=Import["orbit.txt","Table"];
x2=Interpolation[orbitdata[[All,{1,2}]],Method->"Spline",InterpolationOrder->6];
x3=Interpolation[orbitdata[[All,{1,3}]],Method->"Spline",InterpolationOrder->6];
x4=Interpolation[orbitdata[[All,{1,4}]],Method->"Spline",InterpolationOrder->6];
x5=Interpolation[orbitdata[[All,{1,5}]],Method->"Spline",InterpolationOrder->6];
x6=Interpolation[orbitdata[[All,{1,6}]],Method->"Spline",InterpolationOrder->6];
x7=Interpolation[orbitdata[[All,{1,7}]],Method->"Spline",InterpolationOrder->6];


pm=Import["pm.txt","Table"];
x8=Interpolation[pm[[All,{1,2}]],Method->"Spline",InterpolationOrder->6];
x9=Interpolation[pm[[All,{1,3}]],Method->"Spline",InterpolationOrder->6];
x10=Interpolation[pm[[All,{1,4}]],Method->"Spline",InterpolationOrder->6];
x11=Interpolation[pm[[All,{1,5}]],Method->"Spline",InterpolationOrder->6];
x12=Interpolation[pm[[All,{1,6}]],Method->"Spline",InterpolationOrder->6];
x13=Interpolation[pm[[All,{1,7}]],Method->"Spline",InterpolationOrder->6];


pmu=Import["pmu.txt","Table"];
x14=Interpolation[pmu[[All,{1,2}]],Method->"Spline",InterpolationOrder->6];
x15=Interpolation[pmu[[All,{1,3}]],Method->"Spline",InterpolationOrder->6];
x16=Interpolation[pmu[[All,{1,4}]],Method->"Spline",InterpolationOrder->6];
x17=Interpolation[pmu[[All,{1,5}]],Method->"Spline",InterpolationOrder->6];
x18=Interpolation[pmu[[All,{1,6}]],Method->"Spline",InterpolationOrder->6];
x19=Interpolation[pmu[[All,{1,7}]],Method->"Spline",InterpolationOrder->6];


pp=Import["pp.txt","Table"];
x20=Interpolation[pp[[All,{1,2}]],Method->"Spline",InterpolationOrder->6];
x21=Interpolation[pp[[All,{1,3}]],Method->"Spline",InterpolationOrder->6];
x22=Interpolation[pp[[All,{1,4}]],Method->"Spline",InterpolationOrder->6];
x23=Interpolation[pp[[All,{1,5}]],Method->"Spline",InterpolationOrder->6];
x24=Interpolation[pp[[All,{1,6}]],Method->"Spline",InterpolationOrder->6];
x25=Interpolation[pp[[All,{1,7}]],Method->"Spline",InterpolationOrder->6];


pe=Import["pe.txt","Table"];
x26=Interpolation[pe[[All,{1,2}]],Method->"Spline",InterpolationOrder->6];
x27=Interpolation[pe[[All,{1,3}]],Method->"Spline",InterpolationOrder->6];
x28=Interpolation[pe[[All,{1,4}]],Method->"Spline",InterpolationOrder->6];
x29=Interpolation[pe[[All,{1,5}]],Method->"Spline",InterpolationOrder->6];
x30=Interpolation[pe[[All,{1,6}]],Method->"Spline",InterpolationOrder->6];
x31=Interpolation[pe[[All,{1,7}]],Method->"Spline",InterpolationOrder->6];


pd=Import["pd.txt","Table"];
x32=Interpolation[pd[[All,{1,2}]],Method->"Spline",InterpolationOrder->6];
x33=Interpolation[pd[[All,{1,3}]],Method->"Spline",InterpolationOrder->6];
x34=Interpolation[pd[[All,{1,4}]],Method->"Spline",InterpolationOrder->6];
x35=Interpolation[pd[[All,{1,5}]],Method->"Spline",InterpolationOrder->6];
x36=Interpolation[pd[[All,{1,6}]],Method->"Spline",InterpolationOrder->6];
x37=Interpolation[pd[[All,{1,7}]],Method->"Spline",InterpolationOrder->6];


tmaxarray={orbitdata[[-1,1]],pd[[-1,1]],pe[[-1,1]],pm[[-1,1]],pmu[[-1,1]],pp[[-1,1]]};
duration=Min[tmaxarray];
deltat=10;
number=IntegerPart[duration/deltat]-1;
timetable=Table[i,{i,0,number-1}];
timelist=deltat*timetable;

\[Beta]0=0;
\[Alpha]0=0;
phit=\[Beta]0+(2Pi t)/YRSIDSI;
costheta=(1/2)Cos[\[Theta]s]-((\[Sqrt]3)/2)Sin[\[Theta]s]Cos[phit-\[Phi]s];
\[Phi]=\[Alpha]0+phit+ArcTan[(\[Sqrt]3Cos[\[Theta]s]+Sin[\[Theta]s]Cos[phit-\[Phi]s])/(2Sin[\[Theta]s]Sin[phit-\[Phi]s])];
tanpsi=((1/2)Cos[\[Theta]1]-((\[Sqrt]3)/2)Sin[\[Theta]1]Cos[phit-\[Phi]1]-(Cos[\[Theta]1]Cos[\[Theta]s]+Sin[\[Theta]1]Sin[\[Theta]s]Cos[\[Phi]1-\[Phi]s])costheta)/((1/2)Sin[\[Theta]1]Sin[\[Theta]s]Sin[\[Phi]1-\[Phi]s]-((\[Sqrt]3)/2)Cos[phit](Cos[\[Theta]1]Sin[\[Theta]s]Sin[\[Phi]s]-Cos[\[Theta]s]Sin[\[Theta]1]Sin[\[Phi]1])-((\[Sqrt]3)/2)Sin[phit](Cos[\[Theta]s]Sin[\[Theta]1]Cos[\[Phi]1]-Cos[\[Theta]1]Sin[\[Theta]s]Cos[\[Phi]s]));
Fplus1=((1+costheta^2)/2)Cos[2\[Phi]]((1-tanpsi^2)/(1+tanpsi^2))-costheta Sin[2\[Phi]]((2tanpsi)/(1+tanpsi^2));
Fcross1=((1+costheta^2)/2)Cos[2\[Phi]]((2tanpsi)/(1+tanpsi^2))+costheta Sin[2\[Phi]]((1-tanpsi^2)/(1+tanpsi^2));
Fplus2=((1+costheta^2)/2)Cos[2\[Phi]-Pi/2]((1-tanpsi^2)/(1+tanpsi^2))-costheta Sin[2\[Phi]-Pi/2]((2tanpsi)/(1+tanpsi^2));
Fcross2=((1+costheta^2)/2)Cos[2\[Phi]-Pi/2]((2tanpsi)/(1+tanpsi^2))+costheta Sin[2\[Phi]-Pi/2]((1-tanpsi^2)/(1+tanpsi^2));


anfinal=0;
bnfinal=0;
cnfinal=0;
gam=Phiphi-Phir;
cos2gam=Cos[2gam];
sin2gam=Sin[2gam];
For[n=1,n<=Max[Round[30*0.2],4],n++,
fn=n Phirdot/(2Pi)+(Phiphidot-Phirdot)/Pi;
Doppler=2*Pi*fn*AUsec/(M MTSUNSI)*Sin[\[Theta]s]*Cos[(2*Pi*t)/YRSIDSI-\[Phi]s];
nPhi=n*Phir+Doppler;
J0=BesselJ[n-2,n*e];
J1=BesselJ[n-1,n*e];
J2=BesselJ[n,n*e];
J3=BesselJ[n+1,n*e];
J4=BesselJ[n+2,n*e];
an=-n*(J0-2*e*J1+2/n*J2+2*e*J3-J4)*Cos[nPhi];
bn=-n*Sqrt[1-e*e]*(J0-2*J2+J4)*Sin[nPhi];
cn=2*J2*Cos[nPhi];
anfinal=anfinal+an;
bnfinal=bnfinal+bn;
cnfinal=cnfinal+cn;
];


Ldotn=Cos[\[Theta]1]Cos[\[Theta]s]+Sin[\[Theta]1]Sin[\[Theta]s]Cos[\[Phi]1-\[Phi]s];
Aplus=-(1+Ldotn^2)*(anfinal*cos2gam-bnfinal*sin2gam)+cnfinal*(1-Ldotn^2);
Acros=2*Ldotn*(bnfinal*cos2gam+anfinal*sin2gam);
Amp=Phiphidot^(2/3)*\[Mu]/dl MRSUNSI;
det1=((\[Sqrt]3)/2)(Fplus1*Amp*Aplus+Fcross1*Amp*Acros);
det2=((\[Sqrt]3)/2)(Fplus2*Amp*Aplus+Fcross2*Amp*Acros);
parametersvalues={M->10^6,\[Mu]->10,YRSIDSI->31558149.763545603,\[Theta]s->Pi/3,\[Phi]s->Pi/2,\[Theta]1->Pi/4,\[Phi]1->Pi/4,dl->3.0856775814913674*10^25,MTSUNSI->4.92549102587*10^-6,AUsec->499.004783836156412,MRSUNSI->1476.6250615036158,Gpc->3.0856775814913674*10^25};


parameters={M,\[Mu],p0,e0,d};
orbitparameters={e,Phiphi,Phir,Phiphidot,Phirdot};
dss=Table[D[det1,i],{i,orbitparameters}];
dss2=Table[D[det1,i],{i,parameters}];
resultdet1[i_]:=Module[{rule3v,pmatrix,va1,va2},
rule3v={p->x2[t],e->x3[t],Phiphi->x4[t],Phir->x5[t],Phiphidot->x6[t],Phirdot->x7[t]}/.t->i;
pmatrix=({
 {x9[t], x15[t], x21[t], x27[t], x33[t]},
 {x10[t], x16[t], x22[t], x28[t], x34[t]},
 {x11[t], x17[t], x23[t], x29[t], x35[t]},
 {x12[t], x18[t], x24[t], x30[t], x36[t]},
 {x13[t], x19[t], x25[t], x31[t], x37[t]}
})/.t->i;
va1=dss/.parametersvalues/.rule3v/.t->i;
va2=dss2/.parametersvalues/.rule3v/.t->i;
va1.pmatrix+va2
];

d2ss=Table[D[det2,i],{i,orbitparameters}];
d2ss2=Table[D[det2,i],{i,parameters}];
resultdet2[i_]:=Module[{rule3v,pmatrix,va1,va2},
rule3v={p->x2[t],e->x3[t],Phiphi->x4[t],Phir->x5[t],Phiphidot->x6[t],Phirdot->x7[t]}/.t->i;
pmatrix=({
 {x9[t], x15[t], x21[t], x27[t], x33[t]},
 {x10[t], x16[t], x22[t], x28[t], x34[t]},
 {x11[t], x17[t], x23[t], x29[t], x35[t]},
 {x12[t], x18[t], x24[t], x30[t], x36[t]},
 {x13[t], x19[t], x25[t], x31[t], x37[t]}
})/.t->i;
va1=d2ss/.parametersvalues/.rule3v/.t->i;
va2=d2ss2/.parametersvalues/.rule3v/.t->i;
va1.pmatrix+va2
];


signaldetector1=ParallelTable[resultdet1[i],{i,timelist}];
Export["scalar_signal_detector1_five.txt",signaldetector1,"Table"];
signaldetector2=ParallelTable[resultdet2[i],{i,timelist}];
Export["scalar_signal_detector2_five.txt",signaldetector2,"Table"];
Export["scalar_timelist_five.txt",timelist,"Table"];
