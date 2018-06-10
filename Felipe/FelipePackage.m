(* ::Package:: *)

BeginPackage["FelipePackage`"];

Needs["Combinatorica`"]

setPaths[patientnumber_,expgroup_,badepochsrow_]:=(
maintag=patientnumber<>expgroup;
datapath="/Users/drwho/Documents/Research/Felipe/";
savepath="/Users/drwho/Documents/Research/Felipe/PSDmetadata/";
savename="patient"<>maintag<>"_HFonly.mx";
outputfilename = maintag<>"_pcorHFvsEEGband";
);

setConstants[]:=(
EEGsamplerate = 500;
RRresamplerate = 4;
minepoch = 15;
);

importData[]:=(
EKGraw=Flatten[Import[datapath<> maintag<>"_EKG_T_RR.tsv","Data"]];
EEGraw=Import[datapath<>maintag<>"_EEGpreprocessed_500.tsv","Data"];
EEGname=Import[datapath<>"CMIT_eeg_channels.csv","CSV"];
badtimes=Import[datapath<>"CMIT_eeg_bad_epochs.csv","CSV"];
EEGbadsamples=badtimes[[badepochsrow,All]]*EEGsamplerate;
);

getRRtimeseries[]:=(
RR = Differences[EKGraw];
RRts=Table[{Sum[RR[[i]],{i,1,j}],RR[[j]]},{j,1,Length[RR]}];
);

RRPower[timeseries_,resamplerate_]:=Module[{interpolated,resampled,lineartrend,detrended,PSD,PSDf,m,x,b,\[Omega]},
interpolated=Interpolation[timeseries,Method->"Spline",InterpolationOrder->3];
resampled=Table[(interpolated[x]-Mean[timeseries][[2]])/Mean[timeseries][[2]],{x,timeseries[[1]][[1]],timeseries[[-1]][[1]],1/resamplerate}];
lineartrend=FindFit[resampled,m x+b,{m,b},x];
detrended=Table[resampled[[i]]-(m i+b)/.lineartrend,{i,1,Length[resampled]}];
PSD=PowerSpectralDensity[detrended,\[Omega],HannWindow];
PSDf=PSD/.\[Omega]->2\[Pi] f/resamplerate;
Return[PSDf];
];

RRBandPSD[patienttime_,punchIn_,punchOut_]:=Module[{plot,VLF,LF,HF,BandPowers,PSD},
PSD = RRPower[RRts[[punchIn;;punchOut]],RRresamplerate];
plot = Plot[PSD,{f,0,.4},PlotRange->All,Frame->True,FrameLabel->{"frequency (Hz)","power spectral density"}];
VLF=NIntegrate[PSD,{f,.0033,.04},Method->{"Trapezoidal","SymbolicProcessing"->False},AccuracyGoal->2];
LF=NIntegrate[PSD,{f,.04,.15},Method->{"Trapezoidal","SymbolicProcessing"->False},AccuracyGoal->2];
HF=NIntegrate[PSD,{f,.15,.4},Method->{"Trapezoidal","SymbolicProcessing"->False},AccuracyGoal->2];
BandPowers = {VLF,LF,HF};
Return[{plot,BandPowers}]
];

nextindex[in_]:=in+Floor[BinarySearch[Table[RRts[[in+i,1]]-RRts[[in,1]],{i,1,15*5}],minepoch]];

baddata[entry_]:=Not[MatchQ[entry,{_?NumberQ,_?NumberQ}]]

getRRPunches[]:=(
elementend=(Length[Select[badtimes[[badepochsrow,All]],IntegerQ]]/2)-1;
punchInOut = Table[{2i,1+2i}+1,{i,1,elementend}];
punchtimes=Table[-(badtimes[[badepochsrow,i]]-badtimes[[badepochsrow,i+1]]),{i,punchInOut[[1,1]],punchInOut[[Length[punchInOut],1]],2}];
longepochs = Flatten[Position[punchtimes,_?(#>minepoch&)]];
punchFiltered = punchInOut[[longepochs]];
punchIn = Table[Floor[BinarySearch[RRts[[All,1]],badtimes[[badepochsrow,punchFiltered[[i,1]]]]]+1],{i,1,Length[punchFiltered]}];
punchOut = Table[Floor[BinarySearch[RRts[[All,1]],badtimes[[badepochsrow,punchFiltered[[i,2]]]]]],{i,1,Length[punchFiltered]}];
Clear[punches];
punchesall=Array[punches,Length[punchIn]];
For[i=1,i<=Length[punchIn],i++,
punches[i]={punchIn[[i]]};
While[Last[punches[i]]<=punchOut[[i]],
next=If[Last[punches[i]]+15*5<=punchOut[[i]],nextindex[Last[punches[i]]]];
AppendTo[punches[i],next];
];
];
RRpunches=Flatten[Table[{punches[j][[i]],punches[j][[i+1]]},{j,1,Length[punchIn]},{i,1,Length[punches[j]]-1}],1];
RRpunches=DeleteCases[RRpunches,_?baddata];
);

getEEGPunches[]:=(
EEGtimepunches=Floor[Table[{RRts[[RRpunches[[i,1]]]][[1]],RRts[[RRpunches[[i,2]]]][[1]]},{i,1,Length[RRpunches]}]*EEGsamplerate];
);

getRRPSD[]:=(
PSDresultRR=Table[RRBandPSD[badepochsrow,RRpunches[[i,1]],RRpunches[[i,2]]],{i,1,Length[RRpunches]}];
);

EEGPower[timeseries_,samplerate_,cutoff_]:=Module[{interpolated,resampled,lineartrend,detrended,PSD,PSDf,m,x,b,\[Omega]},
lineartrend=FindFit[timeseries,m x+b,{m,b},x];
detrended=Table[timeseries[[i]]-(m i+b)/.lineartrend,{i,1,Length[timeseries]}];
PSD=PowerSpectralDensity[detrended,\[Omega],{cutoff,HannWindow}];
PSDf=PSD/.\[Omega]->2\[Pi] f/samplerate;
Return[PSDf];
];

EEGBandPSD[electrode_,insample_,outsample_]:=Module[{plot,delta,theta,alpha,beta,highbeta,gamma,BandPowers,PSD},
PSD = EEGPower[EEGraw[[electrode,insample;;outsample]],EEGsamplerate,3100];
plot = Plot[Log[PSD],{f,0,12},PlotRange->All,Frame->True,FrameLabel->{"frequency (Hz)","log power spectral density"}];
delta=NIntegrate[PSD,{f,1.5,4},Method->{"Trapezoidal","SymbolicProcessing"->False},AccuracyGoal->1];
theta=NIntegrate[PSD,{f,4,8},Method->{"Trapezoidal","SymbolicProcessing"->False},AccuracyGoal->1];
alpha=NIntegrate[PSD,{f,8,12},Method->{"Trapezoidal","SymbolicProcessing"->False},AccuracyGoal->1];
beta=NIntegrate[PSD,{f,12,20},Method->{"Trapezoidal","SymbolicProcessing"->False},AccuracyGoal->1];
highbeta=NIntegrate[PSD,{f,20,40},Method->{"Trapezoidal","SymbolicProcessing"->False},AccuracyGoal->1];
gamma=NIntegrate[PSD,{f,40,75},Method->{"Trapezoidal","SymbolicProcessing"->False},AccuracyGoal->1];
BandPowers = {delta,theta,alpha,beta,highbeta,gamma};
Return[{plot,BandPowers}]
];

getEEGPSD[]:=(
PSDresultEEG=Table[EEGBandPSD[electrode,EEGtimepunches[[i,1]],EEGtimepunches[[i,2]]],{electrode,3,40},{i,1,Length[EEGtimepunches]}]
);

getEEGpowerCor[]:=(
pearsonEEGPower = Table[Correlation[Flatten[PSDresultEEG[[i,All,2]]],Flatten[PSDresultEEG[[j,All,2]]]],{i,1,Dimensions[PSDresultEEG][[1]]},{j,1,Dimensions[PSDresultEEG][[1]]}];
);

(*HF ONLY*)
getRREEGCor[]:=(
pearsonRREEG=Table[Correlation[PSDresultRR[[All,2,RRband]],Flatten[PSDresultEEG[[electrode,All,2,EEGband]]]],{RRband,3,3},{EEGband,1,6},{electrode,1,Dimensions[PSDresultEEG][[1]]}];
);

initBandNames[]:=(
bandRR = {1,2,3};
bandRRname={"RR VLF","RR LF", "RR HF"};
bandEEG={1,2,3,4,5,6};
bandEEGname={"EEG \[Delta]","EEG \[Theta]","EEG \[Alpha]","EEG \[Beta]","EEG high\[Beta]","EEG \[Gamma]"};
electrode=Range[Dimensions[PSDresultEEG][[1]]];
);

RhoPlot[RRindex_,EEGindex_,electrode_]:=ListPlot[Table[{PSDresultRR[[All,2,bandRR[[RRindex]]]][[i]],PSDresultEEG[[electrode,All,2,bandEEG[[EEGindex]]]][[i]]},
{i,1,Length[PSDresultRR[[All,2,bandRR[[RRindex]]]]]}],
Frame->True,
PlotMarkers->{Automatic,Medium},
PlotLabel->
StringJoin["electrode=",ToString[EEGname[[electrode+2]][[1]]]," \[Rho]=",ToString[SpearmanRho[PSDresultRR[[All,2,bandRR[[RRindex]]]],PSDresultEEG[[electrode,All,2,bandEEG[[EEGindex]]]]]],"  p=",ToString[SpearmanRankTest[PSDresultRR[[All,2,bandRR[[RRindex]]]],PSDresultEEG[[electrode,All,2,bandEEG[[EEGindex]]]]]]],
FrameLabel->{bandRRname[[RRindex]],bandEEGname[[EEGindex]]},
LabelStyle->Medium]

saveSession[]:=(
DumpSave[savepath<>savename,{PSDresultRR,PSDresultEEG,pearsonEEGPower,pearsonRREEG}];
electrodeorder={14,24,18,20,23,25,22,26,19,8,10,7,11,13,15,12,16,3,4,6,9,34,36,35,29,31,39,40,37,38,30,17,21,28,32};
Table[Export[savepath<>outputfilename<>ToString[i]<>".csv",Table[{EEGname[[electrodeorder]][[j]][[1]],pearsonRREEG[[1,i,electrodeorder-2]][[j]],0,0,0},{j,1,35}]],{i,1,6}];
);

runFelipe[patnum_,group_,badepochs_]:=(
patientnumber=patnum;
expgroup=group;
badepochsrow=badepochs;
setPaths[patientnumber,expgroup,badepochsrow];
setConstants[];
importData[];
getRRtimeseries[];
getRRPunches[];
getEEGPunches[];
getRRPSD[];
getEEGPSD[];
getEEGpowerCor[];
getRREEGCor[];
initBandNames[];
saveSession[];
);

EndPackage[]



