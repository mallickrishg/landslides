(* ::Package:: *)

(* ::Package:: *)
(**)


BeginPackage["LandslideProcessor`"];

processSheet::usage = "processSheet[sheet] processes a single Excel sheet of landslide data and returns rescaled quantities.";
processWorkbook::usage = "processWorkbook[file] processes all sheets in an Excel workbook and returns results.";

xPosRoot::usage =
  "xPosRoot[\[Alpha], f, b0, k] finds the positive root position x where H=0.";

xNegRoot::usage =
  "xNegRoot[\[Alpha], f, b0, k] finds the negative root position x where H=0.";

locusPosFun::usage =
  "locusPosFun[f, b0, k] computes \[Alpha] on the positive side stability locus.";

locusNegFun::usage =
  "locusNegFun[f, b0, k] computes \[Alpha] on the negative side stability locus.";

HSolution::usage =
  "HSolution[\[Alpha], f, b0, k,{xmin,xmax}] gives the landslide thickness as a function of x";

solveH::usage =
  "solveH[\[Alpha], f, b0, k,{xmin,xmax}] solves the governing ODE for landslide thickness";
   
Begin["`Private`"];

(* Private helper: solve the H[x] ODE for given parameters *)
solveH[\[Alpha]v_?NumericQ, fv_?NumericQ, b0_?NumericQ, k_?NumericQ, {xmin_, xmax_}] :=
  Quiet @ Check[
    NDSolve[{
      -D[H[x], x] -
        0.5*(\[Alpha]v - 2 b0 x)/(1 + (2 b0 x)^2) * (2 b0)/(1 + 2 b0 x \[Alpha]v) *
          H[x] == fv + 2 b0 x - (\[Alpha]v - 2 b0 x)/(1 + 2 b0 x \[Alpha]v),
      H[0] == k Abs[b0]
    },
    H, {x, xmin, xmax}],
    $Failed
  ];

processSheet[sheet_] := 
  Module[{ptsB, ptsH, xData, zData, linFit, \[Alpha], \[Alpha]Error,residuals, quadFit, 
    a, b, c, x0, Lx, \[Delta]b, b0, BRescaled, HRescaled, qF},
   
   ptsB = Select[sheet[[2 ;;, 1 ;; 2]], ! AllTrue[#, # === "" || # === Null &] &];
   ptsH = Select[sheet[[2 ;;, 3 ;; 4]], ! AllTrue[#, # === "" || # === Null &] &];
   
   xData = Join[ptsB[[All, 1]], ptsH[[All, 1]]];
   zData = Join[ptsB[[All, 2]], ptsH[[All, 2]]];
   
   (* Linear regression using LinearModelFit *)
   linFit = LinearModelFit[Transpose[{xData, zData}], x, x];
   \[Alpha]Error = linFit["ParameterErrors"][[2]];
   linFit = Fit[Transpose[{xData, zData}], {1, x}, x];
   \[Alpha] = Coefficient[linFit, x, 1];
   
   residuals = ptsB[[All, 2]] - (linFit /. x -> ptsB[[All, 1]]);
   
   quadFit = Fit[Transpose[{ptsB[[All, 1]], residuals}], {1, x, x^2}, x];
   a = Coefficient[quadFit, x, 0];
   b = Coefficient[quadFit, x, 1];
   c = Coefficient[quadFit, x, 2];
   
   x0 = -b/(2 c);
   Lx = Max[Abs[xData - x0]];
   \[Delta]b = quadFit /. x -> (x0 + Lx);
   b0 = b^2/(4 c) - a + \[Delta]b;
   
   BRescaled = Transpose[{-Sign[\[Alpha]] (ptsB[[All, 1]] - x0)/Lx, (residuals - \[Delta]b)/Lx}];
   residuals = ptsH[[All, 2]] - (linFit /. x -> ptsH[[All, 1]]);
   HRescaled = Transpose[{-Sign[\[Alpha]] (ptsH[[All, 1]] - x0)/Lx, (residuals - \[Delta]b)/Lx}];
   
   qF[x_] := (x^2 - 1);
   
   <|
     "BRescaled" -> BRescaled,
     "HRescaled" -> HRescaled,
     "qF" -> qF,
     "b0/Lx" -> b0/Lx,
     "\[Alpha]" -> \[Alpha],
     "\[Alpha]\[Epsilon]"-> \[Alpha]Error,
     "Lx"-> Lx
   |>
];

processWorkbook[file_String] := Module[{data},
   data = Import[file];
   Map[processSheet, data]
];

(* Function to calculate landslide thickness from model parameters*)
HSolution[\[Alpha]v_?NumericQ, fv_?NumericQ, b0_?NumericQ, k_?NumericQ, {xmin_, xmax_}] :=
  Module[{sol},
    sol = solveH[\[Alpha]v, fv, b0, k, {xmin, xmax}];
    If[sol === $Failed, Missing["NoSolution"], H /. sol]
  ];

(* === Root finding functions === *)
xPosRoot[\[Alpha]v_?NumericQ, fv_?NumericQ, b0_?NumericQ, k_?NumericQ] :=
  Module[{sol, res},
    sol = solveH[\[Alpha]v, fv, b0, k, {0, 1}];
    If[sol === $Failed, 1.0,
      Quiet @ Check[
        res = x /. Last[FindRoot[(H[x] /. sol) == 0, {x, 0.9, 0, 1}]],
        0
      ]
    ]
  ];

xNegRoot[\[Alpha]v_?NumericQ, fv_?NumericQ, b0_?NumericQ, k_?NumericQ] :=
  Module[{sol, res},
    sol = solveH[\[Alpha]v, fv, b0, k, {-1, 0}];
    If[sol === $Failed, -1.0,
      Quiet @ Check[
        res = x /. Last[FindRoot[(H[x] /. sol) == 0, {x, -0.9, -1, 0}]],
        0
      ]
    ]
  ];

(* === Stability locus functions for \[Alpha] in terms of f,b0,k === *)
locusPosFun[f_?NumericQ, b0_?NumericQ, k_?NumericQ] :=
  Module[{res},
    res = Quiet@FindRoot[
      xPosRoot[\[Alpha]sol, f, b0, k] == 0.99, {\[Alpha]sol, 0.0, 0, 10}];
    \[Alpha]sol /. res
  ];

locusNegFun[f_?NumericQ, b0_?NumericQ, k_?NumericQ] :=
  Module[{res},
    res = Quiet@FindRoot[
      xNegRoot[\[Alpha]sol, f, b0, k] == -0.99, {\[Alpha]sol, 1.2, 0, 10}];
    \[Alpha]sol /. res
  ];

End[];
EndPackage[];
