#!/usr/bin/python

import os.path
import re
import glob
import ROOT as root
from ROOT import gStyle

root.gROOT.SetBatch(True)
tools = root.TMVA.Tools.Instance();

def setStyle():
  gStyle.SetCanvasColor(0)
  gStyle.SetCanvasBorderSize(10)
  gStyle.SetCanvasBorderMode(0)
  gStyle.SetCanvasDefH(700)
  gStyle.SetCanvasDefW(700)

  gStyle.SetPadColor       (0)
  gStyle.SetPadBorderSize  (10)
  gStyle.SetPadBorderMode  (0)
  gStyle.SetPadBottomMargin(0.13)
  gStyle.SetPadTopMargin   (0.08)
  gStyle.SetPadLeftMargin  (0.15)
  gStyle.SetPadRightMargin (0.05)
  gStyle.SetPadGridX       (0)
  gStyle.SetPadGridY       (0)
  gStyle.SetPadTickX       (1)
  gStyle.SetPadTickY       (1)

  gStyle.SetFrameFillStyle ( 0)
  gStyle.SetFrameFillColor ( 0)
  gStyle.SetFrameLineColor ( 1)
  gStyle.SetFrameLineStyle ( 0)
  gStyle.SetFrameLineWidth ( 1)
  gStyle.SetFrameBorderSize(10)
  gStyle.SetFrameBorderMode( 0)

  gStyle.SetNdivisions(505)
 
  gStyle.SetLineWidth(2)
  gStyle.SetHistLineWidth(2)
  gStyle.SetFrameLineWidth(2)
  gStyle.SetLegendFillColor(root.kWhite)
  gStyle.SetLegendFont(42)
  gStyle.SetMarkerSize(1.2)
  gStyle.SetMarkerStyle(20)

  gStyle.SetLabelSize(0.040,"X")
  gStyle.SetLabelSize(0.040,"Y")

  gStyle.SetLabelOffset(0.010,"X")
  gStyle.SetLabelOffset(0.010,"Y")

  gStyle.SetLabelFont(42,"X")
  gStyle.SetLabelFont(42,"Y")

  gStyle.SetTitleBorderSize(0)
  gStyle.SetTitleFont(42)
  gStyle.SetTitleFont(42,"X")
  gStyle.SetTitleFont(42,"Y")

  gStyle.SetTitleSize(0.045,"X")
  gStyle.SetTitleSize(0.045,"Y")

  gStyle.SetTitleOffset(1.4,"X")
  gStyle.SetTitleOffset(1.4,"Y")

  gStyle.SetTextSize(0.055)
  gStyle.SetTextFont(42)

  gStyle.SetOptStat(0)

setStyle()

gStyle.SetPaintTextFormat("3g")
gStyle.SetPadBottomMargin   (0.20)
gStyle.SetPadLeftMargin  (0.22)
gStyle.SetPadRightMargin  (0.10)

def printSeparations(infilename):
  separations = {}
  fileNameMatch = re.match("TMVA_(.+)_(.+).root",infilename)
  infile = root.TFile(infilename)
  for fileKey in infile.GetListOfKeys():
    fileKeyName = fileKey.GetName()
    #match = re.match("Method_(.+)",fileKeyName)
    match = re.match("Method_BDT",fileKeyName)
    if match:
      folder1 = fileKey.ReadObj()
      for folder1Key in folder1.GetListOfKeys():
        folder2 = folder1Key.ReadObj()
        filePathName = folder1.GetName()+"/"+folder2.GetName()
        varNames = set()
        for folder2Key in folder2.GetListOfKeys():
          keyName =  folder2Key.GetName()
          matchS = re.match(r"(.+)__Signal",keyName)
          #matchB = re.match(r"(.+)__Background",keyName)
          if matchS:
            varNames.add(matchS.group(1))
        for varName in varNames:
          sHist = folder2.Get(varName+"__Signal")
          bHist = folder2.Get(varName+"__Background")
          sep = tools.GetSeparation(sHist,bHist)
          separations[varName] = sep
        sHist = folder2.Get("MVA_BDT_S")
        bHist = folder2.Get("MVA_BDT_B")
        sep = tools.GetSeparation(sHist,bHist)
        separations["BDT"] = sep
        
  fileTitle = "Unknown file!"
  if fileNameMatch:
    if fileNameMatch.group(1) == "vbf":
      fileTitle = "VBF"
    elif fileNameMatch.group(1) == "inclusive":
      fileTitle = "Non-VBF"
    energyStr = fileNameMatch.group(2)
    fileTitle += " "+energyStr.replace("TeV"," TeV")
  print(fileTitle)
  print("-"*40)
  for i in reversed(sorted(separations.keys(),key=lambda x: separations[x])):
    print("%-20s: %g" % (i,separations[i]))
  print
  canvas = root.TCanvas()
  fileTitleSave = fileTitle.replace(" ","")
  fileTitleSave = fileTitleSave.replace("-","")
  for i,iLong in zip(["S","B"],["Signal","Background"]):
    objName = "CorrelationMatrix"+i
    corMat = infile.Get(objName)
    corMat.SetTitle(fileTitle+" "+iLong)
    #corMat.Scale(0.01)
    corMat.Draw("coltext")
    canvas.SaveAs(fileTitle+"_"+objName+".png")
    canvas.SaveAs(fileTitle+"_"+objName+".pdf")
    #canvas.SaveAs(fileTitle+"_"+objName+".root")
    #canvas.SaveAs(fileTitle+"_"+objName+".eps")

if __name__ == "__main__":
  for i in glob.glob("TMVA*.root"):
    printSeparations(i)
