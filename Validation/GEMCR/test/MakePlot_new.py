from ROOT import *
import os, sys, copy, shutil, csv


gROOT.SetBatch(1)


def getNameCh(ch, ly): return "chamber_%i_layer_%i_"%(ch, ly)


headLaTeX = """
\documentclass{beamer}
\usefonttheme[onlylarge]{structurebold}
\setbeamerfont*{frametitle}{size=\\normalsize,series=\\bfseries}
\setbeamertemplate{navigation symbols}{}
\setbeamertemplate{itemize items}[circle]
\setbeamerfont{page number in head/foot}{size=\small}
\setbeamertemplate{footline}[frame number]
\usepackage[english]{babel}
\usepackage[latin1]{inputenc}
\usepackage{times}
\usepackage{graphicx}
\usepackage[T1]{fontenc}
\usepackage{alltt}

\\newcommand{\\baseLoc}{./}

\\newcommand{\imageOne}[1]{
\scalebox{0.3}{
\includegraphics{\\baseLoc#1}
}
}

\\newcommand{\imageTwo}[2]{
\scalebox{0.26}{
\includegraphics{\\baseLoc#1}
\includegraphics{\\baseLoc#2}
}
}
\\newcommand{\imageThree}[3]{
\scalebox{0.18}{
\includegraphics{\\baseLoc#1}
\includegraphics{\\baseLoc#2}
\includegraphics{\\baseLoc#3}
}
}

\\newcommand{\imageFour}[4]{
\scalebox{0.18}{
\includegraphics{\\baseLoc#1}
\includegraphics{\\baseLoc#2}
}
\\
\scalebox{0.18}{
\includegraphics{\\baseLoc#3}
\includegraphics{\\baseLoc#4}
}
}

\\newcommand{\imageFive}[5]{
\scalebox{0.18}{
\includegraphics{\\baseLoc#1}
\includegraphics{\\baseLoc#2}
\includegraphics{\\baseLoc#3}
}
\\
\scalebox{0.18}{
\includegraphics{\\baseLoc#4}
\includegraphics{\\baseLoc#5}
}
}


\\newcommand{\imageSix}[6]{
\scalebox{0.18}{
\includegraphics{\\baseLoc#1}
\includegraphics{\\baseLoc#2}
\includegraphics{\\baseLoc#3}
}
\\
\scalebox{0.18}{
\includegraphics{\\baseLoc#4}
\includegraphics{\\baseLoc#5}
\includegraphics{\\baseLoc#6}
}
}

\\newcommand{\imageEight}[8]{
\scalebox{0.14}{
\includegraphics{\\baseLoc#1}
\includegraphics{\\baseLoc#2}
\includegraphics{\\baseLoc#3}
\includegraphics{\\baseLoc#4}
}
\\
\scalebox{0.14}{
\includegraphics{\\baseLoc#5}
\includegraphics{\\baseLoc#6}
\includegraphics{\\baseLoc#7}
\includegraphics{\\baseLoc#8}
}
}

\\newcommand{\imageFiveInOneLine}[5]{
\scalebox{0.11}{
\includegraphics{\\baseLoc#1}
\includegraphics{\\baseLoc#2}
\includegraphics{\\baseLoc#3}
\includegraphics{\\baseLoc#4}
\includegraphics{\\baseLoc#5}
}
}


\\begin{document}
"""


strFilename = sys.argv[ 1 ]
nRunNum = int(strFilename.split("_")[ 2 ].replace("R", ""))
strMainDir = "DQMData/Run %d/MuonGEMRecHitsV/Run summary/GEMRecHitsTask"%(nRunNum)

fDQM = TFile(strFilename)
dirMain = fDQM.Get(strMainDir)

strSaveDir = "Plots_%09i"%(nRunNum)
if not os.path.exists(strSaveDir): os.makedirs(strSaveDir)
os.chdir(strSaveDir)

c1 = TCanvas("tmp", "", 600, 600)

strLaTeX = headLaTeX


def saveHisto(nCh, nLy, strName, option = ""): 
  #strNameFront = strMainDir + "/chamber_%i_layer_%i_"%(nCh, nLy)
  strNameFront = "chamber_%i_layer_%i_"%(nCh, nLy)
  h1 = copy.deepcopy(dirMain.Get(strNameFront + strName))
  
  strSubDir = "%i_%i"%(nCh, nLy)
  if not os.path.exists(strSubDir): os.makedirs(strSubDir)
  
  strFilename = strSubDir + "/" + strSubDir + "_" + strName + ".png"
  h1.Draw(option)
  c1.SaveAs(strFilename)
  
  return strFilename


def drawEff(nCh, nLy, option = "", xtitle = ""): 
  hInit = TH1F("hFirst%i%i"%(nCh, nLy), "Efficiency vs VFAT, Chamber %i Layer %i"%(nCh, nLy), 24, 0 - 0.5, 24 - 0.5)
  
  hInit.GetXaxis().SetTitle(xtitle)
  hInit.GetYaxis().SetTitle("Efficiency")
  hInit.GetYaxis().SetTitleOffset(1.25)
  
  hInit.SetMinimum(0.90)
  hInit.SetMaximum(1.00)
  
  hInit.SetStats(0)
  hInit.Draw()
  
  strNameFront = "chamber_%i_layer_%i_"%(nCh, nLy)
  strNameDen = "th2D_eff"
  strNameNum = "recHit_efficiency"
  
  h2Den = copy.deepcopy(dirMain.Get(strNameFront + strNameDen))
  h2Num = copy.deepcopy(dirMain.Get(strNameFront + strNameNum))
  
  h1Den = TH1F("hDen", "", 24, 0 - 0.5, 24 - 0.5)
  h1Num = TH1F("hNum", "", 24, 0 - 0.5, 24 - 0.5)
  
  for i in range(24):
    nCol = i % 8 + 1
    nRow = i / 8 + 1
    
    h1Den.SetBinContent(i + 1, h2Den.GetBinContent(nRow, nCol))
    h1Num.SetBinContent(i + 1, h2Num.GetBinContent(nRow, nCol))
  
  h1 = TEfficiency(h1Num, h1Den)
  
  strSubDir = "%i_%i"%(nCh, nLy)
  if not os.path.exists(strSubDir): os.makedirs(strSubDir)
  
  strFilename = strSubDir + "/" + strSubDir + "_efficiency_plot.png"
  h1.Draw(option + "same")
  
  line98 = TLine(0 - 0.5, 0.98, 24 - 0.5, 0.98)
  line98.SetLineColor(2)
  line98.Draw()
  
  c1.Modified()
  c1.Update()
  c1.SaveAs(strFilename)
  
  return strFilename


arrBlacklist = []

for i in range(1, 31):
  nCh = i
  nLy = 1
  
  if i % 2 == 0:
    nCh = i - 1
    nLy = 2
  
  if [nCh, nLy] in arrBlacklist: continue
  
  strEffi = drawEff(nCh, nLy, xtitle = "VFAT")
  strRecHit = saveHisto(nCh, nLy, "recHit_efficiency", "colz text")
  strTH2D = saveHisto(nCh, nLy, "th2D_eff", "colz text")
  strRecHitStrip = saveHisto(nCh, nLy, "recHit", "colz")
  strGemDigi = saveHisto(nCh, nLy, "recHit_gemDigi", "colz")
  strHitMul = saveHisto(nCh, nLy, "hit_mul")
  
  strLaTeX += """
\\begin{frame}[plain]{chamber %i layer %i}
\\imageSix{%s}{%s}{%s}{%s}{%s}{%s}
\\end{frame}
  """%(nCh, nLy, strEffi, strRecHit, strTH2D, strRecHitStrip, strGemDigi, strHitMul)
  
strLaTeX += "\n\\end{document}"
open("res_%09i.tex"%nRunNum, "w").write(strLaTeX)
os.system("pdflatex res_%09i.tex"%nRunNum)
