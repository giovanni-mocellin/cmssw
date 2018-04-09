import ROOT
import copy
import sys


def makeCanvas(name, doRatio):
    H_ref = 600;
    if doRatio:
        H_ref = 800;
    W_ref = 800;
    canvas = ROOT.TCanvas(name,name,50,50,W_ref,H_ref)    
    #canvas = ROOT.TCanvas(name,name,W_ref,H_ref)    
    return canvas


def setMargins(canvas, doRatio):
    H_ref = 600;
    if doRatio:
        H_ref = 800;
    W_ref = 800;
    W = W_ref
    H  = H_ref
    T = 0.08*H_ref
    B = 0.12*H_ref 
    L = 0.12*W_ref
    R = 0.04*W_ref
    canvas.SetFillColor(0)
    canvas.SetBorderMode(0)
    canvas.SetFrameFillStyle(0)
    canvas.SetFrameBorderMode(0)
    canvas.SetLeftMargin( L/W )
    canvas.SetRightMargin( R/W )
    canvas.SetTopMargin( T/H )
    canvas.SetBottomMargin( B/H )
    canvas.SetTickx(1)
    canvas.SetTicky(1)
    return canvas


#strName3D = "gemcrCf_g"
#if len(sys.argv) >= 3 and sys.argv[ 2 ] == "ALL": strName3D = "gemcr_g"
#if len(sys.argv) >= 3: strName3D = sys.argv[ 2 ]

strSrc1 = sys.argv[ 1 ]
strSrc2 = sys.argv[ 1 ]
#strSrc1 = "DQM_V0001_R017091404__Global__CMSSW_X_Y_Z__RECO.root"
#strSrc2 = "DQM_V0001_R017091901__Global__CMSSW_X_Y_Z__RECO.root"

strRunNum1 = strSrc1.split("_")[ 2 ].replace("R", "")
strRunNum2 = strSrc2.split("_")[ 2 ].replace("R", "")
nRunNum1 = int(strRunNum1)
nRunNum2 = int(strRunNum2)

strHistTmp1 = "DQMData/Run %i/MuonGEMRecHitsV/Run summary/GEMRecHitsTask/"%nRunNum1 + "chamber_%i_layer_%i_"
strHistTmp2 = "DQMData/Run %i/MuonGEMRecHitsV/Run summary/GEMRecHitsTask/"%nRunNum2 + "chamber_%i_layer_%i_"

fInput1 = ROOT.TFile(strSrc1)
fInput2 = ROOT.TFile(strSrc2)

nNoChFirst = 1
nLayerFirst = 1
strPrefix = ""

if len(sys.argv) >= 3: 
  nLayerFirst = 2 - ( int(sys.argv[ 2 ]) + 1 ) % 2
  nNoChFirst = int(sys.argv[ 2 ]) / 2 * 2 + 1
  strPrefix = "_%02i"%int(sys.argv[ 2 ])

hBkg1 = copy.deepcopy(fInput1.Get(strHistTmp1%(nNoChFirst, nLayerFirst) + "th2D_eff"))
hSig1 = copy.deepcopy(fInput1.Get(strHistTmp1%(nNoChFirst, nLayerFirst) + "recHit_efficiency"))

hBkg2 = copy.deepcopy(fInput2.Get(strHistTmp2%(nNoChFirst, nLayerFirst) + "th2D_eff_scint"))
hSig2 = copy.deepcopy(fInput2.Get(strHistTmp2%(nNoChFirst, nLayerFirst) + "recHit_efficiency_scint"))
#hBkg2 = copy.deepcopy(fInput2.Get(strHistTmp2%(nNoChFirst, nLayerFirst) + "th2D_eff"))
#hSig2 = copy.deepcopy(fInput2.Get(strHistTmp2%(nNoChFirst, nLayerFirst) + "recHit_efficiency"))

strHistCutflow1 = "DQMData/Run %i/MuonGEMRecHitsV/Run summary/GEMRecHitsTask/events_withtraj"%nRunNum1
strHistCutflow2 = "DQMData/Run %i/MuonGEMRecHitsV/Run summary/GEMRecHitsTask/events_withtraj"%nRunNum2

hCutflow1 = copy.deepcopy(fInput1.Get(strHistCutflow1))
hCutflow2 = copy.deepcopy(fInput2.Get(strHistCutflow2))

if len(sys.argv) < 3: 
  for i in range(30):
    if i == 0: continue
    
    if i % 10 == 0 or i % 10 == 1 or i % 10 == 8 or i % 10 == 9: continue
    
    nLayer = 2 - ( i + 1 ) % 2
    nNoCh = i / 2 * 2 + 1
    
    hBkgAdd1 = fInput1.Get(strHistTmp1%(nNoCh, nLayer) + "th2D_eff")
    hSigAdd1 = fInput1.Get(strHistTmp1%(nNoCh, nLayer) + "recHit_efficiency")

    hBkgAdd2 = fInput2.Get(strHistTmp2%(nNoCh, nLayer) + "th2D_eff_scint")
    hSigAdd2 = fInput2.Get(strHistTmp2%(nNoCh, nLayer) + "recHit_efficiency_scint")
    #hBkgAdd2 = fInput2.Get(strHistTmp2%(nNoCh, nLayer) + "th2D_eff")
    #hSigAdd2 = fInput2.Get(strHistTmp2%(nNoCh, nLayer) + "recHit_efficiency")
    
    hBkg1.Add(hBkgAdd1)
    hSig1.Add(hSigAdd1)

    hBkg2.Add(hBkgAdd2)
    hSig2.Add(hSigAdd2)

hBkg1D1 = ROOT.TH1F("hBkg1D1", "", 8, 0 - 0.5, 8 - 0.5)
hSig1D1 = ROOT.TH1F("hSig1D1", "", 8, 0 - 0.5, 8 - 0.5)

hBkg1D2 = ROOT.TH1F("hBkg1DD2", "", 8, 0 - 0.5, 8 - 0.5)
hSig1D2 = ROOT.TH1F("hSig1DD2", "", 8, 0 - 0.5, 8 - 0.5)

nNumBkg1 = 0
nNumSig1 = 0
nNumBkg2 = 0
nNumSig2 = 0

"""
for i in range(24): 
  nCol = i % 8 + 1
  nRow = i / 8 + 1
  
  hBkg1D1.SetBinContent(i + 1, hBkg1.GetBinContent(nRow, nCol))
  hSig1D1.SetBinContent(i + 1, hSig1.GetBinContent(nRow, nCol))
  
  hBkg1D2.SetBinContent(i + 1, hBkg2.GetBinContent(nRow, nCol))
  hSig1D2.SetBinContent(i + 1, hSig2.GetBinContent(nRow, nCol))
"""
for i in range(8): 
  nCol = i + 1
  
  hBkg1D1.SetBinContent(i + 1, hBkg1.GetBinContent(1, nCol) + hBkg1.GetBinContent(2, nCol) + hBkg1.GetBinContent(3, nCol))
  hSig1D1.SetBinContent(i + 1, hSig1.GetBinContent(1, nCol) + hSig1.GetBinContent(2, nCol) + hSig1.GetBinContent(3, nCol))
  
  hBkg1D2.SetBinContent(i + 1, hBkg2.GetBinContent(1, nCol) + hBkg2.GetBinContent(2, nCol) + hBkg2.GetBinContent(3, nCol))
  hSig1D2.SetBinContent(i + 1, hSig2.GetBinContent(1, nCol) + hSig2.GetBinContent(2, nCol) + hSig2.GetBinContent(3, nCol))
  
  nNumBkg1 += hBkg1D1.GetBinContent(i + 1)
  nNumSig1 += hSig1D1.GetBinContent(i + 1)
  
  nNumBkg2 += hBkg1D2.GetBinContent(i + 1)
  nNumSig2 += hSig1D2.GetBinContent(i + 1)

nNumCF1 = hCutflow1.GetNbinsX() / 2
nNumCF2 = hCutflow2.GetNbinsX() / 2

dHitL4_1  = hCutflow1.Integral(1, nNumCF1)
dNumTraj1 = hCutflow1.Integral(nNumCF1 + 1, 2 * nNumCF1)

dHitL4_2  = hCutflow2.Integral(1, nNumCF2)
dNumTraj2 = hCutflow2.Integral(nNumCF2 + 1, 2 * nNumCF2)

print "num events - hit >= 4 / reconstructed trajs :", dHitL4_1, "/", dNumTraj1, "(%lf)"%(dNumTraj1 / dHitL4_1)
print "num (bkg1) :", nNumBkg1
print "num (sig1) :", nNumSig1
print "num events - hit >= 4 / reconstructed trajs :", dHitL4_2, "/", dNumTraj2, "(%lf)"%(dNumTraj2 / dHitL4_2)
print "num (bkg2) :", nNumBkg2
print "num (sig2) :", nNumSig2

hEff1 = ROOT.TEfficiency(hSig1D1, hBkg1D1)
hEff2 = ROOT.TEfficiency(hSig1D2, hBkg1D2)

hEff1.SetLineColor(1)
hEff1.SetMarkerColor(1)
hEff1.SetMarkerStyle(20)

hEff2.SetLineColor(4)
hEff2.SetMarkerColor(4)
hEff2.SetMarkerStyle(34)

#Set canvas
canv = makeCanvas("canv", False)
setMargins(canv, False)

hInit = ROOT.TH1F("hFirst", "Efficiency vs Eta partition", 8, 0 - 0.5, 8 - 0.5)

hInit.GetXaxis().SetTitle("Eta partition")
hInit.GetYaxis().SetTitle("Efficiency")
hInit.GetYaxis().SetTitleOffset(1.25)

hInit.SetMinimum(0.95)
hInit.SetMaximum(1.00)

hInit.SetStats(0)
hInit.Draw()

hEff1.Draw("same")
hEff2.Draw("same")

leg = ROOT.TLegend(0.36, 0.33, 0.59, 0.18)

leg.AddEntry(hEff1, "no scintillators", "pl")
leg.AddEntry(hEff2, "with scintillators", "pl")
#leg.AddEntry(hEff1, "along y-axis", "pl")
#leg.AddEntry(hEff2, "random direction", "pl")

leg.SetTextFont(62)
leg.SetTextSize(0.03)
leg.SetBorderSize(0)
leg.Draw("same")

canv.Modified()
canv.Update()

canv.Draw()

#hEff.GetPaintedGraph().GetXaxis().SetTitle("Eta partition")
#hEff.GetPaintedGraph().GetYaxis().SetTitle("Efficiency")

#canv.SaveAs("ROCCurves_%s_%s_dPV%s.png"%(strPUTitle.replace(" ", ""), id, strDPV))
canv.SaveAs("Eff_%s%s.png"%(strRunNum1, strPrefix))



