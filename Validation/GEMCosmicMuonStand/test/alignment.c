{

cout<<"start"<<endl;
TString run  = "332";
TString run0 = "331";
int printOn=1;
int calcDx2=1;
int calcRz2=1;
int calcAlignmentFromTree=0;
int fixedSeedOn=0;
int viewOnlyOneColumn = 0;
int iColumn = 0;
int viewEtaPlots = 1;
TString trueVar = "tDx";
int viewResidualAllCh = 0;
TString v = "x";
TString V = "X";
if(trueVar=="tDz") v = "z";
if(trueVar=="tDz") V = "Z";

TString srcDir = "/afs/cern.ch/work/j/jslee/QC8";
TString savePngDir = "temp";
TString filename = "temp_out_reco_MC_"+run;
TFile *f = new TFile(srcDir+"/"+filename+".root");

FILE *inFile = NULL;
//inFile = fopen("table_input_run"+run+".txt","r");
inFile = fopen("table/table_output_run"+run0+".txt","r");
//inFile = fopen("table_output_run"+run0+".txt","r");

const int maxSch = 15; //number of Sch
const int maxNchamber = 30;

float dxIn[maxSch] = {0,};
float rzIn[maxSch] = {0,};

if(inFile == NULL) {
        cout<<"Please, set a input file!!!"<<endl;
}
else {
        float num[maxNchamber];
        float inNum;
	for(int i=0;i<maxNchamber;i++) if(!feof(inFile)) {
		fscanf(inFile, "%f", &inNum);
		num[i] = inNum;

		if(i<maxSch) {
			int j = (i%3+1)*5-1 - i/3;
			dxIn[j] = num[i];
			//cout<<"i "<<i<<", j "<<j<<", dxIn "<<dxIn[j]<<endl;
		}
		else {
			int j = (i%3+1)*5-1 - (i-maxSch)/3;
			rzIn[j] = num[i];
			//cout<<"i "<<i<<", j "<<j<<", rzIn "<<rzIn[j]<<endl;
		}
	}
/*
	cout<<endl;
	for(int i=0;i<maxSch;i++) {
		cout<<"dxIn "<<dxIn[i]<<", rzIn "<<rzIn[i]<<endl;
	}
	cout<<endl;
	for(int i=0;i<maxSch;i++) {
		int j = (i%3+1)*5-1 - i/3;
		cout<<"i "<<i<<", j "<<j<<", dxIn "<<dxIn[j]<<", rzIn "<<rzIn[j]<<endl;
	}
	cout<<endl;
	for(int i=0;i<maxSch;i++) {
		int j = (i%3+1)*5-1 - i/3;
		cout<<dxIn[j]<<" ";
		if(i%3==2) cout<<endl;
	}
	cout<<endl;
	for(int i=0;i<maxSch;i++) {
		int j = (i%3+1)*5-1 - i/3;
		cout<<rzIn[j]<<" ";
		if(i%3==2) cout<<endl;
	}
*/
}
fclose(inFile);
	

//TString srcDir = "/afs/cern.ch/work/g/gmocelli/PhD/alignment/src/Validation/GEMCosmicMuonStand/test/res_"+run;
//TString filename = "out";
//TFile *f = new TFile(srcDir+"/"+filename+".root");

//f->cd("gemcrValidation");
f->cd("AlignmentValidationQC8");
cout<<"total triggered ev : "<<hev->GetBinContent(1)<<endl;
cout<<"total good ev      : "<<tree->GetEntries()<<endl;

TH1F *hResidualX[maxNchamber];
TH1F *hResidualXE[maxNchamber];
double residualXMean[maxNchamber]={0,};
double residualXMeanE[maxNchamber]={0,};
double residualXSigma[maxNchamber]={0,};
double residualXSigmaE[maxNchamber]={0,};

double nTrajHit[maxNchamber]={0,};
double nTrajHitRatio[maxNchamber]={0,};
double nTrajHitRatioE[maxNchamber]={0,};
TH1F *heff[maxNchamber];
double eff[maxNchamber]={0,};
double effE[maxNchamber]={0,};

TGraphErrors *gavgChi2[maxNchamber];
TGraphErrors *gnTrajHitRatio[maxNchamber];
TGraphErrors *gResidualXMean[maxNchamber];
TGraphErrors *geff[maxNchamber];

TF1 *fPol1 = new TF1("fPol1","[0]*x+[1]");
TF1 *fPol2 = new TF1("fPol2","[0]*(x-[1])*(x-[1])+[2]");
TF1 *fPol1R = new TF1("fPol1","x/[0]+[1]");
double xParamMax = 1000;
fPol1R->SetParLimits(0,-0.1,0.1);
fPol1R->SetParLimits(1,-xParamMax,xParamMax);
TF1 *ga1 = new TF1("ga1","([0]/([2]*2.50663))*exp(-0.5*(x-[1])*(x-[1])/([2]*[2]))");
ga1->SetParName(1,"mean");
TF1 *ga = new TF1("ga","([0]/([2]*2.50663))*exp(-0.5*(x-[1])*(x-[1])/([2]*[2])) + ([3]/([4]*2.50663))*exp(-0.5*(x-[1])*(x-[1])/([4]*[4]))");
ga->SetParName(1,"mean");
ga->SetParLimits(0,1,100000);
ga->SetParLimits(2,0,10);
ga->SetParLimits(3,1,100000);
ga->SetParLimits(4,0,10);
TF1 *ga3 = new TF1("ga3","([0]/([2]*2.50663))*exp(-0.5*(x-[1])*(x-[1])/([2]*[2])) + ([3]/([4]*2.50663))*exp(-0.5*(x-[1])*(x-[1])/([4]*[4])) + ([5]/([6]*2.50663))*exp(-0.5*(x-[1])*(x-[1])/([6]*[6]))");
ga3->SetParName(1,"mean");
ga3->SetParameters(100,0,0.01,100,0.1,100,0.5);
ga3->SetParLimits(0,10,10000);
ga3->SetParLimits(2,0.001,1);
ga3->SetParLimits(3,10,10000);
ga3->SetParLimits(4,0.01,1);
ga3->SetParLimits(5,10,1000);
ga3->SetParLimits(6,0.00001,0.1);

TF1 *Breit_Wigner = new TF1("Breit_Wigner","[0]/((x-[1])*(x-[1])+[2]*[2]/4)");
Breit_Wigner->SetParName(1,"mean");
Breit_Wigner->SetParName(2,"sigma");

TH2F *hResidualX2D;

TString hname;
TString htitle;
TString png;
TString png0;
TString var;
TString cut;

TString TB[2] = {"B","T"};

const int NSch = 5;
//const int NSch = 15;
const int Nrow = 10;
double row[Nrow] = {22.66,26.13,45.66,49.13,68.66,72.13,91.66,95.13,114.66,118.13}; // z position of recHits
double rowE[Nrow] = {0,};

double z[NSch];
double zE[NSch] = {0,};
for(int i=0;i<NSch;i++) z[i] = (row[i*2] + row[i*2+1])/2.;
	
TGraphErrors *gResidualXMeanVSch;

const int Nrow2 = 8;
double row2[NSch][Nrow] = {0,};
double rowE2[NSch][Nrow] = {0,};
//double row2[NSch][maxNchamber] = {0,};
//double rowE2[NSch][maxNchamber] = {0,};
double bestNchi2 = {0,};

TGraphErrors *gResidualXMeanVSch2[NSch];

const int Nrow3 = 6;
double row3[NSch][NSch][Nrow] = {0,};
double rowE3[NSch][NSch][Nrow] = {0,};
double bestNchi2_2 = {0,};

TGraphErrors *gResidualXMeanVSch3[NSch][NSch];

TCanvas *c1 = new TCanvas("c1","c1",700,500);
TCanvas *c2 = new TCanvas("c2","c2",700,500);
c2->SetRightMargin(0.15);

TCanvas *c5[maxNchamber];

const int nColumn = 3;
TH1F *htDx[maxNchamber];
TH1F *hsDx[maxNchamber];
double tDx[nColumn][NSch] = {0,};
double sDx[nColumn][NSch] = {0,};
TH1F *htRz[maxNchamber];
TH1F *hsRz[maxNchamber];
double tRz[nColumn][NSch] = {0,};
double sRz[nColumn][NSch] = {0,};

for(int it=0;it<nColumn;it++) 
	for(int is=0;is<NSch;is++)
	{
		int j = it*5 + is;
		sDx[it][is] = dxIn[j];
		sRz[it][is] = rzIn[j];
	}

TH1F *hZenithAngle[nColumn];
double zenithAngle[nColumn] = {0,};
double zenithAngleE[nColumn] = {0,};

const int Neta = 8;
const int nY = 9;
////double etaY[nY] = {65,40,20,5,-10,-30,-37.5,-55,-65};
////double recY[Neta] = {50.7,31.3,13.5,-2.67,-17.5,-31.0,-43.4,-54.8};
//double etaY[nY] = {87.5,62.5,45,27.5,12.5,-2.5,-15,-27.5,-45};
//double recY[Neta] = {72.7,53.5,35.5,19.3,4.50,-9.01,-21.4,-32.8};
//double recYE[Neta] = {0,};
////double etaYcenter = (etaY[0]+etaY[Neta-1])/2.;
////double recYcenter = (recY[0]-recY[Neta-1])/2;
//double recYcenter = (recY[0]+recY[Neta-1])/2.;

double nMergedFile = hnMergedFile->GetBinContent(1);
hEtaY->Scale(1./nMergedFile);
hEtaYRange->Scale(1./nMergedFile);

double etaY[nColumn][nY] = {0,};
double recY[nColumn][Neta] = {0,};
double recYE[nColumn][Neta] = {0,};
double recYcenter[nColumn] = {0,};
for(int i=0;i<nColumn;i++) 
{
	for(int iy=0;iy<nY;iy++)
	{
		int ibin = i*10*10 + iy + 1;
		double recYtemp = 1000;
		if(iy<Neta) 
		{
			recYtemp = hEtaY->GetBinContent(ibin);
			recY[i][iy] = hEtaY->GetBinContent(ibin);
		}
		etaY[i][iy] = hEtaYRange->GetBinContent(ibin);
		cout<<"i "<<i<<", iy "<<iy<<", etaY "<<etaY[i][iy]<<", recY "<<recYtemp<<endl;
	}
	//recYcenter[ic] = (recY[ic][0]+recY[ic][Neta-1])/2.;
	recYcenter[i] = (etaY[i][0]+etaY[i][nY-1])/2.;
}

TH1F *hZenithAngleVSieta[nColumn][Neta];
double zenithAngleVSieta[nColumn][Neta] = {0,};
double zenithAngleEVSieta[nColumn][Neta] = {0,};

TGraphErrors *gr_zenithAngleVSieta[nColumn];

double rzVSieta[nColumn][NSch][Neta] = {0,};

TH1D *hchEtaResidualX[maxNchamber][Neta];

double resMean[maxNchamber][Neta];
double resMeanE[maxNchamber][Neta];

double rz1_ch[maxNchamber] = {0,};
double dx1_ch[maxNchamber] = {0,};

double rz1[nColumn][NSch], rz1E[nColumn][NSch];
double rz2[nColumn][NSch], rz2E[nColumn][NSch];
double dx11[nColumn][NSch], dx11E[nColumn][NSch];

TGraphErrors *gr_rz[maxNchamber];

TF1 *ga2 = new TF1("ga2","([0]/([2]*2.50663))*exp(-0.5*(x-[1])*(x-[1])/([2]*[2])) + ([3]/([5]*2.50663))*exp(-0.5*(x-[4])*(x-[4])/([5]*[5]))");
ga2->SetParameters(300,0,0.1,300,0,0.5);
ga2->SetParName(1,"mean");
ga2->SetParName(4,"mean2");
ga2->SetParLimits(0,0,100000);
ga2->SetParLimits(1,-1.5,1.5);
ga2->SetParLimits(2,0.001,0.5);
ga2->SetParLimits(3,0,100000);
ga2->SetParLimits(4,-1.5,1.5);
ga2->SetParLimits(5,0.01,2);

TCanvas *c11[maxNchamber][Neta];
TCanvas *c02[maxNchamber];

int emptySch[nColumn][NSch] = {0,};
int emptyColumn[nColumn] = {0,};
for(int i=0;i<nColumn;i++) emptyColumn[i] = 1;

c1->cd();

for(int ic=0;ic<maxNchamber;ic++)
//for(int ic=20;ic<22;ic++)
//int ic=21;
if( (viewOnlyOneColumn&&ic/10==iColumn) || !viewOnlyOneColumn )
{
	if(viewResidualAllCh) c5[ic] = new TCanvas(Form("c5_%d",ic),Form("c5_%d",ic),700,500);
	TString chTitle = Form("%d/%d/",(ic/2)%5+1,ic/10+1) + TB[ic%2];

	hname = Form("htDx_%d",ic);
	htDx[ic] = new TH1F(hname,"",5000,-5,5);
	var = Form("tDx[%d]",ic);
	cut="";
	tree->Project(hname,var,cut);
	tDx[ic/10][(ic%10)/2] = htDx[ic]->GetMean();

	//hname = Form("hsDx_%d",ic);
	//hsDx[ic] = new TH1F(hname,"",5000,-5,5);
	//var = Form("dx[%d]",ic);
	//cut="";
	//tree->Project(hname,var,cut);
	//sDx[ic/10][(ic%10)/2] = hsDx[ic]->GetMean();

	hname = Form("htRz_%d",ic);
	htRz[ic] = new TH1F(hname,"",5000,-5,5);
	var = Form("tRz[%d]",ic);
	cut="";
	tree->Project(hname,var,cut);
	tRz[ic/10][(ic%10)/2] = htRz[ic]->GetMean();

	//hname = Form("hsRz_%d",ic);
	//hsRz[ic] = new TH1F(hname,"",5000,-5,5);
	//var = Form("rz[%d]",ic);
	//cut="";
	//tree->Project(hname,var,cut);
	//sRz[ic/10][(ic%10)/2] = hsRz[ic]->GetMean();

	//double nev0 = tree->GetEntries(trueCut);
	double nev0 = tree->GetEntries();
	cut = Form("chI[%d]",ic);
	nTrajHit[ic] = tree->GetEntries(cut);
	if(nev0>0) nTrajHitRatio[ic] = nTrajHit[ic]/nev0;
	if(nev0>0) nTrajHitRatioE[ic] = sqrt(nTrajHitRatio[ic]*(1-nTrajHitRatio[ic]))/nev0;

	hname = Form("hResidualX_%d",ic);
	//htitle = "Residual"+V+" (trajHit"+V+"-recHit"+V+", d"+v+", "+chTitle+")";
	htitle = "Residual"+V+" (trajHit"+V+"-recHit"+V+", "+chTitle+")";
	hResidualX[ic] = new TH1F(hname,htitle,100,-2,2);
	//if(trueVar=="tDx") hResidualX[ic] = new TH1F(hname,htitle,100,-2,2);
	//if(trueVar=="tDz") hResidualX[ic] = new TH1F(hname,htitle,60,-10,10);
	if(trueVar=="tDx") var = Form("chTrajHitX[%d]-chRecHitX[%d]",ic,ic);
	//if(trueVar=="tDx") var = Form("chTrajHitX[%d]-chRecHitX[%d]",ic,ic);
	////if(trueVar=="tDx") var = Form("chCHitX[%d]-chRecHitX[%d]",ic,ic);
	if(trueVar=="tDz") var = Form("chPHitZ[%d]-chRecHitZ[%d]",ic,ic);
	cut += Form("&&chF[%d]",ic);
	//cout<<"it "<<it<<", ic "<<ic<<", cut "<<cut<<endl;

	double numerator = tree->GetEntries(cut);
	double denominator = nTrajHit[ic];
	cout<<"ic "<<ic<<", numerator "<<numerator<<", denominator "<<denominator<<endl;

	bool emptyData = false;
	if(numerator!=0 || denominator!=0) emptyColumn[ic/10] = 0;
	if(numerator==0 && denominator==0) 
	{
		emptySch[ic/10][(ic%10)/2] = 1;
		tDx[ic/10][(ic%10)/2] = 0;
		tRz[ic/10][(ic%10)/2] = 0;
		continue;
	}

	if(denominator>0) eff[ic] = numerator/denominator;
	if(denominator>0) effE[ic] = sqrt(eff[ic]*(1-eff[ic]))/denominator;

	//if(ic/10==0) cut += Form("&&chI[%d]&&chI[%d]",9,0);
	//if(ic/10==1) cut += Form("&&chI[%d]&&chI[%d]",19,10);
	//if(ic/10==2) cut += Form("&&chI[%d]&&chI[%d]",29,20);
	cut += Form("&&chI[%d]&&chI[%d]",(ic/10)*10,(ic/10)*10+9);
	tree->Project(hname,var,cut);

	//ga->SetParLimits(1,-0.55,0.55);
	double area = hResidualX[ic]->GetEntries() * hResidualX[ic]->GetBinWidth(1);
	residualXMean[ic] = hResidualX[ic]->GetMean();
	residualXMeanE[ic] = hResidualX[ic]->GetRMS();
	ga->SetParameters(area, residualXMean[ic], residualXMeanE[ic]/4., area, residualXMeanE[ic]);
	
	hResidualX[ic]->Draw("e");
	hResidualX[ic]->Fit("ga");
	residualXMean[ic] = ga->GetParameter(1);
	residualXMeanE[ic] = ga->GetParError(1);
	double sigma = ga->GetParameter(2);
	//ga->SetParLimits(1,residualXMean[it][ic]-sigma, residualXMean[it][ic]+sigma);

	double Nchi2 = 100;
	Nchi2 = ga->GetChisquare()/ga->GetNDF();
	cout<<"ic "<<ic<<", nev "<<nTrajHit[ic]<<", nev0 "<<nev0<<", nTR "<<nTrajHitRatio[ic]<<" resXm "<<residualXMean[ic]<<" +- "<<residualXMeanE[ic]<<", Nchi2 "<<Nchi2<<endl;

	if(Nchi2>20) residualXMean[ic] = hResidualX[ic]->GetMean();
	if(Nchi2>20) residualXMeanE[ic] = hResidualX[ic]->GetMeanError();

	//hResidualX2D->SetBinContent(ic/10+1,ic%10+1,residualXMean[ic]);
	//hResidualX2D->SetBinError(ic/10+1,ic%10+1,residualXMeanE[ic]);

	cout<<"ic "<<ic<<", nev "<<nTrajHit[ic]<<", nev0 "<<nev0<<", nTR "<<nTrajHitRatio[ic]<<" resXm "<<residualXMean[ic]<<" +- "<<residualXMeanE[ic]<<", Nchi2 "<<Nchi2<<", resXs "<<residualXSigma[ic]<<", eff "<<eff[ic]<<endl;

	//hChi2[ic][ix][0][0]->Draw();
	//png = savePngDir+"/run"+run+Form("_NChi2_ch%d_ix%d.png",ic,ix);
	//if(printOn) c1->Print(png);

	c1->cd();
	hResidualX[ic]->Draw("e");
	png = savePngDir+"/run"+run+"_residual"+V+Form("_ch%d.png",ic);
	if(printOn) c1->Print(png);

	for(int ie=0;ie<Neta;ie++)
	{
		c1->cd();
                hname = Form("c11_%d_%d",ic,ie);
                if(viewEtaPlots) c11[ic][ie] = new TCanvas(hname,hname,700,500);
                hname = Form("hchEtaResidualX_%d_%d",ic,ie);
		htitle = "Residual"+V+" (trajHit"+V+"-recHit"+V+", "+chTitle+Form(", i#eta%d)",ie+1);
                if(calcAlignmentFromTree) hchEtaResidualX[ic][ie] = new TH1D(hname,htitle,200,-3,3);
                //if(calcAlignmentFromTree && (calcDx2 || calcRz2)) hchEtaResidualX[ic][ie] = new TH1D(hname,htitle,60,-1.5,1.5);
                var = Form("chTrajHitX[%d]-chRecHitX[%d]",ic,ic);
		int downCh = (ic/10)*10;
		int upCh = (ic/10)*10+9;
                //cut = Form("chI[%d]&&chI[%d]&&chF[%d]&&%f<chTrajHitY[%d]&&chTrajHitY[%d]<%f",downCh, upCh, ic, etaY[ie+1],ic,ic,etaY[ie]);
                cut = Form("chI[%d]&&chI[%d]&&chF[%d]&&%f<chTrajHitY[%d]&&chTrajHitY[%d]<%f",downCh, upCh, ic, etaY[ic/10][ie+1],ic,ic,etaY[ic/10][ie]);
		cut += "&&fabs(atan(trajPy/trajPz))<0.1";
		int downSch = downCh/2;
		int upSch = upCh/2;
		if(fixedSeedOn) cut += Form("&&SchSeed[0]==%d&&SchSeed[1]==%d",downSch,upSch);

		//if(calcDx2 || calcRz2)
		//{
			if(calcAlignmentFromTree) tree->Project(hname,var,cut);
			else hchEtaResidualX[ic][ie] = (TH1D*)f->Get("AlignmentValidationQC8/"+hname);
		//}
		if(hchEtaResidualX[ic][ie]->GetEntries()==0) 
		{
			cout<<"Empty data!!!!!!!!!!!"<<endl;
			emptyData = true;
			continue;
			//break;
		}
                //ga2->SetParameters(100,0,0.1,100,0,0.5);
                //hchEtaResidualX[ic][ie]->Fit(ga2);
                //resMean[ic][ie] = ga2->GetParameter(1);
                //resMeanE[ic][ie] = ga2->GetParError(1);
		double areaFit = hchEtaResidualX[ic][ie]->GetEntries() * hchEtaResidualX[ic][ie]->GetBinWidth(1);
		double meanFit = hchEtaResidualX[ic][ie]->GetMean();
		double sigmaFit = hchEtaResidualX[ic][ie]->GetRMS();
		ga1->SetParameters(areaFit, meanFit, sigmaFit);
                hchEtaResidualX[ic][ie]->Fit(ga1);
                if(viewEtaPlots) hchEtaResidualX[ic][ie]->Draw("e");
                resMean[ic][ie] = ga1->GetParameter(1);
                resMeanE[ic][ie] = ga1->GetParError(1);
                cout<<"ev "<<hchEtaResidualX[ic][ie]->GetEntries()<<", Nchi2 "<<ga1->GetChisquare()/ga1->GetNDF()<<", "<<var<<", "<<cut<<", mean "<<hchEtaResidualX[ic][ie]->GetMean()<<"(histo) "<<resMean[ic][ie]<<"(fit)"<<endl;
		png = savePngDir+"/run"+run+"_residual"+V+Form("_ch%d_ieta%d.png",ic,ie+1);
		if(printOn && viewEtaPlots) c11[ic][ie]->Print(png);
	}
	//if(emptyData)
	//{
	//	cout<<"Empty data!!!!!!!!!!!"<<endl;
	//	break;
	//}

	//double recY0[Neta] = {0,};
	//double recYE0[Neta] = {0,};
	//for(int ie=0;ie<Neta;ie++)
	//{
	//	recY0[ie] = recY[ic][ie];
	//	recYE0[ie] = recYE[ic][ie];
	//	cout<<"ic "<<ic<<", ie "<<ie<<", resMean "<<resMean[ic][ie]<<", recY "<<recY[ic][ie]<<endl;
	//}

        hname = Form("c02_%d",ic);
        c02[ic] = new TCanvas(hname,hname,700,500);
        //gr_rz[ic] = new TGraphErrors(Neta, recY, resMean[ic], recYE, resMeanE[ic]);
        //gr_rz[ic] = new TGraphErrors(Neta, resMean[ic], recY, resMeanE[ic], recYE);
        //gr_rz[ic] = new TGraphErrors(Neta, resMean[ic], recY[ic], resMeanE[ic], recYE[ic]);
        //gr_rz[ic] = new TGraphErrors(Neta, resMean[ic], recY0, resMeanE[ic], recYE0);
        gr_rz[ic] = new TGraphErrors(Neta, resMean[ic], recY[ic/10], resMeanE[ic], recYE[ic/10]);
        gr_rz[ic]->SetTitle("ResidualX mean VS i#eta");
        gr_rz[ic]->GetXaxis()->SetTitle("dx' [cm]");
        gr_rz[ic]->GetYaxis()->SetTitle("y(i#eta) [cm]");

        //double slope = (recY[Neta-1]-recY[0])/(resMean[ic][Neta-1]-resMean[ic][0]);
	//double slope = (resMean[ic][Neta-1]-resMean[ic][0])/(recY[Neta-1]-recY[0]);
	//double slope = (resMean[ic][Neta-1]-resMean[ic][0])/(recY[ic][Neta-1]-recY[ic][0]);
	double slope = (resMean[ic][Neta-1]-resMean[ic][0])/(recY[ic/10][Neta-1]-recY[ic/10][0]);
        //fPol1->SetParameter(0,slope);
        //gr_rz[ic]->Fit(fPol1);
        //double p0 = fPol1->GetParameter(0);
        //double p1 = fPol1->GetParameter(1);
        //rz1_ch[ic] = atan(-1./p0)*180./3.141592;
        //dx1_ch[ic] = -p1/p0;
        fPol1R->SetParameters(0.0001,0);
        gr_rz[ic]->Fit(fPol1R);
	Nchi2 = fPol1R->GetChisquare()/fPol1R->GetNDF();
        double p0 = fPol1R->GetParameter(0);
        double p1 = fPol1R->GetParameter(1);
        if(Nchi2>100 || fabs(p1)>=fabs(xParamMax)) 
	{
		//fPol1R->SetParameters(slope,0);
		//fPol1R->SetParameters(-p0,0);
		//y = x/slope + y0, y0 = y - x/slope
		int indexCenter = Neta/2;
		//double y0 = recY[indexCenter] - resMean[ic][indexCenter]/slope;
		//double y0 = recY[ic][indexCenter] - resMean[ic][indexCenter]/slope;
		double y0 = recY[ic/10][indexCenter] - resMean[ic][indexCenter]/slope;
		fPol1R->SetParameters(slope, y0);
	        gr_rz[ic]->Fit(fPol1R);
        	Nchi2 = fPol1R->GetChisquare()/fPol1R->GetNDF();
	        p0 = fPol1R->GetParameter(0);
        	p1 = fPol1R->GetParameter(1);
	}
        rz1_ch[ic] = atan(-p0)*180./3.141592;
        dx1_ch[ic] = -p1*p0;


	if(ic%2==1)
	{
		rz1[ic/10][(ic%10)/2]  = (rz1_ch[ic] + rz1_ch[ic-1])/2.;
		dx11[ic/10][(ic%10)/2] = (dx1_ch[ic] + dx1_ch[ic-1])/2.;
	}
        cout<<"rz'1 = "<<rz1_ch[ic]<<", dx'1 = "<<dx1_ch[ic]<<", Nchi2 "<<Nchi2<<endl;
	htitle = "ResidualX mean VS i#eta ("+chTitle+Form(", dx'_{1}=%1.3f, rz'_{1}=%1.3f)",dx1_ch[ic],rz1_ch[ic]);
	gr_rz[ic]->SetTitle(htitle);
        gr_rz[ic]->Draw("ap");
	png = savePngDir+"/run"+run+"_residual"+V+Form("_VS_ieta_ch%d.png",ic);
	if(printOn) c02[ic]->Print(png);
}


c2->cd();
for(int i=0;i<nColumn;i++) if( (viewOnlyOneColumn&&i==iColumn) || !viewOnlyOneColumn ) if(calcDx2 || calcRz2) if(!emptyColumn[i])
{
	hname = Form("hZenithAngle_%d",i);
	htitle = Form("zenith angle (column%d)",i+1);
	hZenithAngle[i] = new TH1F(hname,htitle,40,-0.05,0.05);
	//hZenithAngle[i] = new TH1F(hname,htitle,40,-0.5,0.5);
	cut = Form("chI[%d]&&chI[%d]",i*10,i*10+9);
	tree->Project(hname,"atan(trajPx/trajPz)",cut);

	ga->SetParameters(20,0,0.01,20,0.1);
	ga->SetParLimits(1,-0.04,0.04);
	hZenithAngle[i]->Fit("ga");
	zenithAngle[i] = ga->GetParameter(1);
	zenithAngleE[i] = ga->GetParError(1);
	//ga3->SetParameters(100,0,0.01,100,0.1,100,0.5);
	//hZenithAngle[i]->Fit("ga3");
	//zenithAngle[i] = ga3->GetParameter(1);
	//zenithAngleE[i] = ga3->GetParError(1);
	hZenithAngle[i]->Draw("e");
	png = savePngDir+"/run"+run+Form("_zenith_angle_column%d_.png",i+1);
	if(printOn) c2->Print(png);

	for(int ie=0;ie<Neta;ie++)
	//int ie=0;
	{
		hname = Form("hZenithAngleVSieta_%d_%d",i,ie);
		htitle = Form("zenith angle (column%d, i#eta%d)",i+1,ie+1);
		hZenithAngleVSieta[i][ie] = new TH1F(hname,htitle,40,-0.05,0.05);
		//hZenithAngleVSieta[i][ie] = new TH1F(hname,htitle,40,-0.5,0.5);
                //cut = Form("chI[%d]&&chI[%d]&&%f<chTrajHitY[%d]&&chTrajHitY[%d]<%f&&%f<chTrajHitY[%d]&&chTrajHitY[%d]<%f",i*10,i*10+9, etaY[ie+1],i*10,i*10,etaY[ie], etaY[ie+1],i*10+9,i*10+9,etaY[ie]);
                cut = Form("chI[%d]&&chI[%d]&&%f<chTrajHitY[%d]&&chTrajHitY[%d]<%f&&%f<chTrajHitY[%d]&&chTrajHitY[%d]<%f",i*10,i*10+9, etaY[i][ie+1],i*10,i*10,etaY[i][ie], etaY[i][ie+1],i*10+9,i*10+9,etaY[i][ie]);
		cut += "&&fabs(atan(trajPy/trajPz))<0.1";
		tree->Project(hname,"atan(trajPx/trajPz)",cut);

		ga->SetParameters(20,0,0.01,20,0.1);
		hZenithAngleVSieta[i][ie]->Fit("ga");
		zenithAngleVSieta[i][ie] = ga->GetParameter(1);
		zenithAngleEVSieta[i][ie] = ga->GetParError(1);
		if(fabs(zenithAngleVSieta[i][ie])*10 < fabs(zenithAngleEVSieta[i][ie]))
			zenithAngleEVSieta[i][ie] = 0;
		//ga3->SetParameters(100,0,0.01,100,0.1,100,0.5);
		//hZenithAngleVSieta[i][ie]->Fit("ga3");
		//zenithAngleVSieta[i][ie] = ga3->GetParameter(1);
		//zenithAngleEVSieta[i][ie] = ga3->GetParError(1);
		hZenithAngleVSieta[i][ie]->Draw("e");
		png = savePngDir+"/run"+run+Form("_zenith_angle_column%d_eta%d.png",i+1,ie+1);
		if(printOn) c2->Print(png);

//		for(int is=0;is<NSch;is++)
//		{
//			//rzVSieta[i][is][ie] = atan((z[is]-z[2])*tan(zenithAngleVSieta[i][ie])/recY[ie])*180/3.141593;
//			//rzVSieta[i][is][ie] = atan((z[is]-z[2])*TMath::Tan(zenithAngleVSieta[i][ie])/recY[ie])*180/3.141593;
//			double zAngle = tan(zenithAngleVSieta[i][ie]);
//			//rzVSieta[i][is][ie] = atan((z[is]-z[2])*zAngle/recY[ie])*180/3.141593;
//			rzVSieta[i][is][ie] = atan((z[is]-z[2])*zAngle/(recY[ie]-recYcenter))*180/3.141593;
//			if(ie==0) rz2[i][is] = rzVSieta[i][is][ie];
//		}
	}

	//gr_zenithAngleVSieta[i] = new TGraphErrors(Neta,zenithAngleVSieta[i],recY,zenithAngleEVSieta[i],recYE);
	gr_zenithAngleVSieta[i] = new TGraphErrors(Neta,zenithAngleVSieta[i],recY[i],zenithAngleEVSieta[i],recYE[i]);
	gr_zenithAngleVSieta[i]->Draw("ap");
	gr_zenithAngleVSieta[i]->Fit(fPol1R);

        //double slope = (recY[Neta-1]-recY[0])/(resMean[ic][Neta-1]-resMean[ic][0]);
	//double slope = (recY[Neta-1]-recY[0])/(zenithAngleVSieta[i][Neta-1]-zenithAngleVSieta[i][0]);
	double slope = (recY[i][Neta-1]-recY[i][0])/(zenithAngleVSieta[i][Neta-1]-zenithAngleVSieta[i][0]);
        //fPol1->SetParameter(0,slope);
        //gr_rz[ic]->Fit(fPol1);
        //double p0 = fPol1->GetParameter(0);
        //double p1 = fPol1->GetParameter(1);
        //rz1_ch[ic] = atan(-1./p0)*180./3.141592;
        //dx1_ch[ic] = -p1/p0;
        //fPol1R->SetParameters(0.0001,0);
        //gr_rz[ic]->Fit(fPol1R);
	double Nchi2 = fPol1R->GetChisquare()/fPol1R->GetNDF();
	cout<<"column"<<i+1<<" Nchi2 of zenith angle correlation plot : "<<Nchi2<<endl;
        if(Nchi2>20) 
	{
		fPol1R->SetParameters(slope,0);
 		gr_zenithAngleVSieta[i]->Fit(fPol1R);
        	Nchi2 = fPol1R->GetChisquare()/fPol1R->GetNDF();
		cout<<"column"<<i+1<<" Nchi2 of zenith angle correlation plot (new fit) : "<<Nchi2<<endl;
	}
	png = savePngDir+"/run"+run+Form("_zenith_angle_VS_ieta_of_column%d.png",i+1);
	if(printOn) c2->Print(png);
	double p0 = fPol1R->GetParameter(0);
	double p1 = fPol1R->GetParameter(1);
	double p0E = fPol1R->GetParError(0);
	double p1E = fPol1R->GetParError(1);
	//zenithAngle[i] = p0*(etaYcenter-p1); // zenith angle of i-column
	//zenithAngle[i] = p0*(recYcenter-p1); // zenith angle of i-column
	//zenithAngleE[i] = sqrt( (recYcenter-p1)*(recYcenter-p1)*p0E*p0E + p0*p0*p1E*p1E );
	zenithAngle[i] = p0*(recYcenter[i]-p1); // zenith angle of i-column
	zenithAngleE[i] = sqrt( (recYcenter[i]-p1)*(recYcenter[i]-p1)*p0E*p0E + p0*p0*p1E*p1E );

	for(int ie=0;ie<Neta;ie++)
	{
		//double zenithAngleVSieta_fit = p0*(recY[ie]-p1);
		//double zenithAngleVSieta_fit = p0*(recY[ie]-p1) - zenithAngle[i]; 
		double zenithAngleVSieta_fit = p0*(recY[i][ie]-p1) - zenithAngle[i]; 
		// zenith angle of ieta of i-column without zenithAngle[i]
		for(int is=0;is<NSch;is++)
		{
			double zAngle = tan(zenithAngleVSieta_fit);
			//rzVSieta[i][is][ie] = atan((z[is]-z[2])*zAngle/(recY[ie]-etaYcenter))*180/3.141593;
			//rzVSieta[i][is][ie] = atan((z[is]-z[2])*zAngle/(recY[ie]-recYcenter))*180/3.141593;
			rzVSieta[i][is][ie] = atan((z[is]-z[2])*zAngle/(recY[i][ie]-recYcenter[i]))*180/3.141593;
			if(ie==0 && calcRz2 && !emptySch[i][is]) rz2[i][is] = rzVSieta[i][is][ie];
		}
	}
}


cout<<endl;
double dx1[nColumn][NSch]={0,}, dx1E[nColumn][NSch]={0,};
double dx2[nColumn][NSch]={0,}, dx2E[nColumn][NSch]={0,};
double dxT[nColumn][NSch]={0,}, dxTE[nColumn][NSch]={0,};
double sdx1[nColumn][NSch]={0,}, sdx1E[nColumn][NSch]={0,};
double sdxT[nColumn][NSch]={0,}, sdxTE[nColumn][NSch]={0,};
double dxP[nColumn][NSch]={0,}, dxPE[nColumn][NSch]={0,};
double dxM[nColumn][NSch]={0,}, dxME[nColumn][NSch]={0,};

double rzT[nColumn][NSch]={0,}, rzTE[nColumn][NSch]={0,};
double srz1[nColumn][NSch]={0,}, srz1E[nColumn][NSch]={0,};
double srzT[nColumn][NSch]={0,}, srzTE[nColumn][NSch]={0,};
double rzP[nColumn][NSch]={0,}, rzPE[nColumn][NSch]={0,};
double rzM[nColumn][NSch]={0,}, rzME[nColumn][NSch]={0,};

cout<<"Calculated d"<<v<<"1 per chamber : "<<endl; 
for(int ic=0;ic<maxNchamber;ic++) if( (viewOnlyOneColumn&&ic/10==iColumn) || !viewOnlyOneColumn ) {cout<<residualXMean[ic]<<" "; if(ic%10==9) cout<<endl;}
	
cout<<"Calculated d"<<v<<"E1 per chamber : "<<endl; 
for(int ic=0;ic<maxNchamber;ic++) if( (viewOnlyOneColumn&&ic/10==iColumn) || !viewOnlyOneColumn ) {cout<<residualXMeanE[ic]<<" "; if(ic%10==9) cout<<endl;}
cout<<endl;
	
cout<<"zenith angle per column (degree) : "; for(int i=0;i<nColumn;i++) if( (viewOnlyOneColumn&&i==iColumn) || !viewOnlyOneColumn ) cout<<zenithAngle[i]*180/3.141592<<", "; cout<<endl<<endl;

cout<<"Translated d"<<v<<" per Sch (true value) : "<<endl;
for(int it=0;it<nColumn;it++) if( (viewOnlyOneColumn&&it==iColumn) || !viewOnlyOneColumn ) {for(int is=0;is<NSch;is++) cout<<tDx[it][is]<<" "; cout<<endl;} cout<<endl;
	
cout<<"Shifted d"<<v<<" per Sch (calculated dx in previous iteration) : "<<endl;
for(int it=0;it<nColumn;it++) if( (viewOnlyOneColumn&&it==iColumn) || !viewOnlyOneColumn ) {for(int is=0;is<NSch;is++) cout<<sDx[it][is]<<" "; cout<<endl;} cout<<endl;

cout<<"Calculated d"<<v<<"1 (=d"<<v<<"'1) per Sch by ch residuals : "<<endl; 
for(int it=0;it<nColumn;it++) if( (viewOnlyOneColumn&&it==iColumn) || !viewOnlyOneColumn ) if(!emptyColumn[it])
{
	for(int is=0;is<NSch;is++)
	{
		//dx1[it][is] = (residualXMean[it][is] + residualXMean[it][is-1])/2.;
		//dx1E[it][is] = sqrt(residualXMeanE[it][is]*residualXMeanE[it][is] + residualXMeanE[it][is-1]*residualXMeanE[it][is-1])/2.;
		int idx = it*10+is*2;
		dx1[it][is] = (residualXMean[idx] + residualXMean[idx+1])/2.;
		dx1E[it][is] = sqrt(residualXMeanE[idx]*residualXMeanE[idx] + residualXMeanE[idx+1]*residualXMeanE[idx+1])/2.;
		//if(calcDx2) dx2[it][is] = zenithAngle[it] * (z[2]-z[is]);
		//dx2E[it][is] = zenithAngleE[it] * (z[2]-z[is]);
		if(calcDx2 && !emptySch[it][is])
		{
			dx2[it][is] = sin(zenithAngle[it]) * (z[2]-z[is]);
			dx2E[it][is] = sin(zenithAngleE[it]) * (z[2]-z[is]);
		}
		//dx2[it][is] = sin(zenithAngle[it]) * (z[2]-z[is]);
		//dx2E[it][is] = sin(zenithAngleE[it]) * (z[2]-z[is]);
		cout<<dx1[it][is]<<", ";

		dx1[it][is] = dx11[it][is];

		dxT[it][is] = dx1[it][is] + dx2[it][is];
		dxTE[it][is] = sqrt(dx1[it][is]*dx1[it][is] + dx2[it][is]*dx2[it][is]);
		sdx1[it][is] = tDx[it][is] - dx1[it][is];
		sdx1E[it][is] = dx1E[it][is];
		sdxT[it][is] = tDx[it][is] - dxT[it][is];
		sdxTE[it][is] = dxTE[it][is];
		dxP[it][is] = sDx[it][is] + dx1[it][is] + dx2[it][is];
		dxM[it][is] = tDx[it][is] - dxP[it][is];

		rzT[it][is] =  rz1[it][is] + rz2[it][is];
		srz1[it][is] = tRz[it][is] - rz1[it][is];
		srzT[it][is] = tRz[it][is] - rzT[it][is];

	}
	cout<<endl;
	//is = 0;
	//cout<<"Calculated d"<<v<<"E1 per Sch : "<<endl; 
	//for(int is=0;is<NSch;is++) if(is%2==1 && is/10==0)
	//{
	//	cout<<dx1E[is]<<", ";
	//	if(is%10==9) cout<<endl;
	//	is++;
	//}
}
cout<<endl;

cout<<"calculated d"<<v<<"1 per Sch using rz plots : "<<endl;
for(int it=0;it<nColumn;it++) if( (viewOnlyOneColumn&&it==iColumn) || !viewOnlyOneColumn ) {for(int is=0;is<NSch;is++) cout<<dx11[it][is]<<" "; cout<<endl;} cout<<endl;
	
cout<<"Calculated d"<<v<<"2 (=d"<<v<<"'2) per Sch by zenith angle : "<<endl; 
for(int it=0;it<nColumn;it++) if( (viewOnlyOneColumn&&it==iColumn) || !viewOnlyOneColumn ) {for(int is=0;is<NSch;is++) cout<<dx2[it][is]<<", "; cout<<endl;} cout<<endl;

cout<<"d"<<v<<"'1 + d"<<v<<"'2 per Sch : "<<endl; 
for(int it=0;it<nColumn;it++) if( (viewOnlyOneColumn&&it==iColumn) || !viewOnlyOneColumn ) {for(int is=0;is<NSch;is++) cout<<dxT[it][is]<<", "; cout<<endl;} cout<<endl;

cout<<"modulation (= d"<<v<<" - d"<<v<<"'1) per Sch : "<<endl; 
for(int it=0;it<nColumn;it++) if( (viewOnlyOneColumn&&it==iColumn) || !viewOnlyOneColumn ) {for(int is=0;is<NSch;is++) cout<<sdx1[it][is]<<", "; cout<<endl;} cout<<endl;

cout<<"modulation (= d"<<v<<" - d"<<v<<"'1 - d"<<v<<"'2) per Sch : "<<endl; 
for(int it=0;it<nColumn;it++) if( (viewOnlyOneColumn&&it==iColumn) || !viewOnlyOneColumn ) {for(int is=0;is<NSch;is++) cout<<sdxT[it][is]<<", "; cout<<endl;} cout<<endl;

//cout<<"modulation (= d"<<v<<" - sd"<<v<<" - d"<<v<<"'1 - d"<<v<<"'2) per Sch : "<<endl; 
//for(int it=0;it<nColumn;it++) if( (viewOnlyOneColumn&&it==iColumn) || !viewOnlyOneColumn ) {for(int is=0;is<NSch;is++) cout<<sdxT[it][is]-sDx[it][is]<<", "; cout<<endl;} cout<<endl;

cout<<"dx' with all iterations per Sch : "<<endl; 
for(int it=0;it<nColumn;it++) if( (viewOnlyOneColumn&&it==iColumn) || !viewOnlyOneColumn ) {for(int is=0;is<NSch;is++) cout<<dxP[it][is]<<", "; cout<<endl;} cout<<endl;

cout<<"modulation = dx-dx' with all iterations per Sch : "<<endl; 
for(int it=0;it<nColumn;it++) if( (viewOnlyOneColumn&&it==iColumn) || !viewOnlyOneColumn ) {for(int is=0;is<NSch;is++) cout<<dxM[it][is]<<", "; cout<<endl;} cout<<endl;


TLine *line1 = new TLine(0,15,0,125);
line1->SetLineWidth(1);

TCanvas *cx0 = new TCanvas("cx0","cx0",900,500);
TCanvas *cx3 = new TCanvas("cx3","cx3",900,500);
TCanvas *cx1 = new TCanvas("cx1","cx1",900,500);
TCanvas *cx2 = new TCanvas("cx2","cx2",900,500);
TCanvas *cx4 = new TCanvas("cx4","cx4",900,500);
TCanvas *cx5 = new TCanvas("cx5","cx5",900,500);
TCanvas *cx6 = new TCanvas("cx6","cx6",900,500);
cx0->Divide(3,1,0,0);
cx3->Divide(3,1,0,0);
cx1->Divide(3,1,0,0);
cx2->Divide(3,1,0,0);
cx4->Divide(3,1,0,0);
cx5->Divide(3,1,0,0);
cx6->Divide(3,1,0,0);

TGraph *g0[nColumn];
TGraph *g[nColumn];
TGraph *g3[nColumn];
TGraph *g1[nColumn];
TGraph *g2[nColumn];
TGraph *g4[nColumn];
TGraph *g5[nColumn];
TGraph *g6[nColumn];

TLegend *leg1[nColumn];
TLegend *leg2[nColumn];
TLegend *leg3[nColumn];
TLegend *leg5[nColumn];

double posX1=0.15, posY1=0.80;
double posX2=0.38, posY2=0.93;

for(int it=0;it<nColumn;it++) if( (viewOnlyOneColumn&&it==iColumn) || !viewOnlyOneColumn ) if(!emptyColumn[it])
{
	cx0->cd(it+1);
	if(it==0) htitle = Form("dx and dx' : column%d",it+1);
	else htitle = Form("column%d",it+1);
	g0[it] = new TGraph(NSch, tDx[it], z);
	g0[it]->SetTitle(htitle);
	g0[it]->GetXaxis()->SetTitle("x [cm]");
	g0[it]->GetYaxis()->SetTitle("z [cm]     ");
	g0[it]->GetXaxis()->SetLimits(-0.7,0.7);
	g0[it]->Draw("ap");
	g[it] = new TGraph(NSch, dx1[it], z);
	g[it]->SetMarkerColor(2);
	g[it]->Draw("p");
	if(it==0)
	{
		leg1[it] = new TLegend(posX1,posY1,posX2,posY2);
		leg1[it]->AddEntry(g0[it],"dx","p");
		leg1[it]->AddEntry(g[it],"dx'_{1}","p");
		leg1[it]->Draw();
	}
	line1->Draw();


	cx3->cd(it+1);
	if(it==0) htitle = Form("Modulation = dx-dx'_{1} : column%d",it+1);
	else htitle = Form("column%d",it+1);
	g3[it] = new TGraph(NSch, sdx1[it], z);
	g3[it]->SetTitle(htitle);
	g3[it]->GetXaxis()->SetTitle("x [cm]");
	g3[it]->GetYaxis()->SetTitle("z [cm]     ");
	g3[it]->GetXaxis()->SetLimits(-0.7,0.7);
	g3[it]->Draw("ap");
	line1->Draw();
	

	cx1->cd(it+1);
	if(it==0) htitle = Form("dx and dx' : column%d",it+1);
	else htitle = Form("column%d",it+1);
	g1[it] = new TGraph(NSch, dx2[it], z);
	g1[it]->SetTitle(htitle);
	g1[it]->GetXaxis()->SetTitle("x [cm]");
	g1[it]->GetYaxis()->SetTitle("z [cm]     ");
	g1[it]->GetXaxis()->SetLimits(-0.7,0.7);
	g1[it]->SetMarkerColor(3);
	g1[it]->Draw("ap");
	if(it==0)
	{
		leg2[it] = new TLegend(posX1,posY1,posX2,posY2);
		leg2[it]->AddEntry(g1[it],"dx'_{2}","p");
		leg2[it]->Draw();
	}
	line1->Draw();
	
	
	cx2->cd(it+1);
	g0[it]->Draw("ap");
	g2[it] = new TGraph(NSch, dxT[it], z);
	g2[it]->SetMarkerColor(4);
	g2[it]->Draw("p");
	if(it==0)
	{
		leg3[it] = new TLegend(posX1,posY1,posX2,posY2);
		leg3[it]->AddEntry(g0[it],"dx","p");
		leg3[it]->AddEntry(g2[it],"dx'_{1}+dx'_{2}","p");
		leg3[it]->Draw();
	}
	line1->Draw();
	
	
	cx4->cd(it+1);
	if(it==0) htitle = Form("Modulation = dx-dx'_{1}-dx'_{2} : column%d",it+1);
	else htitle = Form("column%d",it+1);
	g4[it] = new TGraph(NSch, sdxT[it], z);
	g4[it]->SetTitle(htitle);
	g4[it]->GetXaxis()->SetTitle("x [cm]");
	g4[it]->GetYaxis()->SetTitle("z [cm]     ");
	g4[it]->GetXaxis()->SetLimits(-0.7,0.7);
	g4[it]->Draw("ap");
	line1->Draw();


	cx5->cd(it+1);
	if(it==0) htitle = Form("dx and dx' : column%d",it+1);
	else htitle = Form("column%d",it+1);
	g5[it] = new TGraph(NSch, dxP[it], z);
	g5[it]->SetTitle(htitle);
	g5[it]->GetXaxis()->SetTitle("x [cm]");
	g5[it]->GetYaxis()->SetTitle("z [cm]     ");
	g5[it]->GetXaxis()->SetLimits(-0.7,0.7);
	g5[it]->SetMarkerColor(2);
	g5[it]->Draw("ap");
	g0[it]->Draw("p");
	if(it==0)
	{
		leg5[it] = new TLegend(posX1,posY1,posX2,posY2);
		leg5[it]->AddEntry(g0[it],"dx","p");
		leg5[it]->AddEntry(g5[it],"dx'","p");
		leg5[it]->Draw();
	}
	line1->Draw();


	cx6->cd(it+1);
	if(it==0) htitle = Form("Modulation = dx-dx' : column%d",it+1);
	else htitle = Form("column%d",it+1);
	g6[it] = new TGraph(NSch, dxM[it], z);
	g6[it]->SetTitle(htitle);
	g6[it]->GetXaxis()->SetTitle("x [cm]");
	g6[it]->GetYaxis()->SetTitle("z [cm]     ");
	g6[it]->GetXaxis()->SetLimits(-0.7,0.7);
	g6[it]->Draw("ap");
	line1->Draw();
}
png = savePngDir+"/run"+run+"_dx_and_dx'1_by_residuals.png";
if(printOn) cx0->Print(png);
png = savePngDir+"/run"+run+"_modulation_by_dx_and_dx'1.png";
if(printOn) cx3->Print(png);
png = savePngDir+"/run"+run+"_dx'2_by_zenith_angle.png";
if(printOn) cx1->Print(png);
png = savePngDir+"/run"+run+"_dx_and_dx'_total.png";
if(printOn) cx2->Print(png);
png = savePngDir+"/run"+run+"_modulation_by_dx_and_dx'1_dx'2.png";
if(printOn) cx4->Print(png);
png = savePngDir+"/run"+run+"_dx_and_dx'.png";
if(printOn) cx5->Print(png);
png = savePngDir+"/run"+run+"_modulation_by_dx_and_dx'.png";
if(printOn) cx6->Print(png);
cout<<endl;


cout<<"calculated dx1 per chamber using rz plots : "<<endl;
for(int ic=0;ic<maxNchamber;ic++) if( (viewOnlyOneColumn&&ic/10==iColumn) || !viewOnlyOneColumn ) {cout<<dx1_ch[ic]<<" "; if(ic%10==9) cout<<endl;} cout<<endl;

cout<<"calculated rz per chamber : "<<endl;
for(int ic=0;ic<maxNchamber;ic++) if( (viewOnlyOneColumn&&ic/10==iColumn) || !viewOnlyOneColumn ) {cout<<rz1_ch[ic]<<" "; if(ic%10==9) cout<<endl;} cout<<endl;

cout<<"calculated d"<<v<<"1 per Sch using rz plots : "<<endl;
for(int it=0;it<nColumn;it++) if( (viewOnlyOneColumn&&it==iColumn) || !viewOnlyOneColumn ) {for(int is=0;is<NSch;is++) cout<<dx11[it][is]<<" "; cout<<endl;} cout<<endl;
	
cout<<"True rz per Sch (true value) : "<<endl;
for(int it=0;it<nColumn;it++) if( (viewOnlyOneColumn&&it==iColumn) || !viewOnlyOneColumn ) {for(int is=0;is<NSch;is++) cout<<tRz[it][is]<<" "; cout<<endl;} cout<<endl;

cout<<"Rotated rz per Sch (calculated rz in previous iteration) : "<<endl;
for(int it=0;it<nColumn;it++) if( (viewOnlyOneColumn&&it==iColumn) || !viewOnlyOneColumn ) {for(int is=0;is<NSch;is++) cout<<sRz[it][is]<<" "; cout<<endl;} cout<<endl;

cout<<"calculated rz (rz'1) per Sch : "<<endl;
for(int it=0;it<nColumn;it++) if( (viewOnlyOneColumn&&it==iColumn) || !viewOnlyOneColumn ) {for(int is=0;is<NSch;is++) cout<<rz1[it][is]<<" "; cout<<endl;} cout<<endl;

cout<<"modulation (=rz-rz'1) per Sch : "<<endl;
for(int it=0;it<nColumn;it++) if( (viewOnlyOneColumn&&it==iColumn) || !viewOnlyOneColumn ) {for(int is=0;is<NSch;is++) cout<<srz1[it][is]<<" "; cout<<endl;} cout<<endl;

cout<<"zenith angle per column per ieta (degree) : "<<endl; 
for(int i=0;i<nColumn;i++) 
	if( (viewOnlyOneColumn&&i==iColumn) || !viewOnlyOneColumn ) 
	{
		for(int j=0;j<Neta;j++) cout<<zenithAngleVSieta[i][j]*180/3.141592<<", "; cout<<endl;
	}
cout<<endl;

cout<<"rz'2 per Sch per ieta (twisted angle) : "<<endl; 
for(int i=0;i<nColumn;i++) 
	if( (viewOnlyOneColumn&&i==iColumn) || !viewOnlyOneColumn ) 
	{
		for(int j=0;j<NSch;j++) 
		{
			cout<<"Sch"<<i*5+j<<" : ";
			for(int k=0;k<Neta;k++) cout<<rzVSieta[i][j][k]<<", "; cout<<endl;
			if(calcRz2 && !emptySch[i][j]) rz2[i][j] = rzVSieta[i][j][0];
			else rz2[i][j] = 0;
			rzP[i][j] = sRz[i][j] + rz1[i][j] + rz2[i][j];
			rzM[i][j] = tRz[i][j] - rzP[i][j];
		}
		cout<<endl;
	}
//cout<<endl;

//cout<<"rz'2 per Sch of ieta0 (twisted angle) : "<<endl; 
cout<<"rz'2 per Sch  : "<<endl; 
for(int it=0;it<nColumn;it++) if( (viewOnlyOneColumn&&it==iColumn) || !viewOnlyOneColumn ) {for(int is=0;is<NSch;is++) cout<<rz2[it][is]<<" "; cout<<endl;} cout<<endl;

cout<<"calculated rz (rz'1+rz'2) per Sch : "<<endl;
for(int it=0;it<nColumn;it++) if( (viewOnlyOneColumn&&it==iColumn) || !viewOnlyOneColumn ) {for(int is=0;is<NSch;is++) cout<<rzT[it][is]<<" "; cout<<endl;} cout<<endl;

cout<<"modulation (=rz-rz'1-rz'2) per Sch : "<<endl;
for(int it=0;it<nColumn;it++) if( (viewOnlyOneColumn&&it==iColumn) || !viewOnlyOneColumn ) {for(int is=0;is<NSch;is++) cout<<srzT[it][is]<<" "; cout<<endl;} cout<<endl;

//cout<<"modulation (=rz-rz'1(1)-rz'2(1)-rz'1(2)-rz'2(2)) per Sch : "<<endl;
//for(int it=0;it<nColumn;it++) if( (viewOnlyOneColumn&&it==iColumn) || !viewOnlyOneColumn ) {for(int is=0;is<NSch;is++) cout<<srzT[it][is]-sRz[it][is]<<" "; cout<<endl;} cout<<endl;

cout<<"rz' with all iterations per Sch : "<<endl; 
for(int it=0;it<nColumn;it++) if( (viewOnlyOneColumn&&it==iColumn) || !viewOnlyOneColumn ) {for(int is=0;is<NSch;is++) cout<<rzP[it][is]<<", "; cout<<endl;} cout<<endl;

cout<<"modulation = rz-rz' with all iterations per Sch : "<<endl; 
for(int it=0;it<nColumn;it++) if( (viewOnlyOneColumn&&it==iColumn) || !viewOnlyOneColumn ) {for(int is=0;is<NSch;is++) cout<<rzM[it][is]<<", "; cout<<endl;} cout<<endl;


cout<<"Next input parameters : "<<endl; 
cout<<"shiftX = [";
for(int it=0;it<nColumn;it++)
{
	if(it>0) cout<<"          ";
	for(int is=0;is<NSch;is++) {cout<<sDx[it][is]+dxT[it][is]; if(it<nColumn-1||is<NSch-1) cout<<", ";}
	if(it<nColumn-1)  cout<<"\\";
	if(it==nColumn-1) cout<<"]";
	cout<<endl;
}

cout<<endl<<"rotationZ = [";
for(int it=0;it<nColumn;it++)
{
	if(it>0) cout<<"             ";
	for(int is=0;is<NSch;is++) {cout<<sRz[it][is]+rzT[it][is]; if(it<nColumn-1||is<NSch-1) cout<<", ";}
	if(it<nColumn-1)  cout<<"\\";
	if(it==nColumn-1) cout<<"]";
	cout<<endl;
}
cout<<endl;


TCanvas *cr0 = new TCanvas("cr0","cr0",900,500);
TCanvas *cr3 = new TCanvas("cr3","cr3",900,500);
TCanvas *cr1 = new TCanvas("cr1","cr1",900,500);
TCanvas *cr2 = new TCanvas("cr2","cr2",900,500);
TCanvas *cr4 = new TCanvas("cr4","cr4",900,500);
TCanvas *cr5 = new TCanvas("cr5","cr5",900,500);
TCanvas *cr6 = new TCanvas("cr6","cr6",900,500);
cr0->Divide(3,1,0,0);
cr3->Divide(3,1,0,0);
cr1->Divide(3,1,0,0);
cr2->Divide(3,1,0,0);
cr4->Divide(3,1,0,0);
cr5->Divide(3,1,0,0);
cr6->Divide(3,1,0,0);

TGraph *gr0[nColumn];
TGraph *gr[nColumn];
TGraph *gr3[nColumn];
TGraph *gr1[nColumn];
TGraph *gr2[nColumn];
TGraph *gr4[nColumn];
TGraph *gr5[nColumn];
TGraph *gr6[nColumn];

TLegend *legr1[nColumn];
TLegend *legr2[nColumn];
TLegend *legr3[nColumn];
TLegend *legr5[nColumn];

//double posX1=0.15, posY1=0.80;
//double posX2=0.38, posY2=0.93;

for(int it=0;it<nColumn;it++) if( (viewOnlyOneColumn&&it==iColumn) || !viewOnlyOneColumn ) if(!emptyColumn[it])
{
	cr0->cd(it+1);
	if(it==0) htitle = Form("rz and rz' : column%d",it+1);
	else htitle = Form("column%d",it+1);
	gr0[it] = new TGraph(NSch, tRz[it], z);
	gr0[it]->SetTitle(htitle);
	gr0[it]->GetXaxis()->SetTitle("[degree]");
	gr0[it]->GetYaxis()->SetTitle("z [cm]     ");
	gr0[it]->GetXaxis()->SetLimits(-1.5,1.5);
	gr0[it]->Draw("ap");
	gr[it] = new TGraph(NSch, rz1[it], z);
	gr[it]->SetMarkerColor(2);
	gr[it]->Draw("p");
	if(it==0)
	{
		legr1[it] = new TLegend(posX1,posY1,posX2,posY2);
		legr1[it]->AddEntry(gr0[it],"rz","p");
		legr1[it]->AddEntry(gr[it],"rz'_{1}","p");
		legr1[it]->Draw();
	}
	line1->Draw();


	cr3->cd(it+1);
	if(it==0) htitle = Form("Modulation = rz-rz'_{1} : column%d",it+1);
	else htitle = Form("column%d",it+1);
	gr3[it] = new TGraph(NSch, srz1[it], z);
	gr3[it]->SetTitle(htitle);
	gr3[it]->GetXaxis()->SetTitle("[degree]");
	gr3[it]->GetYaxis()->SetTitle("z [cm]     ");
	gr3[it]->GetXaxis()->SetLimits(-1.5,1.5);
	gr3[it]->Draw("ap");
	line1->Draw();
	

	cr1->cd(it+1);
	if(it==0) htitle = Form("rz and rz' : column%d",it+1);
	else htitle = Form("column%d",it+1);
	gr1[it] = new TGraph(NSch, rz2[it], z);
	gr1[it]->SetTitle(htitle);
	gr1[it]->GetXaxis()->SetTitle("[degree]");
	gr1[it]->GetYaxis()->SetTitle("z [cm]     ");
	gr1[it]->GetXaxis()->SetLimits(-1.5,1.5);
	gr1[it]->SetMarkerColor(3);
	gr1[it]->Draw("ap");
	if(it==0)
	{
		legr2[it] = new TLegend(posX1,posY1,posX2,posY2);
		legr2[it]->AddEntry(gr1[it],"rz'_{2}","p");
		legr2[it]->Draw();
	}
	line1->Draw();
	
	
	cr2->cd(it+1);
	gr0[it]->Draw("ap");
	gr2[it] = new TGraph(NSch, rzT[it], z);
	gr2[it]->SetMarkerColor(4);
	gr2[it]->Draw("p");
	if(it==0)
	{
		legr3[it] = new TLegend(posX1,posY1,posX2,posY2);
		legr3[it]->AddEntry(gr0[it],"rz","p");
		legr3[it]->AddEntry(gr2[it],"rz'_{1}+rz'_{2}","p");
		legr3[it]->Draw();
	}
	line1->Draw();
	
	
	cr4->cd(it+1);
	if(it==0) htitle = Form("Modulation = rz-rz'_{1}-rz'_{2} : column%d",it+1);
	else htitle = Form("column%d",it+1);
	gr4[it] = new TGraph(NSch, srzT[it], z);
	gr4[it]->SetTitle(htitle);
	gr4[it]->GetXaxis()->SetTitle("[degree]");
	gr4[it]->GetYaxis()->SetTitle("z [cm]     ");
	gr4[it]->GetXaxis()->SetLimits(-1.5,1.5);
	gr4[it]->Draw("ap");
	line1->Draw();


	cr5->cd(it+1);
	if(it==0) htitle = Form("rz and rz' : column%d",it+1);
	else htitle = Form("column%d",it+1);
	gr5[it] = new TGraph(NSch, rzP[it], z);
	gr5[it]->SetTitle(htitle);
	gr5[it]->GetXaxis()->SetTitle("[degree]");
	gr5[it]->GetYaxis()->SetTitle("z [cm]     ");
	gr5[it]->GetXaxis()->SetLimits(-1.5,1.5);
	gr5[it]->SetMarkerColor(2);
	gr5[it]->Draw("ap");
	gr0[it]->Draw("p");
	if(it==0)
	{
		legr5[it] = new TLegend(posX1,posY1,posX2,posY2);
		legr5[it]->AddEntry(gr0[it],"rz","p");
		legr5[it]->AddEntry(gr5[it],"rz'","p");
		legr5[it]->Draw();
	}
	line1->Draw();


	cr6->cd(it+1);
	if(it==0) htitle = Form("Modulation = rz-rz' : column%d",it+1);
	else htitle = Form("column%d",it+1);
	gr6[it] = new TGraph(NSch, rzM[it], z);
	gr6[it]->SetTitle(htitle);
	gr6[it]->GetXaxis()->SetTitle("[degree]");
	gr6[it]->GetYaxis()->SetTitle("z [cm]     ");
	gr6[it]->GetXaxis()->SetLimits(-1.5,1.5);
	gr6[it]->Draw("ap");
	line1->Draw();
}
png = savePngDir+"/run"+run+"_rz_and_rz'1_by_residuals.png";
if(printOn) cr0->Print(png);
png = savePngDir+"/run"+run+"_modulation_by_rz_and_rz'1.png";
if(printOn) cr3->Print(png);
png = savePngDir+"/run"+run+"_rz'2_by_zenith_angle.png";
if(printOn) cr1->Print(png);
png = savePngDir+"/run"+run+"_rz_and_rz'_total.png";
if(printOn) cr2->Print(png);
png = savePngDir+"/run"+run+"_modulation_by_rz_and_rz'1_rz'2.png";
if(printOn) cr4->Print(png);
png = savePngDir+"/run"+run+"_rz_and_rz'.png";
if(printOn) cr5->Print(png);
png = savePngDir+"/run"+run+"_modulation_by_rz_and_rz'.png";
if(printOn) cr6->Print(png);
cout<<endl;

double dxMMean[nColumn] = {0,};
double dxMSigma[nColumn] = {0,};
double dxMTot[nColumn] = {0,}; //total mis-alignment for dx
double rzMMean[nColumn] = {0,};
double rzMSigma[nColumn] = {0,};
double rzMTot[nColumn] = {0,}; //total mis-alignment for rz
double dxMMean0[nColumn] = {0,};
double dxMTot0[nColumn] = {0,}; //total mis-alignment for dx at iteration0
double rzMMean0[nColumn] = {0,};
double rzMTot0[nColumn] = {0,}; //total mis-alignment for rz at iteration0
for(int it=0;it<nColumn;it++)
{
	int NSch2 = 0;
	for(int is=0;is<NSch;is++) if(!emptySch[it][is])
	{
		dxMMean[it] += dxM[it][is];
		rzMMean[it] += rzM[it][is];
		dxMMean0[it] += tDx[it][is];
		rzMMean0[it] += tRz[it][is];
		NSch2++;
	}
	dxMMean[it] /= NSch2;
	rzMMean[it] /= NSch2;
	dxMMean0[it] /= NSch2;
	rzMMean0[it] /= NSch2;
	for(int is=0;is<NSch;is++) if(!emptySch[it][is])
	{
		dxMSigma[it] += (dxMMean[it]-dxM[it][is])*(dxMMean[it]-dxM[it][is]);
		dxMTot[it] += fabs(dxMMean[it]-dxM[it][is]);
		dxMTot0[it] += fabs(dxMMean0[it]-tDx[it][is]);
		rzMSigma[it] += (rzMMean[it]-rzM[it][is])*(rzMMean[it]-rzM[it][is]);
		rzMTot[it] += fabs(rzMMean[it]-rzM[it][is]);
		rzMTot0[it] += fabs(rzMMean0[it]-tRz[it][is]);
	}
	dxMSigma[it] = sqrt(dxMSigma[it]/(NSch2-1));
	rzMSigma[it] = sqrt(rzMSigma[it]/(NSch2-1));
}

for(int it=0;it<nColumn;it++) cout<<"column"<<it+1<<", dxMMean "<<dxMMean[it]<<", dxMSigma "<<dxMSigma[it]<<", dxMTot "<<dxMTot[it]<<", dxMTot0 "<<dxMTot0[it]<<endl;
for(int it=0;it<nColumn;it++) cout<<"column"<<it+1<<", rzMMean "<<rzMMean[it]<<", rzMSigma "<<rzMSigma[it]<<", rzMTot "<<rzMTot[it]<<", rzMTot0 "<<rzMTot0[it]<<endl;
cout<<endl;

double dxMax = 0;
for(int it=0;it<nColumn;it++) for(int is=0;is<NSch;is++) if(dxMax < fabs(dxT[it][is])) dxMax = fabs(dxT[it][is]);
double rzMax = 0;
for(int it=0;it<nColumn;it++) for(int is=0;is<NSch;is++) if(rzMax < fabs(rzT[it][is])) rzMax = fabs(rzT[it][is]);
cout<<"dxMax "<<dxMax<<", rzMax "<<rzMax<<endl<<endl;


FILE *outFile = NULL;
outFile = fopen("table/table_output_run"+run+".txt","w");
//outFile = fopen("table_output_run"+run+".txt","w");

float dxOut[nColumn][NSch];
float rzOut[nColumn][NSch];
for(int it=0;it<nColumn;it++) 
	for(int is=0;is<NSch;is++)
	{
		dxOut[it][is] = sDx[it][is]+dxT[it][is];
		rzOut[it][is] = sRz[it][is]+rzT[it][is];
	}
cout<<"Next dx and rz table : "<<endl;
for(int is=NSch-1;is>=0;is--)
{
	for(int it=0;it<nColumn;it++) 
	{
		float outNum = dxOut[it][is];
		fprintf(outFile, "%.3f ", outNum);
		cout<<outNum<<" ";
	}
	fprintf(outFile,"\n");
	cout<<endl;
}
fprintf(outFile,"\n");
cout<<endl;
for(int is=NSch-1;is>=0;is--)
{
	for(int it=0;it<nColumn;it++) 
	{
		float outNum = rzOut[it][is];
		fprintf(outFile, "%.3f ", outNum);
		cout<<outNum<<" ";
	}
	fprintf(outFile,"\n");
	cout<<endl;
}
fclose(outFile);

}


