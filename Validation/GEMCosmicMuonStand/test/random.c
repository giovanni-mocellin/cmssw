{

gRandom->SetSeed(time(0));
const int nX = 15;
double randDx[nX] = {0,};
double randRz[nX] = {0,};

cout<<"trueDx = [";
for(int i=0;i<nX;i++)
{
	if(i%5==0&&i>0) cout<<"          ";
	randDx[i] = double(int(gRandom->Uniform(-6,6)))/10.;
	cout<<randDx[i]<<", ";
	if(i%5==4&&i<nX-1) cout<<"\\";
	if(i==nX-1) cout<<"]";
	if(i%5==4&&i>0) cout<<endl;
}
cout<<endl;

cout<<"trueRz = [";
for(int i=0;i<nX;i++)
{
	if(i%5==0&&i>0) cout<<"          ";
	randRz[i] = double(int(gRandom->Uniform(-11,11)))/10.;
	cout<<randRz[i]<<", ";
	if(i%5==4&&i<nX-1) cout<<"\\";
	if(i==nX-1) cout<<"]";
	if(i%5==4&&i>0) cout<<endl;
}
cout<<endl;

}

