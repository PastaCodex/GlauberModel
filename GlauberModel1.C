//For use with ROOT library from CERN

//A function to generate a nucleus composed of n nucleons
void genNucleus(const int nNucleons, TRandom3* random, TF1* radius,double* x, double* y, double* z){
  
  for (int point = 0; point<nNucleons; point++){
    double phi = random-> Uniform(0,TMath::TwoPi());
    double costheta = random-> Uniform(-1, 1);
    double theta = acos(costheta);
    double r = radius -> GetRandom();
    x[point] = r*sin(theta)*cos(phi);
    y[point] = r*sin(theta)*sin(phi);
    z[point] = r*cos(theta);
    
  }
  return;
}

//Calculates the centrality of a collision along the x axis
void ImpactParameterShift(const int nNucleons, double* x){
  
  
  TF1* impactp = new TF1("impactp","x*TMath::TwoPi()", 0, 17);
  double b = 0.;
  b = impactp-> GetRandom();

  for (int i = 0; i<nNucleons; i++){
    x[i] += b;
  }
  
  return;
}

//Calculates the positions of all nucleons during impact in the xy-plane with respect to every other nucleon to calculate the number of nucleon-nucleon interactions 
void genCollParti(const int nNucleons, size_t* nCollisions, double* Ax, double* Ay, double* Bx, double* By, bool* ParA, bool* ParB){
 
  size_t c = 1;
  for (size_t l = 0; l<nNucleons; l++){
    for(size_t m = 0; m<nNucleons; m++){
      if (pow(Bx[l] - Ax[m], 2) + pow(By[l] - Ay[m], 2) < 6.5/TMath::Pi()){
        ParA[m] = 1;
        ParB[l] = 1;
	c += 1; 
      }
    }
  }
  nCollisions[0] = c;
  return;
}

//-----------------------------------------------------------------------------------------//

void GlauberModel1(){                               
  
  TRandom3* random = new TRandom3(0); //Generates seed
  
  const int numberOfNucleons(208); //Number of nucleons in a desired nucleus. In this case Pb
  
  double AxCoord[numberOfNucleons] = {0};
  double AyCoord[numberOfNucleons] = {0};
  double AzCoord[numberOfNucleons] = {0};
  double BxCoord[numberOfNucleons] = {0};
  double ByCoord[numberOfNucleons] = {0};
  double BzCoord[numberOfNucleons] = {0};

  //Declaring physical properties of the desired nucleus
  double rhonought = 3.0;
  double nuclearradius = 6.62;
  double skindepth = 0.546;

  TCanvas* RadialProb = new TCanvas("RadialProb", "Radial Probability", 500, 500);

  //Determine the radial position of generated nucleons
  TF1* radialProbability = new TF1("radialProbability", "x*x*[0]/(1 + exp((x-[1])/[2])); Radius (fm); Radial Probability", 0, 11);
  radialProbability -> SetParameter(0, rhonought);
  radialProbability -> SetParameter(1, nuclearradius);
  radialProbability -> SetParameter(2, skindepth);
  radialProbability -> Draw();
  TF1* radialDensity = new TF1("radialDensity", "[0]/(1 + exp((x-[1])/[2])); Radius (fm); Radial Density p(r)", 0, 11);
  TCanvas* Radial = new TCanvas("Radial", "Radial Density", 500, 500);
  radialDensity -> SetParameter(0, rhonought);
  radialDensity -> SetParameter(1, nuclearradius);
  radialDensity -> SetParameter(2, skindepth);
  radialDensity -> Draw();
  
  genNucleus(numberOfNucleons, random, radialProbability, AxCoord, AyCoord, AzCoord);
  genNucleus(numberOfNucleons, random, radialProbability, BxCoord, ByCoord, BzCoord);
  
  ImpactParameterShift(numberOfNucleons, BxCoord);
  
  size_t numberOfCollisions[1];
  size_t numberOfParticipants[1];
  bool ParticipantA[numberOfNucleons] = {0};
  bool ParticipantB[numberOfNucleons] = {0};
  
  genCollParti(numberOfNucleons, numberOfCollisions, AxCoord, AyCoord, BxCoord, ByCoord, ParticipantA, ParticipantB);

  for (n = 0; n<numberOfNucleons; n++){
    cout << "P[" << n << "] = " << ParticipantA[n] << endl; 
  }
  cout << "# of Coll = " << numberOfCollisions[1] << endl;
  //cout << "# of Participants = " << numberOfParticipants << endl;

//----------------------------------------------------------------------------------------//

  TCanvas* NucleusXY = new TCanvas("NucleusXY", "X-Y Collision Graph", 500, 500);
  TGraph* NucleusAxy = new TGraph(numberOfNucleons, AxCoord, AyCoord);
  TGraph* NucleusBxy = new TGraph(numberOfNucleons, BxCoord, ByCoord);
  
  NucleusAxy-> Draw("AP*");
  NucleusAxy-> GetHistogram()-> GetXaxis()-> SetLimits(-15, 15);
  NucleusAxy-> GetHistogram()-> SetMaximum(15);
  NucleusAxy-> GetHistogram()-> SetMinimum(-15);
  NucleusAxy-> GetHistogram()-> SetTitle("X-Y Nucleon Distribution");
  NucleusAxy-> GetXaxis()-> SetTitle("X(fm)");
  NucleusAxy-> GetYaxis()-> SetTitle("Y(fm)");
  NucleusAxy-> SetMarkerStyle(8);
  NucleusAxy-> SetMarkerColor(8);
  NucleusAxy-> SetMarkerSize(3);
  
  NucleusBxy-> Draw("PSame");
  NucleusBxy-> GetHistogram()-> GetXaxis()-> SetLimits(-15, 15);
  NucleusBxy-> GetHistogram()-> SetMaximum(15);
  NucleusBxy-> GetHistogram()-> SetMinimum(-15);
  NucleusBxy-> SetMarkerStyle(8);
  NucleusBxy-> SetMarkerColor(46);
  NucleusBxy-> SetMarkerSize(3);
  
  TCanvas* NucleusZX = new TCanvas("NucleusZX", "X-Z Collision Graph", 500, 500);
  TGraph* NucleusAzx = new TGraph(numberOfNucleons, AzCoord, AxCoord);
  TGraph* NucleusBzx = new TGraph(numberOfNucleons, BzCoord, BxCoord);
  NucleusAzx-> Draw("AP*");
  NucleusAzx-> GetHistogram()-> GetXaxis()-> SetLimits(-15, 15);
  NucleusAzx-> GetHistogram()-> SetMaximum(15);
  NucleusAzx-> GetHistogram()-> SetMinimum(-15);
  NucleusAzx-> GetHistogram()-> SetTitle("X-Z Nucleon Distribution");
  NucleusAzx-> GetXaxis()-> SetTitle("Z(fm)");
  NucleusAzx-> GetYaxis()-> SetTitle("X(fm)");
  NucleusAzx-> SetMarkerStyle(8);
  NucleusAzx-> SetMarkerColor(8);
  NucleusAzx-> SetMarkerSize(3);
  
  NucleusBzx-> Draw("Psame");
  NucleusBzx-> GetHistogram()-> GetXaxis()-> SetLimits(-15, 15);
  NucleusBzx-> GetHistogram()-> SetMaximum(15);
  NucleusBzx-> GetHistogram()-> SetMinimum(-15);
  NucleusBzx-> SetMarkerStyle(8);
  NucleusBzx-> SetMarkerColor(46);
  NucleusBzx-> SetMarkerSize(3);
  
  return;
  
  
}
