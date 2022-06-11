#include "TMath.h"
#include "TCanvas.h"
#include "TSpline.h"
#include "TAxis.h"
// Declare Global Variables, 
// which can be referred from anywhere after their declaration: 
const Long64_t n = 4;
Double_t Point_x[n] = {0, 0.5, 1, 2};
Double_t Point_y[n] = {1, 0.368, 0.135, 0.018};


// main function with the same name of this file
// it will be executed when you type ".x interpolation.C" in ROOT
void interpolation()
{
// Section 1. Draw the points on a canvas
   TCanvas *c1 = new TCanvas("c1","interpolation",0,0,1000,800);

   TGraph * defGraph  = new TGraph(n,Point_x,Point_y);
   TSpline3 * defSplineDefault = new TSpline3("splinedefault", defGraph);
   //https://www.mpp.mpg.de/~jingliu/ECPI/interpolation.C
   // "b2e2" together with the last two "0" means that the second derivatives 
   // of the begin and end points equal to zero
   TSpline3 * defSplineb2e2    = new TSpline3("splineb2e2", Point_x, Point_y, n, "b2e2", 0, 0);

   const int ntest = 100000;

   TGraph * calcGraphLinear   = new TGraph(ntest);
   TGraph * calcGraphSpline   = new TGraph(ntest);
   TGraph * calcSplineDefault = new TGraph(ntest);
   TGraph * calcSplineb2e2    = new TGraph(ntest);

   for(int ii=0; ii<=ntest; ii++){
     const double xx = Point_x[0]+ii*(Point_x[3]-Point_x[0])/ntest;
     calcGraphLinear->SetPoint(ii,xx, defGraph->Eval(xx));
     calcGraphSpline->SetPoint(ii, xx, defGraph->Eval(xx,0x0,"S"));
     calcSplineDefault->SetPoint(ii, xx, defSplineDefault->Eval(xx));
     calcSplineb2e2->SetPoint(ii, xx, defSplineb2e2->Eval(xx));
     printf("ii %d xx %f Graph_Linear %f Graph_Spline %f Spline3_Default %f Spline3_b2e2 %f\n", ii, xx, 
            calcGraphLinear->GetPointY(ii),
            calcGraphSpline->GetPointY(ii),
            calcSplineDefault->GetPointY(ii),
            calcSplineb2e2->GetPointY(ii));
   }

   defGraph->SetMarkerStyle(24);
   defGraph->SetMarkerSize(2);

   calcGraphLinear->SetMarkerStyle(20);
   calcGraphSpline->SetMarkerStyle(21);
   calcSplineDefault->SetMarkerStyle(22);
   calcSplineb2e2->SetMarkerStyle(23);
   
   calcGraphLinear->SetMarkerSize(2);
   calcGraphSpline->SetMarkerSize(2);
   calcSplineDefault->SetMarkerSize(2);
   calcSplineb2e2->SetMarkerSize(2);

   calcGraphLinear->SetMarkerColor(kRed);
   calcGraphSpline->SetMarkerColor(kBlue);
   calcSplineDefault->SetMarkerColor(kMagenta);
   calcSplineb2e2->SetMarkerColor(kBlack);
   
   calcGraphLinear->SetLineColor(kRed);
   calcGraphSpline->SetLineColor(kBlue);
   calcSplineDefault->SetLineColor(kMagenta);
   calcSplineb2e2->SetLineColor(kBlack);

   
   calcGraphLinear->SetLineWidth(1);
   calcGraphSpline->SetLineWidth(2);
   calcSplineDefault->SetLineWidth(1);
   calcSplineb2e2->SetLineWidth(1);
   
   c1->cd();

   
   defGraph->GetYaxis()->SetRangeUser(-0.2, 1.2);
   defGraph->GetXaxis()->SetLimits(-0.5, 2.5);
   defGraph->Draw("ap");
      
   calcGraphLinear->Draw("same l");
   calcGraphSpline->Draw("same l");
   calcSplineDefault->Draw("same l");
   calcSplineb2e2->Draw("same l");

}

int main()
{
  return 0;
}
