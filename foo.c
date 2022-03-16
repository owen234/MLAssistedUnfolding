#include "TMPalette.hh"
TMPalette *hot3Palette = 0;

void testPalettesBad(){
  // Here the Palettes are defined
  gROOT->Macro("Palettes.C");
  
  gStyle->SetOptStat(0);

  TCanvas *c = new TCanvas("c","c",1000,625);
  c->Divide(2,2);
  
  TH2F *h1 = new TH2F("h1","h1",100,-4,4,100,-4,4);
  TH2F *h2 = new TH2F("h2","h2",100,-4,4,100,-4,4);
  TH2F *h3 = new TH2F("h3","h3",100,-4,4,100,-4,4);
  TH2F *h4 = new TH2F("h4","h4",100,-4,4,100,-4,4);
  
  Double_t a,b;
  for (Int_t i=0;i<50000;i++) {
    gRandom->Rannor(a,b);
    h1->Fill(a-1.5,b-1.5);
    h2->Fill(a+1.5,b+1.5);
    h3->Fill(a-1.5,b+1.5);
    h4->Fill(a+1.5,b-1.5);
  }

  hot3Palette = (TMPalette*) gROOT->FindObject("hot3");
  if(!hot3Palette) {
    const Int_t hot3NRGBs = 4;
    const Int_t hot3NCont = 64;
    Double_t hot3Stops[hot3NRGBs] = { 0.00, 0.25, 0.75, 1.00 }; 
    Double_t hot3Red[hot3NRGBs] =   { 0.91, 1.000, 1.000, 0.160 }; 
    Double_t hot3Green[hot3NRGBs] = { 0.91, 0.7, 0.149, 0.160 };
    Double_t hot3Blue[hot3NRGBs] =  { 0.91, 0.000, 0.000, 0.160 };
    hot3Palette = new TMPalette("hot3");
    hot3Palette->CreateGradientColorTable(hot3NRGBs, hot3Stops,hot3Red, hot3Green, hot3Blue, hot3NCont);                                       
  }
  

  TExec *ex1 = new TExec("ex1","electronPalette->cd();");
  TExec *ex2 = new TExec("ex2","hot3Palette->cd();");
  TExec *ex3 = new TExec("ex3","rbowPalette->cd();");
  TExec *ex4 = new TExec("ex4","yellowPalette->cd();");
  
  c->cd(1);
  h1->Draw("axis");
  ex1->Draw();
  h1->Draw("colz same");
  
  c->cd(2);
  h2->Draw("axis");
  ex2->Draw();
  h2->Draw("colz same");
  
  c->cd(3);
  h3->Draw("axis");
  ex3->Draw();
  h3->Draw("colz same");
  
  c->cd(4);
  h4->Draw("axis");
  ex4->Draw();
  h4->Draw("colz same");

  c->Print("./testPalettes.png");
}
