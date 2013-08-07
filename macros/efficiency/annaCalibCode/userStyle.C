/*
 *  userStyle.C
 */
#include "TStyle.h"

void userStyle() {
	TStyle *user = new TStyle("user","Basic TDR Style");
	
// For the canvas:
  user->SetCanvasBorderMode(0);
  user->SetCanvasColor(kWhite);
  user->SetCanvasDefH(600); //Height of canvas
  user->SetCanvasDefW(600); //Width of canvas
  user->SetCanvasDefX(0);   //POsition on screen
  user->SetCanvasDefY(0);

// For the Pad:
  user->SetPadBorderMode(0);
  // user->SetPadBorderSize(Width_t size = 1);
  user->SetPadColor(kWhite);
  user->SetPadGridX(false);
  user->SetPadGridY(false);
  user->SetGridColor(0);
  user->SetGridStyle(3);
  user->SetGridWidth(1);

// For the frame:
  user->SetFrameBorderMode(0);
  user->SetFrameBorderSize(1);
  user->SetFrameFillColor(0);
  user->SetFrameFillStyle(0);
  user->SetFrameLineColor(1);
  user->SetFrameLineStyle(1);
  user->SetFrameLineWidth(2);

// For the histo:
  // user->SetHistFillColor(1);
  // user->SetHistFillStyle(0);
  user->SetHistLineColor(1);
  user->SetHistLineStyle(0);
  user->SetHistLineWidth(2);
  // user->SetLegoInnerR(Float_t rad = 0.5);
  // user->SetNumberContours(Int_t number = 20);

  user->SetEndErrorSize(2);
  //  user->SetErrorMarker(20);
  user->SetErrorX(0.);
  user->SetMarkerStyle(20);

//For the fit/function:
  //user->SetOptFit(1);
  user->SetOptFit(1111);
  user->SetFitFormat("5.4g");
  user->SetFuncColor(2);
  user->SetFuncStyle(1);
  user->SetFuncWidth(2);

//For the date:
  user->SetOptDate(0);
  // user->SetDateX(Float_t x = 0.01);
  // user->SetDateY(Float_t y = 0.01);

// For the statistics box:
  user->SetOptFile(0);
  user->SetOptStat(0); // To display the mean and RMS:   SetOptStat("mr");
  user->SetStatColor(kWhite);
  user->SetStatFont(42);
  user->SetStatFontSize(0.025);
  user->SetStatTextColor(1);
  user->SetStatFormat("6.4g");
  user->SetStatBorderSize(1);
  user->SetStatH(0.1);
  user->SetStatW(0.15);
  // user->SetStatStyle(Style_t style = 1001);
  // user->SetStatX(Float_t x = 0);
  // user->SetStatY(Float_t y = 0);

// Margins:
  user->SetPadTopMargin(0.1);
  user->SetPadBottomMargin(0.15);
  user->SetPadLeftMargin(0.15);
  user->SetPadRightMargin(0.15);

// For the Global title:
  user->SetOptTitle(0);
  user->SetTitleFont(42);
  user->SetTitleColor(1);
  user->SetTitleTextColor(1);
  user->SetTitleFillColor(10);
  user->SetTitleFontSize(0.05);
  // user->SetTitleH(0); // Set the height of the title box
  // user->SetTitleW(0); // Set the width of the title box
  // user->SetTitleX(0); // Set the position of the title box
  // user->SetTitleY(0.985); // Set the position of the title box
  // user->SetTitleStyle(Style_t style = 1001);
  // user->SetTitleBorderSize(2);

// For the axis titles:
  user->SetTitleColor(1, "XYZ");
  user->SetTitleFont(42, "XYZ");
  user->SetTitleSize(0.045, "XYZ");
  // user->SetTitleXSize(Float_t size = 0.02); // Another way to set the size?
  // user->SetTitleYSize(Float_t size = 0.02);
  user->SetTitleXOffset(1.2);
  user->SetTitleYOffset(1.6);
  // user->SetTitleOffset(1.1, "Y"); // Another way to set the Offset

// For the axis labels:
  user->SetLabelColor(1, "XYZ");
  user->SetLabelFont(42, "XYZ");
  user->SetLabelOffset(0.005, "XYZ");
  user->SetLabelSize(0.04, "XYZ");

// For the axis:
  user->SetAxisColor(1, "XYZ");
  user->SetStripDecimals(kTRUE);
  user->SetTickLength(0.03, "XYZ");
  user->SetNdivisions(510, "XYZ");
  user->SetPadTickX(1);  // To get tick marks on the opposite side of the frame
  user->SetPadTickY(1);


// Change for log plots:
  user->SetOptLogx(0);
  user->SetOptLogy(0);
  user->SetOptLogz(0);

// Postscript options:
  user->SetPaperSize(20.,20.);
  // user->SetLineScalePS(Float_t scale = 3);
  // user->SetLineStyleString(Int_t i, const char* text);
  // user->SetHeaderPS(const char* header);
  // user->SetTitlePS(const char* pstitle);

  // user->SetBarOffset(Float_t baroff = 0.5);
  // user->SetBarWidth(Float_t barwidth = 0.5);
  // user->SetPaintTextFormat(const char* format = "g");
  // user->SetPalette(Int_t ncolors = 0, Int_t* colors = 0);
  // user->SetTimeOffset(Double_t toffset);
  // user->SetHistMinimumZero(kTRUE);

  user->cd();

}
