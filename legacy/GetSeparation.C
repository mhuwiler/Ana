/* DO WHAT THE FUCK YOU WANT TO PUBLIC LICENSE
 Version 2, December 2004
 
 Copyright (C) 2017 Conor Fitzpatrick <sam@hocevar.net>
 
 Everyone is permitted to copy and distribute verbatim or modified
 copies of this license document, and changing it is allowed as long
 as the name is changed.
 
 DO WHAT THE FUCK YOU WANT TO PUBLIC LICENSE
 TERMS AND CONDITIONS FOR COPYING, DISTRIBUTION AND MODIFICATION
 
 0. You just DO WHAT THE FUCK YOU WANT TO. */

#include "TH1D.h"
#include "TH1F.h"


Double_t GetSeparation( const TH1D& S, const TH1D& B ){
    // compute "separation" defined as
    // <s2> = (1/2) Int_-oo..+oo { (S^2(x) - B^2(x))/(S(x) + B(x)) dx }
    Double_t separation = 0;
    // sanity checks
    // signal and background histograms must have same number of bins and
    // same limits
    if ((S.GetNbinsX() != B.GetNbinsX()) || (S.GetNbinsX() <= 0)) 
    {
        cout << "<GetSeparation> signal and background"
             << " histograms have different number of bins: "
             << S.GetNbinsX() << " : " << B.GetNbinsX() << endl;
        return -1.; 
    }
    if (S.GetXaxis()->GetXmin() != B.GetXaxis()->GetXmin() ||
        S.GetXaxis()->GetXmax() != B.GetXaxis()->GetXmax() ||
        S.GetXaxis()->GetXmax() <= S.GetXaxis()->GetXmin()) 
    {
        cout << S.GetXaxis()->GetXmin() << " " << B.GetXaxis()->GetXmin()
             << " " << S.GetXaxis()->GetXmax() << " " << B.GetXaxis()->GetXmax()
             << " " << S.GetXaxis()->GetXmax() << " " << S.GetXaxis()->GetXmin() << endl;
        cout << "<GetSeparation> signal and background"
             << " histograms have different or invalid dimensions:" << endl;
        return -1.; 
    }

    Int_t    nbins = S.GetNbinsX();
    Double_t nS    = S.Integral( 0, nbins+1, "width" ); // include under/overflow bins
    Double_t nB    = B.Integral( 0, nbins+1, "width" );

    if (nS > 0 && nB > 0) {
            // include under/overflow bins
            for (Int_t ibin=0; ibin<=nbins+1; ibin++) {
                    Double_t s = S.GetBinContent( ibin )/nS;
                    Double_t b = B.GetBinContent( ibin )/nB;
                    // separation

                    if (s + b > 0) separation += 0.5*(s - b)*(s - b)/(s + b)*S.GetXaxis()->GetBinWidth(ibin);
                    //separation += 0.5*(s - b)*(s - b)/(s + b)*S.GetXaxis()->GetBinWidth(ibin);

                    //   cout << nS << "    " << nB << "    " << separation << endl;
            }
    }
    else {
            cout << "<GetSeparation> histograms with zero entries: "
                    << nS << " : " << nB << " cannot compute separation"
                    << endl;
            separation = 0;
    }
    //cout << separation << endl;
    return separation;
}
