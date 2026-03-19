mkdir -p Thesisplots

VERSION=v7.2

root -l cliPlotting.C ../../data/${VERSION}/B0toDstarDs_tauDNN.root -q -e '_file0->cd(); ntuplizer->cd(); tree->Draw("b_tau_m", "b_B_match"); xLegend("m(D_{s})"); yLegend(""); Style(); gPad->Print("Thesisplots/Dsmasspipipi.pdf");'

#DNN plots
root -l cliPlotting.C ../../data/${VERSION}/Sig_tauDNN.root -q -e '_file0->cd(); ntuplizer->cd(); 
	tree->Draw("b_tau_m", "b_B_match"); xLegend("m(D_{s})"); yLegend(""); Style(); gPad->Print("Thesisplots/Dsmasspipipi.pdf");
	tree->Draw("b_tau_m", "b_B_match"); xLegend("m(D_{s})"); yLegend(""); Style(); gPad->Print("Thesisplots/Dsmasspipipi.pdf");'

#final selection plots
root -l cliPlotting.C ../../data/${VERSION}/dataD1_tauDNN.root -q -e '_file0->cd(); ntuplizer->cd(); tree->Draw("b_B_m:b_tau_m>>htemp(50, 0.5, 1.8, 50, 2.5, 5.2)", "", "COLZ"); xLegend("m(#tau)"); yLegend("m(B)"); Style(); gPad->Print("Thesisplots/mBvstaudata.pdf");'
root -l cliPlotting.C ../../data/${VERSION}/Sig_tauDNN.root -q -e '_file0->cd(); ntuplizer->cd(); tree->Draw("b_B_m:b_tau_m>>htemp(50, 0.5, 1.8, 50, 2.5, 5.2)", "", "COLZ"); xLegend("m(#tau)"); yLegend("m(B)"); Style(); gPad->Print("Thesisplots/mBvstausig.pdf");'
