// elements/TriggerEmbedding.C
//
// Trigger tools and corrections for the analysis. It contains L1 trigger prescale corrections
// as well as code for the application of trigger scale factors. 
//
// Functions provided (all in namespace Ana):
//
//   prescaleWeight(const int ...) : function returning the L1 prescale weight for MC samples
//


namespace Ana {

// Function computing the L1 prescale correction to be applied to all MC samples. It returns an event weight, to factor in to the total event weights on all MC samples. 
float prescaleWeight(const int L1_HTT200er, const int L1_HTT255er, const int L1_HTT280er, const int L1_HTT320er, const int L1_HTT360er, const int L1_HTT400er, const int L1_HTT450er, const int L1_ETT2000, const int L1_SingleJet180, const int L1_SingleJet200, const int L1_DoubleJet30er2p5_Mass_Min250_dEta_Max1p5, const int L1_DoubleJet30er2p5_Mass_Min300_dEta_Max1p5, const int L1_DoubleJet30er2p5_Mass_Min330_dEta_Max1p5) 
{
		float weight = 1.; 
		// If only triggered by L1_HTT200er, which is prescaled by 1600, apply corresponding weight 
		if (L1_HTT200er && !(L1_HTT255er || L1_HTT280er || L1_HTT320er || L1_HTT360er || L1_HTT400er || L1_HTT450er || L1_ETT2000 || L1_SingleJet180 || L1_SingleJet200 || L1_DoubleJet30er2p5_Mass_Min250_dEta_Max1p5 || L1_DoubleJet30er2p5_Mass_Min300_dEta_Max1p5 || L1_DoubleJet30er2p5_Mass_Min330_dEta_Max1p5)) weight = 1./1600.; 

		// If triggered by L1_HTT255er, which is prescaled by 500, and not any seed with lower prescale, apply corresponding weight 
		if (L1_HTT255er && !(L1_HTT280er || L1_HTT320er || L1_HTT360er || L1_HTT400er || L1_HTT450er || L1_ETT2000 || L1_SingleJet180 || L1_SingleJet200 || L1_DoubleJet30er2p5_Mass_Min250_dEta_Max1p5 || L1_DoubleJet30er2p5_Mass_Min300_dEta_Max1p5 || L1_DoubleJet30er2p5_Mass_Min330_dEta_Max1p5)) weight = 1./500.;  


		return weight; 
}

} // namespace Ana