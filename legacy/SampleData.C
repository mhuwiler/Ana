#include <vector>
#include <string>
#include <unordered_map>
#include "TCut.h"
#include "TColor.h"



class SampleData 
{
public: 

	SampleData(int initColor, std::string initLegend, TCut initCut = "", std::string initLatex = "", std::initializer_list<std::string> initRefList = {})
		: color(initColor), legend(initLegend), cut(initCut), latex(initLatex), fileRefs(initRefList.begin(), initRefList.end()) 
	{
		if (initLatex == "") 
		{
			latex = TString::Format("$%s$", TString(initLegend).ReplaceAll("#", "\\").ReplaceAll("rightarrow", "rightarrow ").Data()).Data(); // If no latex string provided, use the legend and make it proper latex
		}
	}; 

	SampleData() = default;

	int color = -999.; 
	std::string legend = ""; 
	TCut cut = ""; 
	std::string latex = ""; 
	std::vector<std::string> fileRefs; 

}; 


template<typename S, typename T>
class RedirectingMap 
{
	public:
		inline RedirectingMap() = default; 
		inline RedirectingMap(std::unordered_map<S, T> list, std::unordered_map<S, S> mapping = {}) : map(list), redirections(mapping) {}; 

		inline T at(S key) 
		{ 
			if (map.find(key) == map.end()) 
			{
				if (redirections.find(key) != redirections.end()) 
				{
					return map.at(redirections.at(key)); 
				}
				else 
				{
					printWarning(key); 
					return SampleData();  
				}
			}
			else 
			{
				return map.at(key); 
			}
		};

		inline T operator[](S key) 
		{
			if (map.find(key) == map.end()) key = redirections.at(key);

			return map[key]; 
		};


		inline void SetRedirections(std::unordered_map<S, S> mapping) { redirections = mapping; };

		inline void printWarning(const std::string key) { std::cout << "Warning: key " << key << " not in map nor redirects. " << std::endl; };

	private:
		std::unordered_map<S, T> map; 

		std::unordered_map<S, S> redirections; 
};



namespace Ana {
	
	RedirectingMap<std::string, SampleData> InitSamples() 
	{
		// https://colorbrewer2.org/?type=diverging&scheme=RdYlBu&n=7#type=diverging&scheme=RdYlBu&n=11
		std::vector<std::vector<int> > colors = 
		{	

			{164, 176, 37}, //{165,0,38},
			{215,48,39}, 
			{244,109,67}, 
			{253,174,97}, 
			{254,224,144}, 
			{255,255,191}, 
			{224,243,248}, 
			{171,217,233}, 
			{116,173,209}, 
			{69,117,180}, 
			{49,54,149}
		}; 

		std::vector<Int_t> mycolors; 
		std::vector<TColor*> rootcolors; 
		mycolors.reserve(colors.size()); 
		rootcolors.reserve(colors.size()); 
	
		Float_t colorintmax = 256.; 
		for (auto color : colors) 
		{
			//std::cout << color->GetNumber() << std::endl; 
			assert(color.size() == 3); 
			TColor *newcolor = new TColor(TColor::GetFreeColorIndex(), static_cast<Float_t>(color[0])/colorintmax, static_cast<Float_t>(color[1])/colorintmax, static_cast<Float_t>(color[2])/colorintmax);
			rootcolors.push_back(newcolor); 
			mycolors.push_back(newcolor->GetNumber()); 
		}

		TCut genMatchCut = "(Dstar_match)&&(b_tau_match)";

		std::unordered_map<std::string, SampleData> samplelist = {
			// B0 decays
			{"Sig", {
				mycolors[0], 
				"signal", 
				genMatchCut, 
				"$B^0\\rightarrow D^{*-}\\tau^+\\nu_\\tau$", 
				{"Sig"}}
			}, 
			{"SigTest", {
				mycolors[0], 
				"signal", 
				genMatchCut, 
				"$B^0\\rightarrow D^{*-}\\tau^+\\nu_\\tau$", 
				{"Sig"}}
			}, 
			{"B0toDstarDs", SampleData(
				mycolors[1], 
				"B^{0}#rightarrowD*D_{s}", 
				genMatchCut, 
				"$B^0\\rightarrow D^{*-}D_{s}^+$", 
				{"B0toDstarDs"})
			},
			{"B0toDstarDsstar", SampleData(
				mycolors[2], 
				"B^{0}#rightarrowD*D*_{s}", 
				genMatchCut, 
				"$B^0\\rightarrow D^{*-}D_{s}^{*+}$", 
				{"B0toDstarDsstar"})
			},
			{"B0toDstarDs1", SampleData(
				mycolors[2], 
				"B^{0}#rightarrowD*D_{s1}", 
				genMatchCut, 
				"$B^{0}\\rightarrow D^{*}D_{s1}$", 
				{"B0toDstarDs1"})
			},
			{"B0toDstarDs0star", SampleData(
				mycolors[2], 
				"B^{0}#rightarrow D*D_{s0}*", 
				genMatchCut, 
				"$B^{0}\\rightarrow D^{*}D_{s0}^{*}$", 
				{"B0toDstarDs0star"})
			},
			{"B0toDstarD", SampleData(
				mycolors[3], 
				"B^{0}#rightarrowD*D", 
				genMatchCut, 
				"", 
				{"B0toDstarD"})
			},
			{"B0toDstarD0K", SampleData(
				mycolors[5], 
				"B^{0}#rightarrowD*D^{0}K", 
				genMatchCut, 
				"$B^0\\rightarrow D^{*}D^{0}K$", 
				{"B0toDstarD0K"})
			},
			{"B0toDstarD0Kstar", SampleData(
				mycolors[6], 
				"B^{0}#rightarrowD*D^{0}K*", 
				genMatchCut, 
				"$B^0\\rightarrow D^{*}D^{0}K^{*}$", 
				{"B0toDstarD0K"})
			},
			{"B0toDstara1", SampleData(
				mycolors[7], 
				"B^{0}#rightarrowD*a_{1}", 
				genMatchCut, 
				"$B^0\\rightarrow D^{*-}a_{1}^+$", 
				{"B0toDstara1"})
			},
			{"B0toDstar3pi", SampleData(
				mycolors[8], 
				"B^{0}#rightarrowD*3#pi", 
				genMatchCut, 
				"$B^{0}\\rightarrow D^{*}3\\pi$", 
				{"B0toDstar3pi"})
			},
			{"B0toDstar3pipi0", SampleData(
				mycolors[9], 
				"B^{0}#rightarrow D*3#pi#pi^{0}", 
				genMatchCut, 
				"$B^{0}\\rightarrow D^{*}3\\pi\\pi^{0}$", 
				{"B0toDstar3pipi0"})
			},
			{"B0toDstar5pi", SampleData(
				mycolors[7], 
				"B^{0}#rightarrowD*5#pi", 
				genMatchCut, 
				"$B^{0}\\rightarrow D^{*}5\\pi$", 
				{"B0toDstar5pi"})
			},
			{"B0toDstarrho0pi", SampleData(
				mycolors[10], 
				"B^{0}#rightarrowD*#rho^{0}#pi", 
				genMatchCut, 
				"$B^{0}\\rightarrow D^{*}\\rho^{0}\\pi$", 
				{"B0toDstarrho0pi"})
			},
			{"B0toDstarKpipi", SampleData(
				mycolors[11], 
				"B^{0}#rightarrowD*K#pi#pi", 
				genMatchCut, 
				"$B^{0}\\rightarrow D^{*}K\\pi\\pi$", 
				{"B0toDstarKpipi"})
			},
			{"B0toDstarKKstar", SampleData(
				mycolors[11], 
				"B^{0}#rightarrowD*KK*", 
				genMatchCut, 
				"$B^{0}\\rightarrow D^{*}KK^{*}$", 
				{"B0toDstarKKstar"})
			},
			// Bu decays 
			{"ButoDstarDK", SampleData(
				mycolors[4], 
				"B^{+}#rightarrowD*DK", 
				genMatchCut, 
				"$B^{+}\\rightarrow D^{*}DK$", 
				{"ButoDstarDK"})
			}, 
			{"ButoDstarpipipipi0", SampleData(
				mycolors[4], 
				"B^{+}#rightarrowD*3#pi#pi^{0}", 
				genMatchCut, 
				"$B^{+}\\rightarrow D^{*0}3\\pi\\pi^{0}$", 
				{"ButoDstarpipipipi0"})
			}, 
			{"ButoDstarpipipi", SampleData(
				mycolors[4], 
				"B^{+}#rightarrowD*3#pi", 
				genMatchCut, 
				"$B^{+}\\rightarrow D^{*0}3\\pi$", 
				{"ButoDstarpipipi"})
			}, 
			{"ButoDstarpipipi0", SampleData(
				mycolors[4], 
				"B^{+}#rightarrowD*#pi#pi#pi^{0}", 
				genMatchCut, 
				"$B^{+}\\rightarrow D^{*}\\pi\\pi\\pi^{0}$", 
				{"ButoDstarpipipi0"})
			}, 
			{"B0toDstarDsX", SampleData(
				mycolors[4], 
				"B^{+}#rightarrowD*D_{s}X", 
				genMatchCut, 
				"", 
				{"B0toDstarDsX"})
			}, 
			{"ButoDstarDK", SampleData(
				mycolors[4], 
				"B^{+}#rightarrowD*DK", 
				genMatchCut, 
				"", 
				{"ButoDstarDK"})
			}, 
			{"ButoDstarXc", SampleData(
				mycolors[4], 
				"B^{+}#rightarrowD*X_{c}", 
				genMatchCut, 
				"", 
				{"ButoDstarXc"})
			}, 
			// To be completed
			{"B0toDD", SampleData(
				mycolors[8], 
				"B^{0}#rightarrowDD", 
				genMatchCut, 
				"$B^0\\rightarrow D_{(s)}^{(*)}D_{(s)}^{(*)}$", 
				{"B0toB0DD"})
			},
			{"BstoDD", SampleData(
				mycolors[8], 
				"B^{s}#rightarrowDD", 
				genMatchCut, 
				"$B_{s}^{0}\\rightarrow D_{(s)}^{(*)}D_{(s)}^{(*)}$", 
				{"BstoDD"})
			},
			{"data", SampleData(
				 mycolors[10], 
				"data", 
				"1", 
				"data", 
				{"dataD1", "dataD2", "dataD3", "dataD4", "dataD5", "dataA1", "dataA2", "dataA3", "dataA4", "dataA5", "dataB1", "dataB2", "dataB3", "dataB4", "dataB5", "dataC1", "dataC2", "dataC3", "dataC4", "dataC5"})
			},
			{"dataD", SampleData(
				 mycolors[10], 
				"data", 
				"1", 
				"data", 
				{"dataD1", "dataD2", "dataD3", "dataD4", "dataD5"})
			},
			{"dataA", SampleData(
				 mycolors[10], 
				"data", 
				"1", 
				"data", 
				{"dataA1", "dataA2", "dataA3", "dataA4", "dataA5"})
			},
			{"dataB", SampleData(
				 mycolors[10], 
				"data", 
				"1", 
				"data", 
				{"dataB1", "dataB2", "dataB3", "dataB4", "dataB5"})
			},
			{"dataC", SampleData(
				 mycolors[10], 
				"data", 
				"1", 
				"data", 
				{"dataC1", "dataC2", "dataC3", "dataC4", "dataC5"})
			},
			{"dataDWS", SampleData(
				 mycolors[10], 
				"data", 
				"1", 
				"data", 
				{"dataD1WS", "dataD2WS", "dataD3WS", "dataD4WS", "dataD5WS"})
			},
			{"WSD", SampleData(
				 mycolors[10], 
				"data", 
				"1", 
				"data", 
				{"dataD1WS", "dataD2WS", "dataD3WS", "dataD4WS", "dataD5WS"})
			},
			{"dataAWS", SampleData(
				 mycolors[10], 
				"data", 
				"1", 
				"data", 
				{"dataA1WS", "dataA2WS", "dataA3WS", "dataA4WS", "dataA5WS"})
			},
			{"dataBWS", SampleData(
				 mycolors[10], 
				"data", 
				"1", 
				"data", 
				{"dataB1WS", "dataB2WS", "dataB3WS", "dataB4WS", "dataB5WS"})
			},
			{"dataCWS", SampleData(
				 mycolors[10], 
				"data", 
				"1", 
				"data", 
				{"dataC1WS", "dataC2WS", "dataC3WS", "dataC4WS", "dataC5WS"})
			},
			// data driven components
			{"WS", SampleData(
				 mycolors[10], 
				"|q_{B}| = 2  WS", 
				"1", 
				"$WS~|q_{B}|}~=~2$", 
				{"dataD1WS", "dataD2WS", "dataD3WS", "dataD4WS", "dataD5WS", "dataA1WS", "dataA2WS", "dataA3WS", "dataA4WS", "dataA5WS", "dataB1WS", "dataB2WS", "dataB3WS", "dataB4WS", "dataB5WS", "dataC1WS", "dataC2WS", "dataC3WS", "dataC4WS", "dataC5WS"})
			},
			{"WSTau", SampleData(
				 mycolors[10], 
				"|q_{#tau}| = 3  WS", 
				"1", 
				"$WS~|q_{#tau}|~=~3$", 
				{"WSTau"})
			},
		};

		RedirectingMap<std::string, SampleData> samples = samplelist; 

		samples.SetRedirections({
			{"dataA1", "data"}, 
			{"dataA2", "data"}, 
			{"dataA3", "data"}, 
			{"dataA4", "data"}, 
			{"dataA5", "data"}, 
			{"dataB1", "data"}, 
			{"dataB2", "data"}, 
			{"dataB3", "data"}, 
			{"dataB4", "data"}, 
			{"dataB5", "data"}, 
			{"dataC1", "data"}, 
			{"dataC2", "data"},
			{"dataC3", "data"}, 
			{"dataC4", "data"}, 
			{"dataC5", "data"}, 
			{"dataD1", "data"}, 
			{"dataD2", "data"}, 
			{"dataD3", "data"}, 
			{"dataD4", "data"}, 
			{"dataD5", "data"}, 
			{"dataA", "data"}, 
			{"dataB", "data"}, 
			{"dataC", "data"}, 
			{"dataD", "data"}, 
			{"data1", "data"}, 
			{"data2", "data"}, 
			{"data3", "data"}, 
			{"data4", "data"}, 
			{"data5", "data"}, 
			{"dataA1WS", "WS"}, 
			{"dataA2WS", "WS"}, 
			{"dataA3WS", "WS"}, 
			{"dataA4WS", "WS"}, 
			{"dataA5WS", "WS"}, 
			{"dataB1WS", "WS"}, 
			{"dataB2WS", "WS"}, 
			{"dataB3WS", "WS"}, 
			{"dataB4WS", "WS"}, 
			{"dataB5WS", "WS"}, 
			{"dataC1WS", "WS"}, 
			{"dataC2WS", "WS"},
			{"dataC3WS", "WS"}, 
			{"dataC4WS", "WS"}, 
			{"dataC5WS", "WS"}, 
			{"dataD1WS", "WS"}, 
			{"dataD2WS", "WS"}, 
			{"dataD3WS", "WS"}, 
			{"dataD4WS", "WS"}, 
			{"dataD5WS", "WS"}, 
			{"dataAWS", "WS"}, 
			{"dataBWS", "WS"}, 
			{"dataCWS", "WS"}, 
			{"dataDWS", "WS"}, 
			{"data1WS", "WS"}, 
			{"data2WS", "WS"}, 
			{"data3WS", "WS"}, 
			{"data4WS", "WS"}, 
			{"data5WS", "WS"}, 
			{"dataAB", "data"}, 
			{"dataABC", "data"}, 
			{"dataBC", "data"},
			{"ABCD", "WS"}
		}); 

		return samples; 
	}
} // namespace Ana

