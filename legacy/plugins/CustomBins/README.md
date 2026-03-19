# CustomBins

CustomBins is a class for handling custom bins for any ROOT histogram class (TH1 or TH2 derived). The class offers nothing special but doing the simple math necessary when using the ROOT histograms with custom bins. It provides the conversion from std::vector to C type array, and ensures the value of Nbins is correct w.r.t. the number of bin edges. It also enables generating uniform bins, following the scheme of TH1 histograms constructor (Nbins, Min, Max) for uniform binning. The aims of this class are to enable developers to always use the custom bin constructor of the ROOT histogram, even in case of uniform binning, for maximum flexibility of the code, and to avoid the pain of handling C type array sizes and number of bins manually. The bin information must be provided at construction under one of the following forms: 

	CustomBins(std::vector<Double_t> binEdges)

Will construct custom bins with the bin edges provided in the vector, the first value being the lower edge of the first bin and the last value the upper edge of the last bin (same behaviour as when constructing a ROOT histogram with a C type array of custom bins). 

	CustomBins(N, Min, Max)

Will construct N uniform bins between Min and Max. Min is the lower edge of the first bin and Max the upper edge of the last bin. 

	CustomBins(Size, Double_t binEdges)

This offers the same behaviour as the first constructor, with binEdges being a C type array and Size its size (and not the number of bins!). This is provided for compatibility reasons, and mainly for use in python where conversion is made from numpy arrays to C type arrays. 


Then, the class is used with its two methods: getNBins() that returns the number of bins (number of bin edges -1), and getBinEdges() which returns the C type array of bin edges. 

	TH1D("myfunkyhisto", "myfunkyhisto", CustomBins.GetNBins(), CustomBins.GetBinEdges())

The class might contain a pointer to the histogram in the future, or become part of a wrapper for histograms if it proves useful... 
