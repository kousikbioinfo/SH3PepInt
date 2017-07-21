#include "Utility.h"
#include "Histogram.h"
#include "BaseGraphClass.h"
#include "GraphClass.h"
#include "NSPDK_FeatureGenerator.h"
#include "wrapper.h"
#include "vectors.h"
#include <ctime>
#include <numeric>
#include <list>


using namespace std;

//---------------------------------------------------------------------------------
const string CREDITS("Name: NeighborhoodSubgraphPaiwiseDistanceKernel\nVersion: 7.4\nProgrammer: Fabrizio Costa\nLast update: June 21 2011");
//---------------------------------------------------------------------------------
ofstream LOGF("log",std::ios_base::app);

FlagsService& The_FlagsService= FlagsService::get_instance();

inline double random01(){
  return (double)rand()/(double)RAND_MAX;
}

inline unsigned randomUnsigned(unsigned aMax){
  return (unsigned)(rand() % aMax);
}

inline int IntHashSimple(int key, int aModulo=RAND_MAX)
{
  key = ~key + (key << 15); // key = (key << 15) - key - 1;
  key = key ^ (key >> 12);
  key = key + (key << 2);
  key = key ^ (key >> 4);
  key = key * 2057; // key = (key + (key << 3)) + (key << 11);
  key = key ^ (key >> 16);
  return key%aModulo;
}


//-----------------------------------------------------------------------------
// MurmurHash2, by Austin Appleby

// Note - This code makes a few assumptions about how your machine behaves -

// 1. We can read a 4-byte value from any address without crashing
// 2. sizeof(int) == 4

// And it has a few limitations -

// 1. It will not work incrementally.
// 2. It will not produce the same results on little-endian and big-endian
//    machines.

unsigned int MurmurHash2 ( const void * key, int len, unsigned int seed )
{
	// 'm' and 'r' are mixing constants generated offline.
	// They're not really 'magic', they just happen to work well.

	const unsigned int m = 0x5bd1e995;
	const int r = 24;

	// Initialize the hash to a 'random' value

	unsigned int h = seed ^ len;

	// Mix 4 bytes at a time into the hash

	const unsigned char * data = (const unsigned char *)key;

	while(len >= 4)
	{
		unsigned int k = *(unsigned int *)data;

		k *= m; 
		k ^= k >> r; 
		k *= m; 
		
		h *= m; 
		h ^= k;

		data += 4;
		len -= 4;
	}
	
	// Handle the last few bytes of the input array

	switch(len)
	{
	case 3: h ^= data[2] << 16;
	case 2: h ^= data[1] << 8;
	case 1: h ^= data[0];
	        h *= m;
	};

	// Do a few final mixes of the hash to ensure the last few
	// bytes are well-incorporated.

	h ^= h >> 13;
	h *= m;
	h ^= h >> 15;

	return h;
} 

inline int IntHash(int key, int aModulo=RAND_MAX, unsigned aSeed=0){
  const double A=sqrt(2)-1;
  return IntHashSimple(key*(aSeed+1)*A,aModulo);
  //return MurmurHash2(&key,aModulo,aSeed);
}


//------------------------------------------------------------------------------------------------
class Timer {
public:
  typedef double diff_type;

  // Same as Timer t; t.begin();
  Timer(): start(std::clock()), elapsed(0) {}
  ~Timer(){cout<<"Elapsed time: "<<end()<<endl;}
  // Last result before a call to begin()
  diff_type last() const { return elapsed; }
  // Reset the timer
  void begin() { start = std::clock(); elapsed = 0; }
  // Save the result
  diff_type end(){
    elapsed = (diff_type)std::clock() - start;
    elapsed /= CLOCKS_PER_SEC;
    return elapsed;
  }

private:
  std::clock_t start;
  diff_type    elapsed;
};

class ProgressBar{
 public:
  ProgressBar():mCounter(0){}
  ~ProgressBar(){cout<<endl<<"Counted "<<mCounter<<" times."<<endl;}
  void Begin(){mCounter=0;}
  void Count(){
    mCounter++;
    if (mCounter<1000){
      if (mCounter%10==0) cout<<"."<<flush;
      if (mCounter%100==0) cout<<"|"<<flush;
    }
    if (mCounter<1000000) 
      if (mCounter%1000==0) cout<<mCounter/(1000)<<"K "<<flush;
    if (mCounter%1000000==0) cout<<mCounter/(1000000)<<"M "<<flush;
  }
  unsigned End(){return mCounter;}  
 private:
  unsigned mCounter;
  Timer mTimer;
};

//---------------------------------------------------------------------------------
class VectorClass{
  friend ostream& operator<<(ostream& out,const VectorClass& aV){aV.Output(out);return out;}
public:
  VectorClass(){}
  VectorClass(unsigned aSize){Init(aSize);}
  void operator=(const VectorClass& aVector){
    mV=aVector.mV;
  }
  VectorClass(const VectorClass& aVector){
    (*this)=aVector;
  }
  VectorClass(const vector<double>& aVector){
    mV=aVector;
  }
  void Init(unsigned aSize){
    mV.clear();
    for (unsigned i=0;i<aSize;++i)
      mV.push_back(-1);
  }
  void Import(const string& aFileName){
    mV.clear();
    ifstream fin;
    fin.open(aFileName.c_str());
    if (!fin){cerr<<"Cannot open file: "<<aFileName<<endl;throw exception();}
    string line;

    //read size
    while(getline(fin,line)){
      if (line!=""){
	stringstream ss;
	ss<<line;
	while(ss.good()){
	  string value_str;
	  ss>>value_str;
	  if (value_str!="") mV.push_back(stream_cast<double>(value_str));
	}
      }
    }
    fin.close();
  }
  void Clear(){mV.clear();}
  unsigned Size()const{return mV.size();}
  ostream& Output(ostream& out)const{
    for (unsigned i=0;i<Size();i++)
	out<<mV[i]<<" ";
    return out;
  }
  void PushBack(double aValue){mV.push_back(aValue);}
  double& operator[](unsigned i){assert(i<Size()); return mV[i];}
  double operator[](unsigned i)const{assert(i<Size()); return mV[i];}
  double Sum()const{
    double avg=0;
    for (unsigned i=0;i<Size();i++) avg+=mV[i];
    return avg;
  }
  double Mean()const{
    return Sum()/Size();
  }
  double StandardDeviation()const{
    double avg=Mean();
    double sd=0;
    for (unsigned i=0;i<Size();i++) sd+=(mV[i]-avg)*(mV[i]-avg);
    sd=sd/(Size()-1);
    sd=sqrt(sd);
    return sd;
  }
  double Order(double aOrder)const{
    vector<double> v(mV);
    sort(v.begin(),v.end());
    if (aOrder==0) return v[0];
    else if (aOrder==1) return v[v.size()-1];
    else return v[(unsigned)((double)Size()*aOrder)];
  }
  double Median()const{
    return Order(.5);
  }
  double MedianAbsoluteDifference()const{
    double median=Median();
    VectorClass v;
    for (unsigned i=0;i<Size();i++) v.PushBack(fabs(mV[i]-median));
    return v.Median();
  }
  double Min()const{
    return Order(0);
  }
  double Max()const{
    return Order(1);
  }
  VectorClass RemoveNulls(){
    VectorClass v;
    for (unsigned i=0;i<Size();i++) if(mV[i]!=-1) v.PushBack(mV[i]);
    return v;
  }
  ostream& OutputStatistics(ostream& out){
    VectorClass v=RemoveNulls();
    if (v.Size()>0){
      out<<"num: "<<v.Size()<<" sum: "<<v.Sum()<<" avg: "<<v.Mean()<<" sd: "<<v.StandardDeviation();
      out<<" min: "<<v.Min()<<" Q.05: "<<v.Order(.05)<<" Q.25: "<<v.Order(.25)<<" Q.5: "<<v.Median()<<" Q.75: "<<v.Order(.75)<<" Q.95: "<<v.Order(.95)<<" max: "<<v.Max();
    } else {}
    return out;
  }


protected:
  vector<double> mV;
};

class ParameterWrapperClass{
public:
  ParameterWrapperClass():mGspanInputFileName(""),
			  mSparseBinaryInputFileName(""),
			  mSparseASCIIInputFileName(""),
			  mSparseASCIITestInputFileName(""),
			  mTrainTargetInputFileName(""),
			  mRadiusMax(1),
			  mDistanceMax(4),
			  mMatchingType("hard"),
			  mVerbose(false),
			  mFeatureBitSize(30),
			  mMinKernel(false),
			  mType("nspdk"),
			  mNormalization(true),
			  mNumNearestNeighbors(10),
			  mNumHashFunctions(30),
			  mSampleSize(10),
			  mNonRedundantFilter(1),
			  mAddAbstractMatchingType(false),
			  mOutputFeatures(false),
			  mOutputFeatureMap(false),
			  mOutputKernel(false),  
			  mOutputApproximateKNN(false),  
			  mOutputFastKNN(false),  
			  mOutputTrueKNN(false),  
			  mOutputCluster("NONE"),  
			  mOutputApproximateKNNPrediction(false),
			  mOutputTrueKNNPrediction(false),
			  mOutputDatasetSplit(false),
			  mOutputHashEncoding(false),
			  mCaching(true),
			  mTrueSort(true),
			  mEccessNeighbourSizeFactor(3),
			  mHashFactor(1),
			  mNumCenters(10),
			  mSizeThreshold(30),
			  mImbalanceTolerance(1.5),
			  mAddSoftMatchingType(false),
			  mAddGappedMatchingType(false),
			  mAddStemInsensitiveMatchingType(false),
			  mAddPPInsensitiveMatchingType(false),
			  mWeightFileName(""),
			  mOutputWeightAnnotation(false),
			  mWhiteListFileName(""),
			  mBlackListFileName(""),
			  mMaxIntersectionSize(0),
			  mDebug(0){}
  
  void Usage(string aCommandName){
    cerr<<"Usage: "<<aCommandName<<endl
	<<"-fg <file name gspan format> for input from file "<<endl
	<<"-fsb <file name sparse format binary> for input from file"<<endl
	<<"-fsa <file name sparse format ascii> for input from file"<<endl
	<<"-fsats <file name sparse format ascii> for input from file for test set"<<endl
	<<"-ftrt <file name> for target for train set"<<endl
	<<"-wl <file name for white listed instances ids (1 based)> for data set"<<endl
	<<"-bl <file name for black listed instances ids (1 based)> for data set"<<endl
	<<"[-of flag to output feature encoding (default: "<< mOutputFeatures<<")]"<<endl
	<<"[-ofm flag to output feature map encoding (default: "<< mOutputFeatureMap<<")]"<<endl
	<<"[-ok flag to output kernel matrix (default: "<< mOutputKernel<<")]"<<endl
	<<"[-oaknn flag to output approximate k-nearest neighburs (default: "<< mOutputApproximateKNN<<")]"<<endl
	<<"[-ofknn flag to output fast k-nearest neighburs (default: "<< mOutputFastKNN<<")]"<<endl
	<<"[-otknn flag to output true (i.e. implies full kernel matrix evaluation) k-nearest neighburs (default: "<< mOutputTrueKNN<<")]"<<endl	
	<<"[-oc DENSITY|GREEDY|SAMPLE to output clusters (default: "<< mOutputCluster<<")]"<<endl
	<<"[-oaknnp flag to output approximate knn prediction (default: "<< mOutputApproximateKNNPrediction<<")]"<<endl
	<<"[-otknnp flag to output approximate knn prediction (default: "<< mOutputTrueKNNPrediction<<")]"<<endl
	<<"[-ods flag to output both training set and test set split (default: "<< mOutputDatasetSplit<<")]"<<endl
	<<"[-ohe flag to output hash encoding (default: "<< mOutputHashEncoding<<")]"<<endl
	<<"[-b <feature space bits size> (default: "<< mFeatureBitSize <<")]"<<endl
	<<"[-R <max radius> (default: "<< mRadiusMax <<")]"<<endl
	<<"[-D <max distance relations> (default: "<< mDistanceMax <<")]"<<endl
	<<"[-T <nspdk|anspdk|mnspdk|nspdk3d|alnspdk|nspdkvp> (default: "<< mType <<")] Note:use mnspdk for speed up in weight annotation"<<endl
	<<"[-t <hard | soft | hard_soft> as neighborhood matching use HARD=exact matching, use SOFT=attribute matching with root identifier as radius 0 neighborhood, use HARD-SOFT=attribute matching with root identifier as full neighborhood encoding (default: "<< mMatchingType <<")]"<<endl
	<<"[-t+a flag to add abstract matching (default: "<< mAddAbstractMatchingType <<")]"<<endl
	<<"[-t+s flag to add soft matching (default: "<< mAddSoftMatchingType <<")]"<<endl
	<<"[-t+g flag to add gapped invariant matching (default: "<< mAddGappedMatchingType <<")]"<<endl
	<<"[-t+si flag to add stem insensitive matching (default: "<< mAddStemInsensitiveMatchingType <<")]"<<endl
	<<"[-t+pp flag to add purine-pyrimidine insensitive matching (default: "<< mAddPPInsensitiveMatchingType <<")]"<<endl
	<<"[-mink flag to set minimum kernel rather than dot product (default: "<< mMinKernel <<")]"<<endl
	<<"[-nn flag to de-acivate normalization (default: "<< !mNormalization <<")]"<<endl
	<<"[-nhf <num hash functions> for the Locality Sensitive Hashing function (default: "<< mNumHashFunctions <<")]"<<endl
	<<"[-hf <hash factor> for the Locality Sensitive Hashing function (default: "<< mHashFactor <<")]"<<endl
	<<"[-knn <num nearest neighbors> (default: "<< mNumNearestNeighbors<<")]"<<endl
	<<"[-ensf <eccess neighbour size factor> (default: "<< mEccessNeighbourSizeFactor<<") (0 to avoid trimming)] "<<endl
	<<"[-ss <sample size> for clustering procedure (default: "<< mSampleSize <<")]"<<endl
	<<"[-nrt <similarity filtering of redundant centers [0,1]> for clustering procedure (default: "<< mNonRedundantFilter <<") (the smaller the less similar the centers)]"<<endl
	<<"[-no-cache flag to deactivate caching of kernel value computation (to minimize memory usage) (default: "<< !mCaching <<")]"<<endl
	<<"[-no-true-sort flag to deactivate sorting approximate neighbours with true kernel computation (default: "<< !mTrueSort <<")]"<<endl
	<<"[-nc num centers (default: "<< mNumCenters <<")]"<<endl
	<<"[-st <size threshold> (default: "<< mSizeThreshold <<")]"<<endl
	<<"[-it <imbalance tolerance> (default: "<< mImbalanceTolerance <<")]"<<endl
	<<"[-fw <file name of weights associated to features>]"<<endl
	<<"[-ow flag to output gspan file for input annotated with weight (default: "<< mOutputWeightAnnotation <<")]"<<endl
	<<"[-mi <maximum number of elements in intersection between center neighborhoods> (default: "<< mMaxIntersectionSize <<")]"<<endl
	<<"[-debug <debug level> (default: "<< mDebug <<")]"<<endl;

    exit(0);
  }

  void Init(int argc, char** argv){
    vector<string> options;
    for (int i=1;i<argc;i++) options.push_back(argv[i]);
    for (vector<string>::iterator it=options.begin();it!=options.end();++it)
      {
	if ((*it)=="-h" || (*it)=="--help") Usage(argv[0]); 
	else if ((*it)=="-fg") mGspanInputFileName=(*(++it)); 
	else if ((*it)=="-fsb") mSparseBinaryInputFileName=(*(++it)); 
	else if ((*it)=="-fsa") mSparseASCIIInputFileName=(*(++it)); 
	else if ((*it)=="-fsats") mSparseASCIITestInputFileName=(*(++it)); 
	else if ((*it)=="-ftrt") mTrainTargetInputFileName=(*(++it)); 
	else if ((*it)=="-wl") mWhiteListFileName=(*(++it)); 
	else if ((*it)=="-bl") mBlackListFileName=(*(++it)); 
	else if ((*it)=="-knn") mNumNearestNeighbors=stream_cast<unsigned>(*(++it)); 
	else if ((*it)=="-ensf") mEccessNeighbourSizeFactor=stream_cast<double>(*(++it)); 
	else if ((*it)=="-hf") mHashFactor=stream_cast<unsigned>(*(++it)); 
	else if ((*it)=="-R") mRadiusMax=stream_cast<double>(*(++it)); 
	else if ((*it)=="-D") mDistanceMax=stream_cast<double>(*(++it)); 
	else if ((*it)=="-t") mMatchingType=(*(++it));
	else if ((*it)=="-t+a") mAddAbstractMatchingType=true;
	else if ((*it)=="-t+s") mAddSoftMatchingType=true;
	else if ((*it)=="-t+g") mAddGappedMatchingType=true;
	else if ((*it)=="-t+si") mAddStemInsensitiveMatchingType=true;
	else if ((*it)=="-t+pp") mAddPPInsensitiveMatchingType=true;
	else if ((*it)=="-v") mVerbose=true;
	else if ((*it)=="-mink") mMinKernel=true;
	else if ((*it)=="-b") mFeatureBitSize=stream_cast<int>(*(++it)); 
	else if ((*it)=="-T") mType=(*(++it)); 
	else if ((*it)=="-of") mOutputFeatures=true;
	else if ((*it)=="-ofm") mOutputFeatureMap=true;
	else if ((*it)=="-ok") mOutputKernel=true;
	else if ((*it)=="-oaknn") mOutputApproximateKNN=true;
	else if ((*it)=="-ofknn") mOutputFastKNN=true;
	else if ((*it)=="-otknn") mOutputTrueKNN=true;
	else if ((*it)=="-oaknnp") mOutputApproximateKNNPrediction=true;
	else if ((*it)=="-otknnp") mOutputTrueKNNPrediction=true;
	else if ((*it)=="-ods") mOutputDatasetSplit=true;
	else if ((*it)=="-ohe") mOutputHashEncoding=true;
	else if ((*it)=="-oc") mOutputCluster=(*(++it)); 
	else if ((*it)=="-nn") mNormalization=false;
	else if ((*it)=="-debug") mDebug=stream_cast<int>(*(++it)); 
	else if ((*it)=="-nhf") mNumHashFunctions=stream_cast<unsigned>(*(++it)); 
	else if ((*it)=="-ss") mSampleSize=stream_cast<unsigned>(*(++it)); 
	else if ((*it)=="-nrt") mNonRedundantFilter=stream_cast<double>(*(++it)); 
	else if ((*it)=="-no-cache") mCaching=false;
	else if ((*it)=="-no-true-sort") mTrueSort=false;
	else if ((*it)=="-nc") mNumCenters=stream_cast<unsigned>(*(++it)); 
	else if ((*it)=="-st") mSizeThreshold=stream_cast<unsigned>(*(++it)); 
	else if ((*it)=="-it") mImbalanceTolerance=stream_cast<double>(*(++it)); 
	else if ((*it)=="-fw") mWeightFileName=(*(++it)); 
	else if ((*it)=="-ow") mOutputWeightAnnotation=true;
	else if ((*it)=="-mi") mMaxIntersectionSize=stream_cast<unsigned>(*(++it)); 
	else {cerr<<"Unrecognized parameter: "<<(*it)<<"."<<endl;throw exception();}
      }
    if (!(mMatchingType=="hard"||mMatchingType=="soft"||mMatchingType=="hard_soft"||mMatchingType=="multiview"||mMatchingType=="mixed")) {cerr<<"Wrong value for parameter: -t: "<<mMatchingType<<endl;throw exception();}
  }

public:
  string mGspanInputFileName;
  string mSparseBinaryInputFileName;
  string mSparseASCIIInputFileName;
  string mSparseASCIITestInputFileName;
  string mTrainTargetInputFileName;
  double mRadiusMax;
  double mDistanceMax;
  string mMatchingType;
  bool mVerbose;
  int mFeatureBitSize;
  bool mMinKernel;
  string mType;
  bool mNormalization;
  unsigned mNumNearestNeighbors;
  unsigned mNumHashFunctions;
  unsigned mSampleSize;
  double mNonRedundantFilter;
  bool mAddAbstractMatchingType;
  bool mOutputFeatures;
  bool mOutputFeatureMap;
  bool mOutputKernel; 
  bool mOutputApproximateKNN; 
  bool mOutputFastKNN; 
  bool mOutputTrueKNN; 
  string mOutputCluster; 
  bool mOutputApproximateKNNPrediction;
  bool mOutputTrueKNNPrediction;
  bool mOutputDatasetSplit;
  bool mOutputHashEncoding;
  bool mCaching;
  bool mTrueSort;
  double mEccessNeighbourSizeFactor;
  unsigned mHashFactor;
  unsigned mNumCenters;
  unsigned mSizeThreshold;
  double mImbalanceTolerance;
  bool mAddSoftMatchingType;
  bool mAddGappedMatchingType;
  bool mAddStemInsensitiveMatchingType;
  bool mAddPPInsensitiveMatchingType;
  string mWeightFileName;
  bool mOutputWeightAnnotation;
  string mWhiteListFileName;
  string mBlackListFileName;
  unsigned mMaxIntersectionSize;
  int mDebug;
};


class MainProcessClass{
protected:
  NSPDK_FeatureGenerator* pmFeatureGenerator;
  NSPDK_FeatureGenerator* pmAbstractFeatureGenerator;
  NSPDK_FeatureGenerator* pmSoftFeatureGenerator;
  NSPDK_FeatureGenerator* pmGappedFeatureGenerator;
  bool mAddAbstractMatchingType;
  bool mAddSoftMatchingType;
  bool mAddGappedMatchingType;
  bool mAddStemInsensitiveMatchingType;
  bool mAddPPInsensitiveMatchingType;
  unsigned mNumHashFunctions;
  unsigned mHashFactor;
  unsigned mNeighborhoodSize;
  unsigned mSampleSize;
  double mNonRedundantFilter;
  string mMode;
  bool mCaching;
  bool mTrueSort;
  double mEccessNeighbourSizeFactor;
  vector<SVector> mDataset;  
  vector<multimap<unsigned,unsigned> > mBinDataStructure;
  map<pair<unsigned,unsigned> , double> mKernelMap;
  multimap <unsigned,unsigned> mInvertedIndex;
  map<unsigned,vector<unsigned> > mSignatureMap;
  bool mOutputHashEncoding;
  string mWhiteListFileName;
  string mBlackListFileName;
  double mAlpha;
  vector<double> mApproximateDensityMap;
  vector<double> mTrueDensityMap;
  map<unsigned,vector<unsigned> > mApproximateNeighborhoodMap;
  vector<unsigned> mIdMap;
 
public:
  MainProcessClass(NSPDK_FeatureGenerator* paFeatureGenerator, 
		   NSPDK_FeatureGenerator* paAbstractFeatureGenerator, 
		   NSPDK_FeatureGenerator* paSoftFeatureGenerator, 
		   NSPDK_FeatureGenerator* paGappedFeatureGenerator, 
		   bool aAddAbstractMatchingType, 
		   bool aAddSoftMatchingType,
		   bool aAddGappedMatchingType,
		   bool aAddStemInsensitiveMatchingType,
		   bool aAddPPInsensitiveMatchingType,
		   unsigned aNumHashFunctions, 
		   unsigned aHashFactor,
		   unsigned aNeighborhoodSize, 
		   unsigned aSampleSize, 
		   double aNonRedundantFilter, 
		   string aMode, 
		   bool aCaching,
		   bool aTrueSort,
		   double aEccessNeighbourSizeFactor,
		   bool aOutputHashEncoding,
		   string aWhiteListFileName,
		   string aBlackListFileName):
    pmFeatureGenerator(paFeatureGenerator), 
    pmAbstractFeatureGenerator(paAbstractFeatureGenerator), 
    pmSoftFeatureGenerator(paSoftFeatureGenerator), 
    pmGappedFeatureGenerator(paGappedFeatureGenerator), 
    mAddAbstractMatchingType(aAddAbstractMatchingType),
    mAddSoftMatchingType(aAddSoftMatchingType),
    mAddGappedMatchingType(aAddGappedMatchingType),
    mAddStemInsensitiveMatchingType(aAddStemInsensitiveMatchingType),
    mAddPPInsensitiveMatchingType(aAddPPInsensitiveMatchingType),
    mNumHashFunctions(aNumHashFunctions),
    mHashFactor(aHashFactor),
    mNeighborhoodSize(aNeighborhoodSize),
    mSampleSize(aSampleSize),
    mNonRedundantFilter(aNonRedundantFilter),
    mMode(aMode),
    mCaching(aCaching),
    mTrueSort(aTrueSort),
    mEccessNeighbourSizeFactor(aEccessNeighbourSizeFactor),
    mOutputHashEncoding(aOutputHashEncoding),
    mWhiteListFileName(aWhiteListFileName),
    mBlackListFileName(aBlackListFileName){}

  void SetMode(string& aMode){mMode=aMode;}
  void Generate(const GraphClass& aG, SVector& oX){
    //create base graph features
    pmFeatureGenerator->generate_feature_vector(aG,oX);

    if (mAddAbstractMatchingType){    
      GraphClass abstract_g(aG);
      //remove vertex label information to make abstract graph i.e. graph based on pure connectivity information
      for (unsigned i=0;i<abstract_g.VertexSize();++i){
	abstract_g.SetVertexSymbolicAttributeList(i,0,"0");
      }
      SVector abstract_x;
      pmAbstractFeatureGenerator->generate_feature_vector(abstract_g,abstract_x);
      //add features to original ones
      oX.add(abstract_x);
    }

    if (mAddStemInsensitiveMatchingType){    
      const string stem_edge_label="s";//NOTE:this is the *hardcoded* label for stem edges!
      GraphClass stem_g(aG);
      for (unsigned u=0;u<stem_g.VertexSize();u++){
	//remove label if at least one incident edge has label=aEdgeLabel
	vector<unsigned> adj=stem_g.GetVertexAdjacentList(u);
	for (unsigned i=0;i<adj.size();++i){
	  if (stem_g.GetEdgeLabel(u,adj[i])==stem_edge_label) {
	    stem_g.SetVertexLabel(u,"S");
	    break;
	  }
	}
      }
      SVector stem_x;
      pmFeatureGenerator->generate_feature_vector(stem_g,stem_x);
      //add features to original ones
      oX.add(stem_x);
    }

    if (mAddPPInsensitiveMatchingType){    
      GraphClass pp_g(aG);
      //Substitute labels with code: R=A|G Y=C,T,U
      for (unsigned i=0;i<pp_g.VertexSize();++i){
	string original_code=aG.GetVertexSymbolicAttributeList(i,0);
	if (original_code=="A" || original_code=="G")
	  pp_g.SetVertexSymbolicAttributeList(i,0,"R");
	if (original_code=="C" || original_code=="U" || original_code=="T")
	  pp_g.SetVertexSymbolicAttributeList(i,0,"Y");
      }
      SVector pp_x;
      pmFeatureGenerator->generate_feature_vector(pp_g,pp_x);
      //add features to original ones
      oX.add(pp_x);
    }

    if (mAddSoftMatchingType){
      SVector soft_x;
      pmSoftFeatureGenerator->generate_feature_vector(aG,soft_x);
      //add features to original ones
      oX.add(soft_x);
    }    

    if (mAddGappedMatchingType){
      SVector gapped_x;
      pmGappedFeatureGenerator->generate_feature_vector(aG,gapped_x);
      //add features to original ones
      oX.add(gapped_x);
    }    

    //re-normalize the overall feature representation
    double norm=dot(oX,oX);
    oX.scale(1/norm);
  }

  void InputStringList(const string& aFileName, vector<string>& oStringList){
    cout<<"Reading "<<aFileName<<endl;
    ifstream fin;
    fin.open(aFileName.c_str());
    if (!fin) throw range_error("Cannot open file:"+aFileName);
    ProgressBar progress_bar;
    while (!fin.eof() && fin.good()){
      string line;
      getline(fin,line);
      stringstream ss;
      ss<<line<<endl;
      while (!ss.eof() && ss.good()){
	string target;
	ss>>target;
	if (target!=""){	    
	  oStringList.push_back(target);
	  progress_bar.Count();
	}
      }
    }
    fin.close();
  }

  void InputIntList(const string& aFileName, vector<int>& oList){
    cout<<"Reading "<<aFileName<<endl;
    ifstream fin;
    fin.open(aFileName.c_str());
    if (!fin) throw range_error("Cannot open file:"+aFileName);
    ProgressBar progress_bar;
    while (!fin.eof() && fin.good()){
      string line;
      getline(fin,line);
      if (line!=""){
	stringstream ss;
	ss<<line;
	int target;
	ss>>target;
	oList.push_back(target);
	progress_bar.Count();
      }
    }
    fin.close();
  }

  void InputWeightList(const string& aFileName, SVector& oW){
    cout<<"Reading "<<aFileName<<endl;
    ifstream fin;
    fin.open(aFileName.c_str());
    if (!fin) throw range_error("Cannot open file:"+aFileName);
    ProgressBar progress_bar;
    while (!fin.eof() && fin.good()){
      string line;
      getline(fin,line);
      stringstream ss;
      ss<<line<<endl;
      while (!ss.eof() && ss.good()){
	unsigned index;
	double value;
	ss>>index>>value;
	if (ss.good()){	    
	  oW.set(index,value);
	  progress_bar.Count();
	}
      }
    }
    fin.close();
  }

  void LoadAndAnnotate(const string& aInputFileName,const string& aInputWeightFileName){
    SVector w;
    InputWeightList(aInputWeightFileName,w);

    string ofname=aInputFileName+".annotated";
    ofstream ofs_f(ofname.c_str());

    cout<<"Reading gspan data and computing annotation"<<endl;
    ifstream fin;
    fin.open(aInputFileName.c_str());
    if (!fin) throw range_error("Cannot open file:"+aInputFileName);
    ProgressBar progress_bar;
    unsigned counter=0;
    vector<unsigned> viewpoint;viewpoint.push_back(0);
    while (!fin.eof()){
      GraphClass G;
      SetGraphFromFileGSPAN(fin, G, mMode);
      pmFeatureGenerator->Clear();
      SVector tmpx;
      pmFeatureGenerator->generate_feature_vector(G,tmpx);
      GraphClass annotated_g(G);
      //for all vertices compute importance as dot product of features of single vertex with weight vector
      for (unsigned i=0;i<G.VertexSize();i++){
	//select only one vertex
	SVector x;
	viewpoint[0]=i;
	pmFeatureGenerator->generate_feature_vector(G,x,viewpoint);
	double score=dot(x,w);
	annotated_g.SetVertexNumericAttributeList(i,0,score);
      }
      pmFeatureGenerator->Clear();
      //output annotated graph
      string gspan_filename=aInputFileName+"_"+stream_cast<string>(counter+1);//NOTE:+1 to be consistent with id numbering starting from 1 rather than 0
      annotated_g.SaveAsGspanFile(gspan_filename);
      progress_bar.Count();
      counter++;
    }
    fin.close();
  }

  void InputAndOutput(const string& aInputFileName){
    Load(aInputFileName,true);
  }

  void Input(const string& aInputFileName){
    Load(aInputFileName,false);
  }

  void Load(const string& aInputFileName, bool aDirectProcess){
    ofstream ofs_f;
    ofstream ofs_fb;
    if (aDirectProcess){
      string ofname=aInputFileName+".feature";
      ofs_f.open(ofname.c_str());
      ofname=aInputFileName+".feature_bin";
      ofs_fb.open(ofname.c_str());
    }

    //read white list
    vector<int> select_list;
    if (mWhiteListFileName!=""){
      InputIntList(mWhiteListFileName,select_list);
    } 
    //read black list
    if (mBlackListFileName!=""){
      InputIntList(mBlackListFileName,select_list);
    }
    set<int> select_list_set;
    select_list_set.insert(select_list.begin(),select_list.end());

    cout<<"Reading gspan data and computing features"<<endl;
    ifstream fin;
    fin.open(aInputFileName.c_str());
    if (!fin) throw range_error("Cannot open file:"+aInputFileName);
    ProgressBar progress_bar;
    int counter=1;
    while (!fin.eof()){
      GraphClass G;
      SetGraphFromFileGSPAN(fin, G, mMode);
      SVector x;

      //only if counter id is consistent with white and black list (if they have been specified) then accept the instance 
      bool accept_flag=true;
      if (mWhiteListFileName!=""){
	if (select_list_set.count(counter)>0) accept_flag=true;
	else accept_flag=false;
      }
      if (mBlackListFileName!=""){
	if (select_list_set.count(counter)>0) accept_flag=false;
	else accept_flag=true;
      }

      if (accept_flag==true){
	Generate(G,x);
	if (aDirectProcess){
	  ofs_f<<x;
	  x.save(ofs_fb);
	} else { 
	  mDataset.push_back(x);
	}
	mIdMap.push_back(counter);
      }
      progress_bar.Count();
      counter++;
    }
    fin.close();
  }

  void InputSparse(const string& aInputFileName, string aMode){
    ifstream fin;
    fin.open(aInputFileName.c_str());
    if (!fin) throw range_error("Cannot open file:"+aInputFileName);
    InputSparse(fin,aMode,mDataset);
    fin.close();
  }

  void InputSparse(const string& aInputFileName, string aMode, vector<SVector>& oDataset){
    ifstream fin;
    fin.open(aInputFileName.c_str());
    if (!fin) throw range_error("Cannot open file:"+aInputFileName);
    InputSparse(fin,aMode,oDataset);
    fin.close();
  }

  void InputSparse(ifstream& aFin, string aMode, vector<SVector>& oDataset){
    //read white list
    vector<int> select_list;
    if (mWhiteListFileName!=""){
      InputIntList(mWhiteListFileName,select_list);
    } 
    //read black list
    if (mBlackListFileName!=""){
      InputIntList(mBlackListFileName,select_list);
    }
    set<int> select_list_set;
    select_list_set.insert(select_list.begin(),select_list.end());

    cout<<"Reading file in "<<aMode<<" mode"<<endl;
    int counter=1;
    ProgressBar progress_bar;
    while (!aFin.eof() && aFin.good()){
      SVector x;
      if (aMode=="binary") x.load(aFin);
      else ParseASCIILine2Vector(aFin,x);
      if (InstanceIsValid(x)==true) {
	//only if counter id is consistent with white and black list (if they have been specified) then accept the instance 
	bool accept_flag=true;
	if (mWhiteListFileName!=""){
	  if (select_list_set.count(counter)>0) accept_flag=true;
	  else accept_flag=false;
	}
	if (mBlackListFileName!=""){
	  if (select_list_set.count(counter)>0) accept_flag=false;
	  else accept_flag=true;
	}
	if (accept_flag==true){
	  oDataset.push_back(x);
	  mIdMap.push_back(counter);
	}
	progress_bar.Count();
	counter++;
      } else {}//discard non valid instances
    }
  }

  inline void ParseASCIILine2Vector(ifstream& aFin, SVector& aX){
    string line;
    getline(aFin,line);
    if (line=="") return;
    stringstream ss;
    ss<<line<<endl;
    while (!ss.eof() && ss.good()){
      string key_value;
      ss>>key_value;
      size_t limit=key_value.find_first_of(":", 0);
      if (limit!=string::npos){//if the delimiter ':' is found then proceed
	string key=key_value.substr(0, limit);
	string value=key_value.substr(limit+1, key_value.size());
	unsigned key_int=stream_cast<unsigned>(key);
	double val_real=stream_cast<double>(value);
	aX.set(key_int,val_real);
      }
    }
  }

  inline  bool InstanceIsValid(SVector& aX){
    bool is_valid=false;
    //if (dot(aX,aX)>0) is_valid=true;
    if (aX.sparse_size()>0) is_valid=true;
    return is_valid;
  }

  void SetGraphFromFileGSPAN(istream& in, GraphClass& oG, string aMode=""){
    map<string,int> index_map_nominal_to_real;
    string line;
    getline(in,line);
    assert(line[0]=='t');//first line must have as first char a 't'
    while(!in.eof() && in.good() && in.peek()!='t' && getline(in,line)){//read until next 't' or end of file
      stringstream ss;
      ss<<line<<endl;
      char code;
      ss>>code;
      if (code=='v'){
	//extract vertex id and make map nominal_id -> real_id
	string nominal_vertex_index;
	ss>>nominal_vertex_index;
	unsigned real_vertex_index=oG.InsertVertex();
	index_map_nominal_to_real[nominal_vertex_index]=real_vertex_index;
	//label
	vector<string> vertex_symbolic_attribute_list;
	string label;
	ss>>label;
	vertex_symbolic_attribute_list.push_back(label);
	if (aMode=="read_absolute_vertex_id"){
	  string absolute_vertex_id;
	  ss>>absolute_vertex_id;
	  vertex_symbolic_attribute_list.push_back(absolute_vertex_id);
	}
	oG.SetVertexSymbolicAttributeList(real_vertex_index, vertex_symbolic_attribute_list);
	if (aMode=="3D"){
	  double x,y,z;
	  ss>>x>>y>>z;
	  vector<double> vertex_numeric_attribute_list;
	  vertex_numeric_attribute_list.push_back(x);
	  vertex_numeric_attribute_list.push_back(y);
	  vertex_numeric_attribute_list.push_back(z);
	  oG.SetVertexNumericAttributeList(real_vertex_index, vertex_numeric_attribute_list);
	}

	//status
	vector<bool> status; 
	 
	status.push_back(true);//kernel point
	status.push_back(true);//kind
	if (aMode=="viewpoint"){//if mode=viewpoint then assume extra vertex label is \in {0,1} representing viewpoint flag
	  bool vp;
	  ss>>vp;
	  status.push_back(vp);//viewpoint
	}	  
	else status.push_back(true);//viewpoint
	status.push_back(false);//dead
	status.push_back(false);//forbidden distance
	status.push_back(false);//forbidden neighborhood
	oG.SetVertexStatusAttributeList(real_vertex_index, status);
      } else if (code=='e'){
	//extract src and dest vertex id 
	string nominal_src_index,nominal_dest_index;
	string label;
	double weight;
	ss>>nominal_src_index>>nominal_dest_index>>label>>weight;
	unsigned real_src_index=index_map_nominal_to_real[nominal_src_index];
	unsigned real_dest_index=index_map_nominal_to_real[nominal_dest_index];
	unsigned edge_index=oG.InsertEdge(real_src_index,real_dest_index);
	unsigned reverse_edge_index=oG.InsertEdge(real_dest_index,real_src_index);
	vector<string> edge_symbolic_attribute_list;
	edge_symbolic_attribute_list.push_back(label);
	vector<double> edge_numeric_attribute_list;
	edge_numeric_attribute_list.push_back(weight);
	oG.SetEdgeSymbolicAttributeList(edge_index,edge_symbolic_attribute_list);
	oG.SetEdgeSymbolicAttributeList(reverse_edge_index,edge_symbolic_attribute_list);
	oG.SetEdgeNumericAttributeList(edge_index,edge_numeric_attribute_list);
	oG.SetEdgeNumericAttributeList(reverse_edge_index,edge_numeric_attribute_list);
      } else {}//NOTE: ignore other markers
    }
  }

  void ComputeBinDataStructure(unsigned aNumHashFunctions){
    string ofname="hash_encoding";
    ofstream of(ofname.c_str());
        
    cout<<"Computing bin data structure..."<<endl;
    ProgressBar progress_bar;

    //init structure
    mBinDataStructure.clear();
    for (unsigned k=0;k<aNumHashFunctions;++k)
      mBinDataStructure.push_back(multimap<unsigned,unsigned>());
    
    //fill structure
    for (unsigned i=0;i<mDataset.size();++i){
      vector<unsigned> min_list=ComputeHashSignature(i,aNumHashFunctions);
      if (mOutputHashEncoding) {for(unsigned j=0;j<min_list.size();j++) of<<min_list[j]<<" ";of<<endl;}
      for (unsigned k=0;k<aNumHashFunctions;++k)
	mBinDataStructure[k].insert(make_pair(min_list[k],i));
      progress_bar.Count();
    }
  }

  inline  vector<unsigned> ComputeHashSignature(unsigned aID, unsigned aNumHashFunctions){
    if (mSignatureMap.count(aID)>0) return mSignatureMap[aID];
    else { 
      vector<unsigned> signature=ComputeHashSignature(mDataset[aID],aNumHashFunctions);
      mSignatureMap[aID]=signature;
      return signature;
    }
  }

  inline  vector<unsigned> ComputeHashSignature(SVector& aX, unsigned aNumHashFunctions){
    unsigned effective_num_hash_functions=aNumHashFunctions*mHashFactor;
    //const double A=sqrt(2)-1;
    const unsigned MAXUNSIGNED=2<<30;
    vector<unsigned> signature;
    //special case for first hash function just consider first feature id
    //signature.push_back(mDataset[aID].extract_component(0).first);
    //for all other hash functions rehash and select the minimum
    for (unsigned k=0;k<effective_num_hash_functions;++k)
      signature.push_back(MAXUNSIGNED);
    unsigned size=(unsigned)aX.sparse_size();
    for (unsigned f=0;f<size;++f){    
      unsigned hash_id=aX.extract_component(f).first;
      if (hash_id==0) throw range_error("Error: Feature ID = 0. Feature ID  has to be strictly > 0");
      for (unsigned k=0;k<effective_num_hash_functions;++k){
	unsigned new_hash=IntHash(hash_id,MAXUNSIGNED,k);
	if (signature[k]>new_hash) signature[k]=new_hash;
      }
    }
    //compact signature
    vector<unsigned> compact_signature;
    for (unsigned i=0;i<signature.size();i=i+mHashFactor){
      unsigned new_hash=0;
      for(unsigned j=0;j<mHashFactor;j++)
	new_hash+=signature[i+j];
      compact_signature.push_back(new_hash);
    }
    return compact_signature;
  }

  void OutputBinDataStructureStatistics()const{
    cout<<"Bin size statistics: ";
    VectorClass stats;
    for (unsigned k=0;k<mBinDataStructure.size();++k){
      for (multimap<unsigned,unsigned>::const_iterator it=mBinDataStructure[k].begin();it!=mBinDataStructure[k].end();++it){
	unsigned key=it->first;
	unsigned bin_size=mBinDataStructure[k].count(key);
	stats.PushBack(bin_size);
      }
    }
    stats.OutputStatistics(cout);
    cout<<endl;
  }

  void ComputeInvertedFeatureInstanceIndex(){
    cout<<"Computing inverted data structure..."<<endl;
    ProgressBar progress_bar;
    for (unsigned i=0;i<mDataset.size();++i){
      for (unsigned k=0;k<(unsigned)mDataset[i].sparse_size();++k){
	unsigned hash_id=mDataset[i].extract_component(k).first;
	mInvertedIndex.insert(make_pair(hash_id,i));
      }
      progress_bar.Count();
    }
  }

  vector<unsigned> ComputeFastNeighborhood(unsigned aID){
    map<unsigned,int> neighborhood;
    for (unsigned k=0;k<(unsigned)mDataset[aID].sparse_size();++k){
      unsigned hash_id=mDataset[aID].extract_component(k).first;
      //return equal range wrt hash_id
      pair<multimap<unsigned,unsigned>::iterator,multimap<unsigned,unsigned>::iterator> erange=mInvertedIndex.equal_range(hash_id);
      //insert into neighborhood
      for (multimap<unsigned,unsigned>::iterator it=erange.first;it!=erange.second;++it){
	unsigned instance_id=it->second;
	//count number of occurences
	if (neighborhood.count(instance_id)>0) neighborhood[instance_id]++;
	else neighborhood[instance_id]=1;
      }
    }
    return TrimNeighborhood(neighborhood);
  }

  vector<unsigned> ComputeApproximateNeighborhood(unsigned aID){
    if (mApproximateNeighborhoodMap.count(aID)==0){ 
      vector<unsigned> hash_signature=ComputeHashSignature(aID,mNumHashFunctions);
      vector<unsigned> neighborhood=ComputeApproximateNeighborhood(hash_signature);
      //select neighborhood under true similarity function on the subset of indiced returned by ComputeApproximateNeighborhood
      vector<unsigned> true_neighborhood=ComputeTrueSubNeighborhood(aID,neighborhood);
      mApproximateNeighborhoodMap[aID]=true_neighborhood;
    }
    return mApproximateNeighborhoodMap[aID];
  }

  vector<unsigned> ComputeApproximateNeighborhood(const vector<unsigned>& aInstanceSignature){
    map<unsigned,int> neighborhood;
    vector<pair<unsigned, double> > vec;
    for (unsigned k=0;k<mNumHashFunctions;++k){
      unsigned hash_id=aInstanceSignature[k];
      //return equal range wrt hash_id
      pair<multimap<unsigned,unsigned>::iterator,multimap<unsigned,unsigned>::iterator> erange=mBinDataStructure[k].equal_range(hash_id);
      //fill neighborhood set counting number of occurrences
      for (multimap<unsigned,unsigned>::iterator it=erange.first;it!=erange.second;++it){
	unsigned instance_id=it->second;
	if (neighborhood.count(instance_id)>0) neighborhood[instance_id]++;
	else neighborhood[instance_id]=1;
      }
    }
    return TrimNeighborhood(neighborhood);
  }

  vector<unsigned> ComputeTrueSubNeighborhood(unsigned aID, vector<unsigned>& aApproximateNeighborhoodList){
    vector<pair<double, unsigned> > rank_list;
    for (unsigned i=0;i<aApproximateNeighborhoodList.size();++i){
      unsigned id_neighbour=aApproximateNeighborhoodList[i];
      double k=Kernel(aID,id_neighbour);
      rank_list.push_back(make_pair(-k,id_neighbour));
    }
    unsigned effective_size=min((unsigned)rank_list.size(),mNeighborhoodSize);
    sort(rank_list.begin(),rank_list.end());
    vector<unsigned> neighbor;
    for (unsigned j=0;j<effective_size;j++){
      neighbor.push_back(rank_list[j].second);
    }
    return neighbor;
  }

  vector<unsigned> ComputeTrueNeighborhood(const SVector& aX, unsigned aSize){
    vector<pair<double, unsigned> > rank_list;
    for (unsigned i=0;i<mDataset.size();++i){
      double k=dot(aX,mDataset[i]);
      rank_list.push_back(make_pair(-k,i));
    }
    unsigned effective_size=min((unsigned)rank_list.size(),aSize);
    partial_sort(rank_list.begin(),rank_list.begin()+effective_size,rank_list.end());
    vector<unsigned> neighbor;
    for (unsigned j=0;j<effective_size;j++)
      neighbor.push_back(rank_list[j].second);
    return neighbor;
  }

  vector<unsigned> TrimNeighborhood(map<unsigned,int>& aNeighborhood){
    const int MIN_BINS_IN_COMMON=2;//Minimum number of bins that two instances have to have in common in order to be considered similar
    //given a list of neighbours with an associated occurences count, return only a fraction of the highest count ones
    vector<unsigned> neighborhood_list;
    if (mEccessNeighbourSizeFactor>0){
      //sort by num occurences
      vector<pair<int,unsigned> > count_list;
      for (map<unsigned,int>::const_iterator it=aNeighborhood.begin();it!=aNeighborhood.end();++it){
	unsigned id=it->first;
	int count=it->second;
	if (count>=MIN_BINS_IN_COMMON)//NOTE: consider instances that have at least MIN_BINS_IN_COMMON 
	  count_list.push_back(make_pair(-count,id));//NOTE:-count to sort from highest to lowest
      }
      unsigned effective_size=min((unsigned)count_list.size(),(unsigned)(mEccessNeighbourSizeFactor*mNeighborhoodSize));
      sort(count_list.begin(),count_list.end());
      for (unsigned i=0;i<effective_size;++i)
	neighborhood_list.push_back(count_list[i].second);
    } else {//if mEccessNeighbourSizeFactor==0 then just consider all the ids in the approximate neighborhood 
      for (map<unsigned,int>::const_iterator it=aNeighborhood.begin();it!=aNeighborhood.end();++it){
	neighborhood_list.push_back(it->first);
      }
    }
    return neighborhood_list;
  }

  double ComputeApproximateDensity(unsigned aID){
    double density=0;
    if (mApproximateDensityMap[aID]==-1){
      vector<unsigned> approximate_neighborhood=ComputeApproximateNeighborhood(aID);
      //compute kernel pairs between aID and all elements in aApproximateNeighborhood
      for (unsigned j=0;j<approximate_neighborhood.size();j++){
	unsigned u=aID;
	unsigned v=approximate_neighborhood[j];
	if (u!=v){
	  double k_uv=Kernel(u,v);
	  density+=k_uv;
	}
      }
      density=density/(approximate_neighborhood.size()-1);
      mApproximateDensityMap[aID]=density;
    } else {
      density=mApproximateDensityMap[aID];
    }
    return density;
  }

  double ComputeTrueDensity(unsigned aID){
    double density=0;
    if (mTrueDensityMap[aID]==-1){
      vector<pair<double,unsigned> > sim_list;
      for (unsigned j=0;j<mDataset.size();j++){
	if (aID!=j){
	  double k_ij=Kernel(aID,j);
	  density+=k_ij;
	}
      }
      density=density/(mDataset.size()-1);
      mTrueDensityMap[aID]=density;
    } else
      density=mTrueDensityMap[aID]; 
    return density;
  }

  vector<unsigned> ComputeRandomSampleHighDensityInstancesGreedy(unsigned aSampleSize){
    cout<<"Starting iterative density center sampling.."<<endl;////
    ProgressBar progress_bar;
    
    set<unsigned> best_workingset;
    double best_similarity=0;
    double best_density=0;
    const unsigned NUM_REPEATS=100;
    for (unsigned repeat=0;repeat<NUM_REPEATS;repeat++){
      //create initial random center list
      unsigned workingset_limit=0;
      vector<unsigned> workingset;
      for (unsigned i=0;i<mDataset.size();++i)
	workingset.push_back(i);
      for (unsigned i=0;i<mDataset.size();++i){
	unsigned j=randomUnsigned(mDataset.size());
	swap(workingset[i],workingset[j]);
      } 

      //keep priority queue of -density for working set
      priority_queue< pair<double, unsigned> > density_pq;   
      set<unsigned> active_workingset;
      for (unsigned i=0;i<aSampleSize;i++){
	unsigned id=workingset[i];
	active_workingset.insert(id);
	workingset_limit++;
	double density=ComputeApproximateDensity(id);  
	density_pq.push(make_pair(-density,id));//NOTE: use -density to prioritize from the least dense
	progress_bar.Count();
      }
      double avg_similarity=ComputeAverageSimilarity(active_workingset);

      for(unsigned i=workingset_limit;i<mDataset.size();++i){
	//take instances in random order
	unsigned id=workingset[i];
	double density=ComputeApproximateDensity(id);
	//if density is higher than current minimum density then proceed further...
	double min_density=-density_pq.top().first;
	unsigned min_density_id=density_pq.top().second;
	if (density>min_density){
	  set<unsigned> tentative_set;
	  for (set<unsigned>::iterator it=active_workingset.begin();it!=active_workingset.end();++it){
	    unsigned current_id=*it;
	    if (current_id!=min_density_id) tentative_set.insert(current_id);
	  }
	  tentative_set.insert(id);
	  double tentative_avg_similarity=ComputeAverageSimilarity(tentative_set);
	  bool random_choice=(random01()<(1-(repeat+1)/NUM_REPEATS))?true:false;
	  if (random_choice==true || tentative_avg_similarity<avg_similarity) {//if similarity is decreased then accept instance and swap out the least dense instance from the active set
	    active_workingset=tentative_set;
	    avg_similarity=tentative_avg_similarity;
	    density_pq.pop();
	    density_pq.push(make_pair(density,id));
	  }
	} 
	progress_bar.Count();
      }
      double active_similarity=ComputeAverageSimilarity(active_workingset);
      double active_density=ComputeAverageDensity(active_workingset);
      if (repeat==0 || ( active_similarity < best_similarity && active_density > best_density)){
	cout<<" sim:"<<active_similarity<<" dens:"<<active_density;
	best_workingset=active_workingset;
	best_similarity=active_similarity;
	best_density=active_density;
      }
    }
    vector<unsigned> result;
    for (set<unsigned>::iterator it=best_workingset.begin();it!=best_workingset.end();++it)
      result.push_back(*it);
    return result;
  }

  double ComputeAverageDensity(set<unsigned>& aSet){
    double density=0;
    for (set<unsigned>::iterator it=aSet.begin();it!=aSet.end();++it){
      unsigned id = *it;
      density+=ComputeApproximateDensity(id);  
    }
    return density;
  }

  vector<unsigned> ComputeMinimallyOverlappingHighDensityCenterList(unsigned aSampleSize, unsigned aMaxIntersectionSize){
    //compute density for each instance in dataset
    cout<<"Computing approximate density information for each instance."<<endl;////
    vector<pair<double,unsigned> > density_list;
    {
      ProgressBar progress_bar;
      for (unsigned i=0;i<mDataset.size();++i){
	double density=ComputeApproximateDensity(i);
	density_list.push_back(make_pair(-density,i));
	progress_bar.Count();
	
      }
    }
    sort(density_list.begin(),density_list.end());
    vector<unsigned> result;
    set<unsigned> active_neighborhood;
    {
      ProgressBar progress_bar;
      cout<<"Computing minimally overlapping high density center list for up to "<<aSampleSize<<" centers."<<endl;////
      for (unsigned i=0;i<density_list.size() && result.size()<aSampleSize;i++){
	unsigned id=density_list[i].second;
	vector<unsigned> neighborhood=ComputeApproximateNeighborhood(id);
	set<unsigned> neighborhood_set;
	neighborhood_set.insert(neighborhood.begin(),neighborhood.end());
	set<unsigned> intersection;
	set_intersection( active_neighborhood.begin(), active_neighborhood.end(), neighborhood_set.begin(), neighborhood_set.end(), inserter(intersection,intersection.begin()) );
	if (i==0 || intersection.size()<=aMaxIntersectionSize){//if the intersection between the neighborhood of the current center and the union of all active neighborhoods is less than a defined constant (eg. 0) then accept the new center in the active set
	  active_neighborhood.insert(neighborhood.begin(),neighborhood.end());
	  result.push_back(id);
	  progress_bar.Count();
	} 
      }
    }
    return result;
  }

  vector<unsigned> ComputeRandomSampleHighDensityInstances(unsigned aSampleSize){
    //compute density for each instance in dataset
    cout<<"Computing approximate density information for each instance."<<endl;////
    vector<double> density_list;
    double max_density_val=ComputeApproximateDensity(0);
    unsigned max_density_id=0;
    double min_density_val=max_density_val;
    unsigned min_density_id=0;
    unsigned non_zero_density_count=0;
    {
      ProgressBar progress_bar;
      progress_bar.Count();
      for (unsigned i=1;i<mDataset.size();++i){
	double density=ComputeApproximateDensity(i);
	density_list.push_back(density);
	progress_bar.Count();
	if (density>0){
	  non_zero_density_count++;
	  if (max_density_val<density) {
	    max_density_val=density;
	    max_density_id=i;
	  }
	  if (min_density_val>density) {
	    min_density_val=density;
	    min_density_id=i;
	  }
	}
      }
    }

    VectorClass density_stats(density_list);
    cout<<"Density statistics: ";
    density_stats.OutputStatistics(cout);
    cout<<endl;

    //sample tentative centers
    cout<<"Sampling high density instances."<<endl;////
    set<unsigned> center_sample_set;
    center_sample_set.insert(max_density_id);
    const unsigned SAMPLE_SIZE_FACTOR=4;//extract a set of tentative high density instances of size SAMPLE_SIZE_FACTOR*aSampleSize 
    const int MAX_NUM_DATASET_SCANS=5000;
    int counter=MAX_NUM_DATASET_SCANS;
    unsigned effective_size_center_sample=min(non_zero_density_count,aSampleSize*SAMPLE_SIZE_FACTOR);
    {
      ProgressBar progress_bar;
      while (counter>0 && center_sample_set.size()<effective_size_center_sample){
	for (unsigned i=0;i<density_list.size();++i){
	  double d=(density_list[i]-min_density_val)/(max_density_val-min_density_val);//Normalize density between 1=max density and 0=min density
	  double r=random01();
	  if (d > r) {//bias current instance selection so that it is not too similar to previously selected instances
	    double avg_sim=ComputeAverageSimilarityToSet(i,center_sample_set);
	    double sim=pow(avg_sim,mNonRedundantFilter);
	    double r=random01();
	    if ( sim < r ){
	      center_sample_set.insert(i);
	      density_list[i]=0;
	      progress_bar.Count();
	    } 
	  }
	}
	counter--;
      }
    }

    //compute the true density for the selected sample and retain only the highest density instances
    cout<<"Sparsifying phase for selected "<<center_sample_set.size()<<" centers."<<endl;////
    vector<pair<double,unsigned> > sim_list;
    {
      ProgressBar progress_bar;
      for (set<unsigned>::iterator it=center_sample_set.begin();it!=center_sample_set.end();++it){
	unsigned id=(*it);
	double value=-ComputeApproximateDensity(id);//NOTE:-density to sort starting from highest true density
	sim_list.push_back(make_pair(value,id));
	progress_bar.Count();
      }
    }
    sort(sim_list.begin(),sim_list.end());
    //sample one every SAMPLE_SIZE_FACTOR entries to guarantee centers diversity
    vector<unsigned> sample;
    for (unsigned k=0;k<sim_list.size() && sample.size()<aSampleSize;++k)
      if (k%SAMPLE_SIZE_FACTOR==0) sample.push_back(sim_list[k].second);
    return sample;
  }

  double ComputeAverageSimilarityToSet(unsigned aID, set<unsigned>& aSet){
    double avg_sim=0;
    for (set<unsigned>::iterator it=aSet.begin();it!=aSet.end();++it){
      unsigned id=(*it);
      avg_sim+=Kernel(aID,id);
    }
    return avg_sim/aSet.size();
  }

  double ComputeAverageSimilarity(set<unsigned>& aSet){
    double avg_sim=0;
    for (set<unsigned>::iterator it=aSet.begin();it!=aSet.end();++it){
      unsigned id=(*it);
      avg_sim+=ComputeAverageSimilarityToSet(id,aSet);
    }
    return avg_sim/aSet.size();
  }

  vector<unsigned> ComputeNeighborhoodRanking(unsigned aID){
    vector<pair<double,unsigned> > sim_list;
    for (unsigned i=0;i<mDataset.size();++i){
      if (i!=aID) {
	double k=Kernel(aID,i);
	sim_list.push_back(make_pair(-k,i));//note: use -k to sort in decreasing order
      }
    }
    //sort and take most similar <effective_size>
    sort(sim_list.begin(),sim_list.end());
    vector<unsigned> neighborhood;
    for (unsigned t=0;t<sim_list.size();++t)
      neighborhood.push_back(sim_list[t].second);
    return neighborhood;
  }

  void OutputCluster(ostream& out, string aOutputCluster, unsigned aMaxIntersectionSize){
    //initialize density cache
    mApproximateDensityMap.clear();
    mTrueDensityMap.clear();
    for (unsigned i=0;i<mDataset.size();++i){
      mApproximateDensityMap.push_back(-1);
      mTrueDensityMap.push_back(-1);
    }

    vector<double> density_list;
    vector<unsigned> density_center_list;
    if (aOutputCluster=="GREEDY") density_center_list=ComputeRandomSampleHighDensityInstancesGreedy(mSampleSize);
    else if (aOutputCluster=="DENSITY") density_center_list=ComputeMinimallyOverlappingHighDensityCenterList(mSampleSize, aMaxIntersectionSize);
    else if (aOutputCluster=="SAMPLE") density_center_list=ComputeRandomSampleHighDensityInstances(mSampleSize);
    else throw range_error("Unknown center mode type:"+aOutputCluster);
    cout<<"Compute (approximate) neighborhood for selected "<<density_center_list.size()<<" cluster centers."<<endl;////
    ProgressBar progress_bar;
    for (unsigned i=0;i<density_center_list.size();++i){
      unsigned id=density_center_list[i];
      vector<unsigned> neighborhood=ComputeApproximateNeighborhood(id);
      density_list.push_back(ComputeApproximateDensity(id));
      progress_bar.Count();
      //unsigned relative_id=id;
      //unsigned absolute_id=mIdMap[relative_id];
      //out<<absolute_id<<"   ";      
      for (unsigned i=0;i<neighborhood.size();i++){
	unsigned relative_id=neighborhood[i];
	unsigned absolute_id=mIdMap[relative_id];
	out<<absolute_id<<" ";       
      }
      out<<endl;
    }


    {
      VectorClass density_stats(density_list);
      cout<<endl<<"Centers density statistics: ";
      density_stats.OutputStatistics(cout);
      cout<<endl;
    }
    {
      VectorClass density_stats(mApproximateDensityMap);
      cout<<endl<<"Global density statistics: ";
      density_stats.OutputStatistics(cout);
      cout<<endl;
    }
  }

  void OutputApproximateKNN(ostream& out){
    cout<<"Compute approximate nearest neighbours."<<endl;////
    ProgressBar progress_bar;
    for (unsigned u=0;u<mDataset.size();++u){
      vector<unsigned> approximate_neighborhood=ComputeApproximateNeighborhood(u);

      if (mTrueSort==true){//sort all instances in neighborhood according to true similarity and take only the mNeighborhoodSize nearest ones
	vector<pair<double,unsigned> > sim_list;
	//compute kernel pairs between aID and all elements in aApproximateNeighborhood
	for (unsigned j=0;j<approximate_neighborhood.size();j++){
	  unsigned v=approximate_neighborhood[j];
	  double k_uv=Kernel(u,v);
	  sim_list.push_back(make_pair(-k_uv,v));//note: use -k to sort in decreasing order
	}
	//sort and take truly most similar 
	sort(sim_list.begin(),sim_list.end());
	unsigned effective_size=min(mNeighborhoodSize,(unsigned)sim_list.size());
	for (unsigned k=0;k<effective_size;++k){
	  out<<sim_list[k].second+1<<" ";//NOTE: numbering starts from 1
	}
	out<<endl;
      } else {//do not sort neighbours according to true similarity, just output them
	for (unsigned t=0;t<approximate_neighborhood.size();++t){
	  out<<approximate_neighborhood[t]+1<<" ";//NOTE: numbering starts from 1
	}
	out<<endl;
      }
      progress_bar.Count();
    }
  }

  void OutputFastKNN(ostream& out){
    cout<<"Compute fast nearest neighbours."<<endl;////
    ProgressBar progress_bar;
    for (unsigned u=0;u<mDataset.size();++u){
      vector<unsigned> neighborhood=ComputeFastNeighborhood(u);
      vector<pair<double,unsigned> > sim_list;
      //compute kernel pairs between aID and all elements in neighborhood
      for (unsigned j=0;j<neighborhood.size();j++){
	unsigned v=neighborhood[j];
	double k_uv=Kernel(u,v);
	sim_list.push_back(make_pair(-k_uv,v));//note: use -k to sort in decreasing order
      }
      //sort and take truly most similar 
      sort(sim_list.begin(),sim_list.end());
      unsigned effective_size=min(mNeighborhoodSize,(unsigned)sim_list.size());
      for (unsigned k=0;k<effective_size;++k){
	out<<sim_list[k].second+1<<" ";//NOTE: numbering starts from 1
      }
      out<<endl;
      progress_bar.Count();
    }
  }

  void OutputTrueKNN(ostream& out){
    cout<<"Compute true nearest neighbours."<<endl;////
    ProgressBar progress_bar;
    for (unsigned u=0;u<mDataset.size();++u){
      vector<pair<double,unsigned> > sim_list;
      //compute kernel pairs between aID and all elements 
      for (unsigned v=0;v<mDataset.size();v++){
	double k_uv=Kernel(u,v);
	sim_list.push_back(make_pair(-k_uv,v));//note: use -k to sort in decreasing order
      }
      //sort and take truly most similar 
      sort(sim_list.begin(),sim_list.end());
      for (unsigned k=0;k<mNeighborhoodSize;++k){
	out<<sim_list[k].second+1<<" ";//NOTE: numbering starts from 1
      }
      out<<endl;
      progress_bar.Count();
    }
  }

  void OutputApproximateKNNPrediction(ostream& out, string& aSparseASCIITestInputFileName, string& aTrainTargetInputFileName){
    cout<<"Compute approximate nearest neighbour prediction of test instances."<<endl;////
    ProgressBar progress_bar;

    //read test instances
    vector<SVector> test_dataset;
    InputSparse(aSparseASCIITestInputFileName, "ascii", test_dataset);

    //read train targets
    vector<string> train_target_list;
    InputStringList(aTrainTargetInputFileName,train_target_list);

    cout<<"Computing k-NN predictions."<<endl;
    //for each test instance
    for (unsigned u=0;u<test_dataset.size();++u){
      //extract signature
      vector<unsigned> hash_signature=ComputeHashSignature(test_dataset[u],mNumHashFunctions);
      //extract knn
      vector<unsigned> approximate_neighborhood=ComputeApproximateNeighborhood(hash_signature);
      string prediction=KNNPredict(approximate_neighborhood,train_target_list);
      out<<prediction<<endl;
      progress_bar.Count();
    }
  }

  void OutputTrueKNNPrediction(ostream& out, string& aSparseASCIITestInputFileName, string& aTrainTargetInputFileName){
    cout<<"Compute true nearest neighbour prediction of test instances."<<endl;////
    ProgressBar progress_bar;

    //read test instances
    vector<SVector> test_dataset;
    InputSparse(aSparseASCIITestInputFileName, "ascii", test_dataset);

    //read train targets
    vector<string> train_target_list;
    InputStringList(aTrainTargetInputFileName,train_target_list);

    cout<<"Computing k-NN predictions on "<<test_dataset.size()<<" instances."<<endl;
    //for each test instance
    for (unsigned u=0;u<test_dataset.size();++u){
      //extract knn
      vector<unsigned> neighborhood=ComputeTrueNeighborhood(test_dataset[u],mNeighborhoodSize);
      //compute majority class
      string prediction=KNNPredict(neighborhood,train_target_list);
      out<<prediction<<endl;
      progress_bar.Count();
    }
  }

  string KNNPredict(const vector<unsigned>& aNeighborhood, const vector<string>& aTargetList)const{
    //compute histogram of targets in neighborhood
    map<string,unsigned> histogram;
    for (unsigned i=0;i<aNeighborhood.size();++i){
      unsigned nn_id=aNeighborhood[i];
      assert(nn_id<aTargetList.size());      
      string predicted_target=aTargetList[nn_id];
      if (histogram.count(predicted_target)==0) histogram[predicted_target]=1;
      else histogram[predicted_target]++;
    }
    //compute majority vote for target
    string max_target=aTargetList[0];//initialization with one arbitrary target
    unsigned max_val=0;
    for (map<string,unsigned>::const_iterator it=histogram.begin();it!=histogram.end();++it){
      string target=it->first;
      unsigned vote=it->second;
      if (max_val<vote) {
	max_val=vote;
	max_target=target;
      }
    }
    return max_target;
  }

  void OutputDatasetSplit(string& aSparseASCIITestInputFileName, string& aTrainTargetInputFileName, unsigned aNumCenters, unsigned aSizeThreshold, double aImbalanceTolerance){
    cout<<"Compute train split."<<endl;////
    //read test instances
    vector<SVector> test_dataset;
    InputSparse(aSparseASCIITestInputFileName, "ascii", test_dataset);

    //read train targets
    vector<string> train_target_list;
    InputStringList(aTrainTargetInputFileName,train_target_list);
    
    //randomly shuffle test instances
    vector<unsigned> test_id_list(test_dataset.size());
    for (unsigned i=0;i<test_dataset.size();i++) test_id_list.push_back(i);
    random_shuffle(test_id_list.begin(),test_id_list.end());
    //iterate until the prespecified number of centers that satisfy the required considtions are found
    vector<vector<double> > valid_center_list;
    cout<<"Computing test centers."<<endl;
    ProgressBar progress_bar;
    for (unsigned i=0; i<test_id_list.size();++i){
      progress_bar.Count();
      //check if it is possible to extract a valid train subset from the selected test instance
      unsigned test_id=test_id_list[i];
      //sort all training set according to distance from selected test intance
      vector<unsigned> neighborhood=ComputeTrueNeighborhood(test_dataset[test_id],train_target_list.size());
      //compute the sorted sequence of targets
      vector<string> neighborhood_target;
      for (unsigned j=0;j<neighborhood.size();j++)
	neighborhood_target.push_back(train_target_list[neighborhood[j]]);
      //return the largest number of neighbour such that the target distribution is uniform
      double positive_class_count=0;
      double negative_class_count=0;
      unsigned ksize=0;
      double kvalue=0;
      double nbalance=0;
      for (unsigned j=0;j<neighborhood_target.size();j++){
	if (neighborhood_target[j]=="1") positive_class_count++;
	else if (neighborhood_target[j]=="-1") negative_class_count++;
	double balance=(negative_class_count>0)?positive_class_count/negative_class_count:0;
	if (balance>1/aImbalanceTolerance && balance<aImbalanceTolerance) {
	  ksize=j;
	  unsigned train_id=neighborhood[j];
	  kvalue=dot(mDataset[train_id],test_dataset[test_id]);
	  nbalance=balance;
	}
      }
      //if a valid neighborhood has been found then store it
      if (ksize>aSizeThreshold) {
	vector<double> item;
	item.push_back((double)test_id);
	item.push_back((double)ksize);
	item.push_back(kvalue);
	item.push_back(nbalance);
	valid_center_list.push_back(item);
      }
    }

    //iterate over all aNumCenters valid centers
    unsigned dataset_counter=0;
    for (unsigned i=0; i<valid_center_list.size() && dataset_counter<aNumCenters;++i){
      assert(valid_center_list[i].size()==4);
      unsigned center_id=(unsigned)valid_center_list[i][0];
      //unsigned ksize=(unsigned)valid_center_list[i][1];
      double kvalue=valid_center_list[i][2];
      double nbalance=valid_center_list[i][3];

      //iterate over test set and collect test instances that are within kvalue from the center 
      vector<unsigned> part_test_list;
      for (unsigned j=0;j<test_dataset.size();++j){
	double k=dot(test_dataset[j],test_dataset[center_id]);
	if (k>=kvalue) part_test_list.push_back(j);
      }
      //if condistions are met, then output train-test datsets
      if (part_test_list.size()>aSizeThreshold){
	//output test set
	string ofnamets=aSparseASCIITestInputFileName+".test_"+stream_cast<string>(dataset_counter);
	ofstream ofsts(ofnamets.c_str());
	for (unsigned t=0;t<part_test_list.size();++t) ofsts<<part_test_list[t]+1<<endl;//+1 as ids range from 1 and not from 0
	ofsts.close();

	//collect training instances that are in the balanced neighborhood
	vector<unsigned> part_train_list;
	for (unsigned j=0;j<mDataset.size();++j){
	  double k=dot(mDataset[j],test_dataset[center_id]);
	  if (k>=kvalue) part_train_list.push_back(j);
	}
	//output training set
	string ofnametr=aSparseASCIITestInputFileName+".train_"+stream_cast<string>(dataset_counter);
	ofstream ofstr(ofnametr.c_str());
	for (unsigned t=0;t<part_train_list.size();++t) {
	  
	  ofstr<<part_train_list[t]+1<<endl;//+1 as ids range from 1 and not from 0
	}
	ofstr.close();

	//logging info
	LOGF<<dataset_counter<<" center_id:"<<center_id+1<<" kval:"<<kvalue<<" #train:"<<part_train_list.size()<<" #test:"<<part_test_list.size()<<" balance:"<<nbalance<<endl;
      
	dataset_counter++;
      }
    }
    cout<<endl<<"Computed "<<dataset_counter<<" valid centers"<<endl;
  }

  
  void OutputKernel(ostream& out){
    cout<<"Compute kernel matrix."<<endl;////
    ProgressBar progress_bar;
    for (unsigned i=0;i<mDataset.size();i++){
      for (unsigned j=0;j<mDataset.size();j++)
	out<<Kernel(i,j)<<" ";
      out<<endl;
      progress_bar.Count();
    }
  }

  void Output(ostream& out){
    for (unsigned i=0;i<mDataset.size();i++)
      out<<mDataset[i];
  }

  void OutputFeatureMap(string aFileName)const{
    {
      string ofname=aFileName+".feature_map";
      ofstream of(ofname.c_str());
      pmFeatureGenerator->OutputFeatureMap(of);
    }
    {
      string ofname=aFileName+".abstract_feature_map";
      ofstream of(ofname.c_str());
      pmAbstractFeatureGenerator->OutputFeatureMap(of);
    }
    {
      string ofname=aFileName+".soft_feature_map";
      ofstream of(ofname.c_str());
      pmSoftFeatureGenerator->OutputFeatureMap(of);
    }
    {
      string ofname=aFileName+".gapped_feature_map";
      ofstream of(ofname.c_str());
      pmGappedFeatureGenerator->OutputFeatureMap(of);
    }
  }

  double Kernel(unsigned aI, unsigned aJ){
    if (mCaching){
      unsigned i=min(aI,aJ);
      unsigned j=max(aI,aJ);
      pair<unsigned,unsigned> key=make_pair(i,j); 
      if (mKernelMap.count(key)==0){
	double value=dot(mDataset[i],mDataset[j]);
	mKernelMap[key]=value;
      }
      return mKernelMap[key];
    } else 
      return dot(mDataset[aI],mDataset[aJ]);
  }
};

  //---------------------------------------------------------------------------------
  int main(int argc, char** argv){
    Timer T;
    srand(time(NULL));
    LOGF<<"--------------------------------------------------------------------------------"<<endl;
    LOGF<<CREDITS<<endl;
    time_t rawtime_start; time ( &rawtime_start ); LOGF<<"Start logging: "<<asctime ( localtime ( &rawtime_start ) )<<endl;
    LOGF<<"Command line: ";for(int i=0;i<argc;i++) LOGF<<stream_cast<string>(argv[i])<<" ";LOGF<<endl;
    try
      {
	ParameterWrapperClass P;
	P.Init(argc,argv);

	string mode="";
	//factory
	NSPDK_FeatureGenerator fg("nspdk");
	ANSPDK_FeatureGenerator afg("anspdk");
	ALNSPDK_FeatureGenerator alfg("alnspdk");
	MNSPDK_FeatureGenerator mfg("mnspdk");
	NSPDK3D_FeatureGenerator fg3d("nspdk3d");
	GNSPDK_FeatureGenerator gfg("gnspdk");

	NSPDK_FeatureGenerator* pfg;
	if (P.mType=="nspdk") pfg=&fg;
	else if (P.mType=="gnspdk") pfg=&gfg;
	else if (P.mType=="anspdk") pfg=&afg;
	else if (P.mType=="alnspdk") pfg=&alfg;
	else if (P.mType=="mnspdk") pfg=&mfg;
	else if (P.mType=="nspdk3d") {
	  pfg=&fg3d;
	  mode="3D";
	}
	else if (P.mType=="nspdkvp") {
	  pfg=&fg;
	  mode="viewpoint";
	}
	else throw range_error("Unknown feature generator type:"+P.mType);

	if (P.mOutputFeatureMap) P.mDebug+=1;//if the output of the feature map is required then the Debug level has to be at least 1 

	pfg->set_flag("radius",stream_cast<string>(P.mRadiusMax));
	pfg->set_flag("distance",stream_cast<string>(P.mDistanceMax));
	pfg->set_flag("match_type",stream_cast<string>(P.mMatchingType));
	pfg->set_flag("hash_bit_size",stream_cast<string>(P.mFeatureBitSize));
	pfg->set_flag("hash_bit_mask",stream_cast<string>((2 << P.mFeatureBitSize)-1));
	pfg->set_flag("verbosity",stream_cast<string>(P.mDebug));
	if (P.mMinKernel) pfg->set_flag("min_kernel","true");
	if (!P.mNormalization) pfg->set_flag("normalization","false");

	//feature generator for additional abstract mode
	const unsigned abstract_radius=2;
	const unsigned abstract_distance=100;
	
	NSPDK_FeatureGenerator abstract_fg("abstract_nspdk");
	NSPDK_FeatureGenerator* p_abstract_fg=&abstract_fg;
	p_abstract_fg->set_flag("radius",stream_cast<string>(abstract_radius));
	p_abstract_fg->set_flag("distance",stream_cast<string>(abstract_distance));
	p_abstract_fg->set_flag("match_type","hard");
	p_abstract_fg->set_flag("hash_bit_size",stream_cast<string>(P.mFeatureBitSize));
	p_abstract_fg->set_flag("hash_bit_mask",stream_cast<string>((2 << P.mFeatureBitSize)-1));
	p_abstract_fg->set_flag("verbosity",stream_cast<string>(P.mDebug));
	if (P.mMinKernel) p_abstract_fg->set_flag("min_kernel","true");
	if (!P.mNormalization) p_abstract_fg->set_flag("normalization","false");

	//feature generator for additional soft mode
	NSPDK_FeatureGenerator soft_fg("soft_nspdk");
	NSPDK_FeatureGenerator* p_soft_fg=&soft_fg;
	p_soft_fg->set_flag("radius",stream_cast<string>(P.mRadiusMax));
	p_soft_fg->set_flag("distance",stream_cast<string>(P.mDistanceMax));
	p_soft_fg->set_flag("match_type","soft");
	p_soft_fg->set_flag("hash_bit_size",stream_cast<string>(P.mFeatureBitSize));
	p_soft_fg->set_flag("hash_bit_mask",stream_cast<string>((2 << P.mFeatureBitSize)-1));
	p_soft_fg->set_flag("verbosity",stream_cast<string>(P.mDebug));
	if (P.mMinKernel) p_soft_fg->set_flag("min_kernel","true");
	if (!P.mNormalization) p_soft_fg->set_flag("normalization","false");

	//feature generator for additional gapped mode
	NSPDK_FeatureGenerator gapped_fg("gapped_nspdk");
	NSPDK_FeatureGenerator* p_gapped_fg=&gapped_fg;
	p_gapped_fg->set_flag("radius",stream_cast<string>(P.mRadiusMax));
	p_gapped_fg->set_flag("distance",stream_cast<string>(P.mDistanceMax));
	p_gapped_fg->set_flag("match_type","hard");
	p_gapped_fg->set_flag("hash_bit_size",stream_cast<string>(P.mFeatureBitSize));
	p_gapped_fg->set_flag("hash_bit_mask",stream_cast<string>((2 << P.mFeatureBitSize)-1));
	p_gapped_fg->set_flag("verbosity",stream_cast<string>(P.mDebug));
	if (P.mMinKernel) p_gapped_fg->set_flag("min_kernel","true");
	if (!P.mNormalization) p_gapped_fg->set_flag("normalization","false");


	string ofname;
	//main process
	MainProcessClass C(pfg,
			   p_abstract_fg,
			   p_soft_fg,
			   p_gapped_fg,
			   P.mAddAbstractMatchingType,
			   P.mAddSoftMatchingType,
			   P.mAddGappedMatchingType,
			   P.mAddStemInsensitiveMatchingType,
			   P.mAddPPInsensitiveMatchingType,
			   P.mNumHashFunctions,
			   P.mHashFactor,
			   P.mNumNearestNeighbors,
			   P.mSampleSize,
			   P.mNonRedundantFilter,
			   mode, 
			   P.mCaching, 
			   P.mTrueSort, 
			   P.mEccessNeighbourSizeFactor,
			   P.mOutputHashEncoding,
			   P.mWhiteListFileName,
			   P.mBlackListFileName);
	if (P.mOutputFeatures) {
	  C.InputAndOutput(P.mGspanInputFileName);
	  if (P.mOutputFeatureMap)
	    C.OutputFeatureMap(P.mGspanInputFileName);
		
	} else if (P.mOutputWeightAnnotation) {
	  string mode="read_absolute_vertex_id";
	  C.SetMode(mode);
	  if (P.mWeightFileName=="") throw range_error("ERROR:No weight file name specified");
	  C.LoadAndAnnotate(P.mGspanInputFileName,P.mWeightFileName);
	} else {
	  if (P.mSparseBinaryInputFileName!="")
	    C.InputSparse(P.mSparseBinaryInputFileName,"binary");
	  else if (P.mSparseASCIIInputFileName!="")
	    C.InputSparse(P.mSparseASCIIInputFileName,"ascii");
	  else if (P.mGspanInputFileName!="")
	    C.Input(P.mGspanInputFileName);
	  else
	    throw range_error("ERROR:No input file name specified");

	  if (P.mOutputCluster!="NONE"){
	    C.ComputeBinDataStructure(P.mNumHashFunctions);
	    C.OutputBinDataStructureStatistics();
	    ofname=P.mGspanInputFileName+P.mSparseASCIIInputFileName+P.mSparseBinaryInputFileName+".fast_cluster";
	    ofstream ofs_fc(ofname.c_str());
	    C.OutputCluster(ofs_fc,P.mOutputCluster,P.mMaxIntersectionSize);
	  }

	  if (P.mOutputApproximateKNN){
	    C.ComputeBinDataStructure(P.mNumHashFunctions);
	    ofname=P.mGspanInputFileName+P.mSparseASCIIInputFileName+P.mSparseBinaryInputFileName+".approx_knn";
	    ofstream ofs_aknn(ofname.c_str());
	    C.OutputApproximateKNN(ofs_aknn);
	  }

	  if (P.mOutputFastKNN){
	    C.ComputeInvertedFeatureInstanceIndex();
	    ofname=P.mGspanInputFileName+P.mSparseASCIIInputFileName+P.mSparseBinaryInputFileName+".fast_knn";
	    ofstream ofs_fknn(ofname.c_str());
	    C.OutputFastKNN(ofs_fknn);
	  }

	  if (P.mOutputTrueKNN){
	    ofname=P.mGspanInputFileName+P.mSparseASCIIInputFileName+P.mSparseBinaryInputFileName+".knn";
	    ofstream ofs_fknn(ofname.c_str());
	    C.OutputTrueKNN(ofs_fknn);
	  }

	  if (P.mOutputKernel){
	    ofname=P.mGspanInputFileName+P.mSparseASCIIInputFileName+P.mSparseBinaryInputFileName+".kernel";
	    ofstream ofs_fk(ofname.c_str());
	    C.OutputKernel(ofs_fk);
	  }

	  if (P.mOutputApproximateKNNPrediction){
	    C.ComputeBinDataStructure(P.mNumHashFunctions);
	    ofname=P.mGspanInputFileName+P.mSparseASCIIInputFileName+P.mSparseBinaryInputFileName+".approx_knn_prediction";
	    ofstream ofs_knnp(ofname.c_str());
	    C.OutputApproximateKNNPrediction(ofs_knnp,P.mSparseASCIITestInputFileName,P.mTrainTargetInputFileName);
	  }

	  if (P.mOutputTrueKNNPrediction){
	    ofname=P.mGspanInputFileName+P.mSparseASCIIInputFileName+P.mSparseBinaryInputFileName+".knn_prediction";
	    ofstream ofs_knnp(ofname.c_str());
	    C.OutputTrueKNNPrediction(ofs_knnp,P.mSparseASCIITestInputFileName,P.mTrainTargetInputFileName);
	  }
	  if (P.mOutputDatasetSplit){
	    C.OutputDatasetSplit(P.mSparseASCIITestInputFileName,P.mTrainTargetInputFileName,P.mNumCenters, P.mSizeThreshold,P.mImbalanceTolerance);
	  }
	}
      }
    catch(exception& e)
      {
	cerr<<e.what();
	LOGF<<e.what()<<endl;
      }
    time_t rawtime_end; time ( &rawtime_end );LOGF<<"End logging: "<<asctime ( localtime ( &rawtime_end ) )<<endl;
    LOGF<<"Time elapsed: "<<T.end()<<endl;
    return 0;
  }
