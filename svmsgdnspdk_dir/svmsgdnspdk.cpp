#include "Utility.h"
#include "Histogram.h"
#include "BaseGraphClass.h"
#include "GraphClass.h"
#include "NSPDK_FeatureGenerator.h"

#include "vectors.h"
#include "gzstream.h"
#include <iostream>
#include <iomanip>
#include <string>
#include <map>
#include <vector>
#include <cassert>
#include <cstdlib>
#include <cmath>

using namespace std;

const string PROG_CREDIT="SVMSGDNSPDK Vers. 7.0 (4 Jan 2012)\nAuthor: Fabrizio Costa costa@informatik.uni-freiburg.de";
const string SEP="--------------------------------------------------------------";
const string TAB="    ";

FlagsService& The_FlagsService= FlagsService::get_instance();

//------------------------------------------------------------------------------------------------------------------------
inline double random01(){
  return (double)rand()/(double)RAND_MAX;
}

inline unsigned randomUnsigned(unsigned aMax){
  return (unsigned)(rand() % aMax);
}
  
//------------------------------------------------------------------------------------------------------------------------
class ParameterWrapperClass{
public:
  string mAction;
  string mDataFileName;
  string mTargetFileName;
  string mModelFileName;
  string mMode;
  string mFileType;
  vector<pair<unsigned,unsigned> > mLevelList;
  bool mThicknessOpt;
  double mLambda;
  int mEpochs;
  double mDataSizeThreshold;
  double mInstancePartsThreshold;
  unsigned mBitSize;
  unsigned mCrossValidationNumFolds;
  unsigned mLearningCurveNumPoints;
  string mThresholdMode;
  bool mFeatureCache;
  unsigned mMaxIntersectionSize;
  unsigned mNumHashFunctions;
  unsigned mNumCompactHashFunctions;
  unsigned mNeighborhoodSize;
  bool mDataFilter;
  unsigned mRandSeed;
  string mKernelType;
  string mGraphType;
  unsigned mSemiSupervisedIteration;
  double mSemiSupervisedThreshold;
  bool mSemiSupervisedInduceOnlyPositive;
  bool mSemiSupervisedInduceOnlyNegative;
  string mPrefix;
  double mOneClassSizeRatio;
  bool mIntraInstanceRegularization;
  double mIntraInstanceRegularizationLambda;
public:
  ParameterWrapperClass(): 			   
    mAction("TRAIN"),
    mDataFileName(""),
    mTargetFileName(""),
    mModelFileName("model"),
    mMode("MEMORY"),
    mFileType("GSPAN"),
    mLevelList(vector<pair<unsigned,unsigned> >()),
    mThicknessOpt(false),
    mLambda(1e-4),
    mEpochs(10),
    mDataSizeThreshold(.5),
    mInstancePartsThreshold(1),
    mBitSize(14),
    mCrossValidationNumFolds(10),
    mLearningCurveNumPoints(10),
    mThresholdMode("PROPORTIONAL"),
    mFeatureCache(false),
    mMaxIntersectionSize(1),
    mNumHashFunctions(200),
    mNumCompactHashFunctions(1),
    mNeighborhoodSize(10),
    mDataFilter(false),
    mRandSeed(1),
    mKernelType("HARD"),
    mGraphType("UNDIRECTED"),
    mSemiSupervisedIteration(3),
    mSemiSupervisedThreshold(1),
    mSemiSupervisedInduceOnlyPositive(false),
    mSemiSupervisedInduceOnlyNegative(false),
    mPrefix(""),
    mOneClassSizeRatio(1),
    mIntraInstanceRegularization(false),
    mIntraInstanceRegularizationLambda(1e-5){}
				    				      
  void Usage(string aCommandName){
    cerr<<"Usage: "<<aCommandName<<endl
	<<"[-a <action TRAIN|TEST|TEST-PART|CROSSVALIDATION|LEARNINGCURVE> (default: "<< mAction <<")]"<<endl
	<<"-d <file name for data in (optionally) gzipped format>"<<endl
	<<"-t <file name for target (single column format of 1|-1)>"<<endl
	<<"-ll <num of levels followed by list of pairs (radius,distance)>"<<endl<<"    eg. 2 1 4 2 5 => 2 levels, the first one with radius=1 and distance=4 and the second with radius=2 and distance=5"<<endl
	<<"[-m <file name for model file> (default:"<<mModelFileName<<")]"<<endl
	<<"[-f <file type GSPAN|SDF> (default: "<< mFileType <<")]"<<endl
	<<"[-gt <graph type UNDIRECTED|DIRECTED> (default: "<< mGraphType <<")]"<<endl
	<<"[-kt <kernel type HARD|SOFT> (default: "<< mKernelType <<")]"<<endl
	<<"[-topt <thickness features option flag> (default: "<< mThicknessOpt <<")]"<<endl
	<<"[-iiropt <intra instance regularization option flag> (default: "<< mIntraInstanceRegularization <<")]"<<endl
	<<"[-iirl <intra instance regularization lambda> (default: "<< mIntraInstanceRegularizationLambda <<")]"<<endl

	<<SEP<<endl<<"    Training options"<<endl<<SEP<<endl
	<<"[-mode <mode MEMORY|FILE> (default: "<< mMode <<")]"<<endl<<"    in FILE mode data is read at each iteration from disk and retained in memory in sparse vector format;"<<endl<<"    in MEMORY mode data is read once from disk and is retained in memory in graph format."<<endl
	<<"[-l <lambda: real> (default: "<< mLambda <<")]"<<endl
	<<"[-e <epochs: integer> (default: "<< mEpochs <<")]"<<endl
	<<"[-b <bit size: integer> (default: "<< mBitSize <<")]"<<endl
	<<"[-cv <cross validation num of folds: integer> (default: "<< mCrossValidationNumFolds <<")]"<<endl
	<<"[-lc <learning curve num of evaluation points: integer> (default: "<< mLearningCurveNumPoints <<")]"<<endl<<"    The total number of instances is randomly divided in lc sets of same size; one set is used for test"<<endl
	<<"[-rs <rand seed: integer> (default: "<< (mRandSeed) <<")]"<<endl
	<<"[-pfx <prefix for output file names: string> (default: \""<< mPrefix <<"\")]"<<endl

	<<SEP<<endl<<"    Fast training options"<<endl<<SEP<<endl
	<<"[-tm <threshold mode PERFORMANCE|PROPORTIONAL> (default: "<< mThresholdMode <<")]"<<endl<<"    in PERFORMANCE mode, parameter dt is used to set a required balanced F-measure (computed on the training set). The classification of low confidence instances is be delegated to the successor level."<<endl<<"    in PROPORTIONAL mode, parameter dt identifies the fraction of instances that are classified at each level."<<endl
	<<"[-dt <data size threshold: double> (default: "<< mDataSizeThreshold <<")]"<<endl
	<<"[-it <instance parts threshold: double> (default: "<< mInstancePartsThreshold <<")]"<<endl<<"    indicates the fraction of vertices that are predicted to be interesting at a given level and that will be processed at successor level."<<endl
	<<"[-fc <feature cache flag> (default: "<< (mFeatureCache) <<")]"<<endl<<"    if set, features associations with vertices of all instances are memorized. In crossvalidation this ensures that each instance is processed only once."<<endl

	<<SEP<<endl<<"    Semi-supervised setting options"<<endl<<SEP<<endl
	<<"[-ssi <semi supervised number of iterations: integer> (default: "<< (mSemiSupervisedIteration) <<")]"<<endl
	<<"[-sst <semi supervised threshold: real> (default: "<< (mSemiSupervisedThreshold) <<")]"<<endl<<"    Only the top and low quantile will be used as positives and negative instances. A threshold of 1 means that all unsupervised instaces are used in the next phase."<<endl
	<<"[-ssp <induce only positive class instances in semi supervised setting flag> (default: "<< (mSemiSupervisedInduceOnlyPositive) <<")]"<<endl
	<<"[-ssn <induce only negative class instances in semi supervised setting flag> (default: "<< (mSemiSupervisedInduceOnlyNegative) <<")]"<<endl

	<<SEP<<endl<<"    Once-class setting options"<<endl<<SEP<<endl
	<<"[-ocsr <one-class size ratio for the set of negative class instances> (default: "<< (mOneClassSizeRatio) <<")]"<<endl

	<<SEP<<endl<<"    Data pre-filtering options"<<endl<<SEP<<endl
	<<"[-df <data filter flag> (default: "<< mDataFilter <<")]"<<endl<<"    if set, features extracted from the first level are used in the min-hash table to compute approximate neighborhoods. The average pairwise similarity in each neighborhood is used to rank instances. The top half of the most dense instances whose neighborhood does not overlap more than the max intersection size (as specifiedd by the flag -mi) is reteined as valid training instances."<<endl
	<<"[-nh <number of hash functions: integer> (default: "<< mNumHashFunctions <<")]"<<endl
	<<"[-nch <number of hash functions compacted into a signle signature: integer> (default: "<< mNumCompactHashFunctions <<")]"<<endl
	<<"[-ns <neighbourhood size: integer> (default: "<< mNeighborhoodSize <<")]"<<endl
	<<"[-mi <max intersection size: integer> (default: "<< mMaxIntersectionSize <<")]"<<endl;

    exit(0);
  }

  void Init(int argc, const char** argv){
    if (argc==1) Usage(argv[0]); 
    vector<string> options;
    for (int i=1;i<argc;i++) options.push_back(argv[i]);
    for (vector<string>::iterator it=options.begin();it!=options.end();++it)
      {
	if ((*it)=="-h" || (*it)=="--help") Usage(argv[0]); 
	else if ((*it)=="-a") mAction=(*(++it)); 
	else if ((*it)=="-d") mDataFileName=(*(++it)); 
	else if ((*it)=="-t") mTargetFileName=(*(++it)); 
	else if ((*it)=="-m") mModelFileName=(*(++it)); 
	else if ((*it)=="-f") mFileType=(*(++it)); 
	else if ((*it)=="-mode") mMode=(*(++it)); 
	else if ((*it)=="-l") mLambda=stream_cast<double>(*(++it)); 
	else if ((*it)=="-e") mEpochs=stream_cast<unsigned>(*(++it)); 
	else if ((*it)=="-dt") mDataSizeThreshold=stream_cast<double>(*(++it)); 
	else if ((*it)=="-it") mInstancePartsThreshold=stream_cast<double>(*(++it)); 
	else if ((*it)=="-b") mBitSize=stream_cast<unsigned>(*(++it)); 
	else if ((*it)=="-cv") mCrossValidationNumFolds=stream_cast<unsigned>(*(++it)); 
	else if ((*it)=="-lc") mLearningCurveNumPoints=1+stream_cast<unsigned>(*(++it)); //NOTE:add 1 fold to reserve material for testing
	else if ((*it)=="-tm") mThresholdMode=(*(++it)); 
	else if ((*it)=="-fc") mFeatureCache=true;
	else if ((*it)=="-mi") mMaxIntersectionSize=stream_cast<unsigned>(*(++it)); 
	else if ((*it)=="-nh") mNumHashFunctions=stream_cast<unsigned>(*(++it)); 
	else if ((*it)=="-nch") mNumCompactHashFunctions=stream_cast<unsigned>(*(++it)); 
	else if ((*it)=="-ns") mNeighborhoodSize=stream_cast<unsigned>(*(++it)); 
	else if ((*it)=="-rs") mRandSeed=stream_cast<unsigned>(*(++it)); 
	else if ((*it)=="-ssi") mSemiSupervisedIteration=stream_cast<unsigned>(*(++it)); 
	else if ((*it)=="-sst") mSemiSupervisedThreshold=stream_cast<double>(*(++it)); 
	else if ((*it)=="-ssp") mSemiSupervisedInduceOnlyPositive=true;
	else if ((*it)=="-ssn") mSemiSupervisedInduceOnlyNegative=true;
	else if ((*it)=="-ocsr") mOneClassSizeRatio=stream_cast<double>(*(++it)); 
	else if ((*it)=="-df") mDataFilter=true;
	else if ((*it)=="-kt") mKernelType=(*(++it)); 
	else if ((*it)=="-gt") mGraphType=(*(++it));
	else if ((*it)=="-topt") mThicknessOpt=true; 
	else if ((*it)=="-pfx") mPrefix=(*(++it)); 
	else if ((*it)=="-iiropt") mIntraInstanceRegularization=true; 
	else if ((*it)=="-iirl") mIntraInstanceRegularizationLambda=stream_cast<double>(*(++it)); 
	else if ((*it)=="-ll") {
	  unsigned num_levels=stream_cast<unsigned>(*(++it)); 
	  for (unsigned i=0;i<num_levels;i++){
	    unsigned radius=stream_cast<unsigned>(*(++it)); 
	    unsigned distance=stream_cast<unsigned>(*(++it)); 
	    mLevelList.push_back(make_pair(radius,distance));
	  }
	}
	else {cerr<<"Unrecognized parameter: "<<(*it)<<"."<<endl;throw exception();}
      }
    if (mDataFileName=="")
      throw range_error("ERROR: -d <file> is missing.");
    if ((mAction=="TRAIN" || mAction=="CROSSVALIDATION") && mTargetFileName=="")
      throw range_error("ERROR: -t <file> is missing.");
    if ((mAction=="TEST" || mAction=="TEST-PART") && mModelFileName=="")
      throw range_error("ERROR: -m <file> is missing.");
    if ((mAction=="TRAIN" || mAction=="CROSSVALIDATION" )&& mLevelList.size()==0)       
      throw range_error("ERROR: Missing -ll parameter.");
  }
} PARAM;


//---------------------------------------------------------------------------------
inline int IntHashSimple(int key, int aModulo=RAND_MAX){
  key = ~key + (key << 15); // key = (key << 15) - key - 1;
  key = key ^ (key >> 12);
  key = key + (key << 2);
  key = key ^ (key >> 4);
  key = key * 2057; // key = (key + (key << 3)) + (key << 11);
  key = key ^ (key >> 16);
  return key%aModulo;
}

inline int IntHash(int key, int aModulo=RAND_MAX, unsigned aSeed=0){
  const double A=sqrt(2)-1;
  return IntHashSimple(key*(aSeed+1)*A,aModulo);
}

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


//------------------------------------------------------------------------------------------------------------------------
class TimerClass {
public:
  typedef double diff_type;
  TimerClass(): start_sec(time(NULL)),start(std::clock()){}
  ~TimerClass(){
    //std::clock_t end=std::clock();
    //std::clock_t elapsed=end-start;
    std::time_t end_sec=time(NULL);
    diff_type elapsed_sec=end_sec-start_sec;
    diff_type elapsed_min=elapsed_sec/60;
    diff_type elapsed_hour=elapsed_min/60;
    //CLOCKS_PER_SEC
    cout<<"Elapsed time: h: "<<floor(elapsed_hour)<<" m: "<<floor(elapsed_min)<<" s: "<<floor(elapsed_sec)<<endl;
  }
private:
  std::time_t start_sec;
  std::clock_t start;
};

//------------------------------------------------------------------------------------------------------------------------
class ProgressBar{
public:
  ProgressBar(unsigned aStep=100):mStep(aStep),mCounter(0){}
  ~ProgressBar(){if (mCounter>1) cout<<endl<<"Counted "<<mCounter<<" times."<<endl;}
  void Begin(){mCounter=0;}
  void Count(){
    mCounter++;
    if (mCounter%mStep==0) cout<<"."<<flush;
    if (mCounter%1000==0) cout<<mCounter/(1000)<<"K"<<flush;
  }
  unsigned End(){return mCounter;}  
private:
  unsigned mStep;
  unsigned mCounter;
  TimerClass mTimer;
};

//------------------------------------------------------------------------------------------------------------------------


  // Available losses
#define HINGELOSS 1
#define SMOOTHHINGELOSS 2
#define SQUAREDHINGELOSS 3
#define LOGLOSS 10
#define LOGLOSSMARGIN 11

// Select loss: NOTE: the selection done in the makefile
//#define LOSS LOGLOSS
//#define LOSS HINGELOSS
//#define LOSS SMOOTHHINGELOSS

  // Zero when no bias
  // One when bias term
#define BIAS 1


  inline 
  double loss(double z)
  {
#if LOSS == LOGLOSS
    if (z > 18)
      return exp(-z);
    if (z < -18)
      return -z;
    return log(1+exp(-z));
#elif LOSS == LOGLOSSMARGIN
    if (z > 18)
      return exp(1-z);
    if (z < -18)
      return 1-z;
    return log(1+exp(1-z));
#elif LOSS == SMOOTHHINGELOSS
    if (z < 0)
      return 0.5 - z;
    if (z < 1)
      return 0.5 * (1-z) * (1-z);
    return 0;
#elif LOSS == SQUAREDHINGELOSS
    if (z < 1)
      return 0.5 * (1 - z) * (1 - z);
    return 0;
#elif LOSS == HINGELOSS
    if (z < 1)
      return 1 - z;
    return 0;
#else
# error "Undefined loss"
#endif
  }

  inline 
  double dloss(double z)
  {
#if LOSS == LOGLOSS
    if (z > 18)
      return exp(-z);
    if (z < -18)
      return 1;
    return 1 / (exp(z) + 1);
#elif LOSS == LOGLOSSMARGIN
    if (z > 18)
      return exp(1-z);
    if (z < -18)
      return 1;
    return 1 / (exp(z-1) + 1);
#elif LOSS == SMOOTHHINGELOSS
    if (z < 0)
      return 1;
    if (z < 1)
      return 1-z;
    return 0;
#elif LOSS == SQUAREDHINGELOSS
    if (z < 1)
      return (1 - z);
    return 0;
#else
    if (z < 1)
      return 1;
    return 0;
#endif
  }

//------------------------------------------------------------------------------------------------------------------------
class DataAccessorClass{
protected:
  string mFileType;
  int mOffset;
public:
  DataAccessorClass(string aFileType):mFileType(aFileType){   
    if (mFileType=="GSPAN") mOffset=1;
    else if (mFileType=="SDF") mOffset=0;
    else throw range_error("ERROR: Unricognized file format:"+mFileType);    
  }
  virtual void GetGraph(int aIndex, GraphClass& oG)=0;

  virtual void Reset()=0;

  virtual unsigned Size()=0;

  void Skip(istream& in){
    if (mFileType=="GSPAN") SkipGSPAN(in);
    else if (mFileType=="SDF") SkipSDF(in);
    else throw range_error("ERROR: Unricognized file format:"+mFileType);
  }
  void SetGraphFromFile(istream& in, GraphClass& oG){
    if (mFileType=="GSPAN") SetGraphFromFileGSPAN(in,oG);
    else if (mFileType=="SDF") SetGraphFromFileSDF(in,oG);
    else throw range_error("ERROR: Unricognized file format:"+mFileType);    
  }
  void SkipGSPAN(istream& in){
    do {
      string line;
      getline(in,line);
      if (line[0]=='t') break;
    } while(!in.eof());
  }

  void SetGraphFromFileGSPAN(istream& in, GraphClass& oG){
    //status
    vector<bool> status; 	 
    status.push_back(true);//kernel point
    status.push_back(true);//kind
    status.push_back(true);//viewpoint
    status.push_back(false);//dead
    status.push_back(false);//forbidden distance
    status.push_back(false);//forbidden neighborhood


    map<string,int> index_map_nominal_to_real;
    string line;
    unsigned line_counter=0;
    do{
      line_counter++;
      getline(in,line);
      if (line=="") break;
      stringstream ss;
      ss<<line<<endl;
      char code;
      ss>>code;
      if (code=='t') break;
      else if (code=='v'){
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
	oG.SetVertexSymbolicAttributeList(real_vertex_index, vertex_symbolic_attribute_list);
	oG.SetVertexStatusAttributeList(real_vertex_index, status);
      } else if (code=='e'){
	//extract src and dest vertex id 
	string nominal_src_index,nominal_dest_index;
	string label;
	ss>>nominal_src_index>>nominal_dest_index>>label;
	assert(index_map_nominal_to_real.count(nominal_src_index)>0);
	assert(index_map_nominal_to_real.count(nominal_dest_index)>0);
	vector<string> edge_symbolic_attribute_list;
	edge_symbolic_attribute_list.push_back(label);
	unsigned real_src_index=index_map_nominal_to_real[nominal_src_index];
	unsigned real_dest_index=index_map_nominal_to_real[nominal_dest_index];
	unsigned edge_index=oG.InsertEdge(real_src_index,real_dest_index);
	oG.SetEdgeSymbolicAttributeList(edge_index,edge_symbolic_attribute_list);

	if (PARAM.mGraphType=="UNDIRECTED"){
	  unsigned reverse_edge_index=oG.InsertEdge(real_dest_index,real_src_index);
	  oG.SetEdgeSymbolicAttributeList(reverse_edge_index,edge_symbolic_attribute_list);
	}
      } 
      else {}//NOTE: ignore other markers
    } while(!in.eof() && in.good());
    if (PARAM.mGraphType=="DIRECTED"){
      unsigned vsize=oG.VertexSize();
      //add a copy of all vertices
      for (unsigned i=0;i<vsize;i++){
	unsigned real_vertex_index=oG.InsertVertex();
	assert(real_vertex_index=i+vsize);
	vector<string> vertex_symbolic_attribute_list=oG.GetVertexSymbolicAttributeList(i);
	for (unsigned t=0;t<vertex_symbolic_attribute_list.size();t++)//prepend a prefix to mark the reverse direction
	  vertex_symbolic_attribute_list[t]="r."+vertex_symbolic_attribute_list[t];
	oG.SetVertexSymbolicAttributeList(real_vertex_index,vertex_symbolic_attribute_list);
	oG.SetVertexStatusAttributeList(real_vertex_index, status);
      }
      //copy all edges swapping src with dest
      for (unsigned i=0;i<vsize;i++){
	//get all edges
	vector<unsigned> adj=oG.GetVertexAdjacentList(i);
	for (unsigned j=0;j<adj.size();j++){
	  unsigned orig_src=i;
	  unsigned orig_dest=adj[j];
	  unsigned reverse_src=orig_dest+vsize;
	  unsigned reverse_dest=orig_src+vsize;
	  unsigned edge_index=oG.InsertEdge(reverse_src,reverse_dest);
	  oG.SetEdgeSymbolicAttributeList(edge_index,oG.GetEdgeSymbolicAttributeList(orig_src,orig_dest));
	}
      }
    }
  }

  void SkipSDF(istream& in){
    do {
      string line;
      getline(in,line);
      if (line=="$$$$") break;
    } while(!in.eof());
  }

  void SetGraphFromFileSDF(istream& in, GraphClass& oG){
    //status
    vector<bool> status; 	 
    status.push_back(true);//kernel point
    status.push_back(true);//kind
    status.push_back(true);//viewpoint
    status.push_back(false);//dead
    status.push_back(false);//forbidden distance
    status.push_back(false);//forbidden neighborhood


    map<string,int> index_map_nominal_to_real;
    string line;
    unsigned line_counter=0;
    string num_vertex_str;
    unsigned num_vertex=0;
    unsigned num_edge=0;
    string num_edge_str;
    unsigned edge_counter=0;
    unsigned vertex_id=0;
    string vertex_label_str;
    string edge_src_id_str;
    string edge_dest_id_str;
    string edge_label_str;
    vector<string> vertex_symbolic_attribute_list;
    vertex_symbolic_attribute_list.push_back("");
    vector<string> edge_symbolic_attribute_list;
    edge_symbolic_attribute_list.push_back("");
    do{
      line_counter++;
      getline(in,line);
      if (line=="$$$$") break;
      if (line_counter==4) {//at 4th line expect num vertices and num edges 
	num_vertex_str=line.substr(0,3);
	num_vertex=stream_cast<unsigned>(num_vertex_str);
	num_edge_str=line.substr(3,3);
	num_edge=stream_cast<unsigned>(num_edge_str);
      }
      if (line_counter>4 && line_counter<=4+num_vertex){ //read vertex info
	vertex_label_str=line.substr(31,3);
	vertex_id++;
	string nominal_vertex_index=stream_cast<string>(vertex_id);
	unsigned real_vertex_index=oG.InsertVertex();
	index_map_nominal_to_real[nominal_vertex_index]=real_vertex_index;
	vertex_symbolic_attribute_list[0]=vertex_label_str;
	oG.SetVertexSymbolicAttributeList(real_vertex_index, vertex_symbolic_attribute_list);
	oG.SetVertexStatusAttributeList(real_vertex_index, status);
      }
      if (line_counter>4+num_vertex && line_counter<=4+num_vertex+num_edge){ //read edge info
	edge_counter++;
	edge_src_id_str=line.substr(0,3);
	edge_dest_id_str=line.substr(3,3);
	edge_label_str=line.substr(6,3);
	string nominal_src_index,nominal_dest_index;
	nominal_src_index=stream_cast<string>(stream_cast<unsigned>(edge_src_id_str));
	nominal_dest_index=stream_cast<string>(stream_cast<unsigned>(edge_dest_id_str));
	assert(index_map_nominal_to_real.count(nominal_src_index)>0);
	assert(index_map_nominal_to_real.count(nominal_dest_index)>0);
	unsigned real_src_index=index_map_nominal_to_real[nominal_src_index];
	unsigned real_dest_index=index_map_nominal_to_real[nominal_dest_index];
	unsigned edge_index=oG.InsertEdge(real_src_index,real_dest_index);
	unsigned reverse_edge_index=oG.InsertEdge(real_dest_index,real_src_index);
	edge_symbolic_attribute_list[0]=edge_label_str;
	oG.SetEdgeSymbolicAttributeList(edge_index,edge_symbolic_attribute_list);
	oG.SetEdgeSymbolicAttributeList(reverse_edge_index,edge_symbolic_attribute_list);
      }
    } while(!in.eof() && in.good());
    if (line_counter>1){
      //check post-conditions
      assert(edge_counter==num_edge && vertex_id==num_vertex);
    } else {}//if only one line has been read and we have arrived at eof then just end
  } 
};

//------------------------------------------------------------------------------------------------------------------------
class ReferenceDataAccessorClass: public DataAccessorClass{
protected:
public:  
  vector<GraphClass> *mpGraphList;
public:
  ReferenceDataAccessorClass(vector<GraphClass>& aGraphList,string aFileType=""):DataAccessorClass(aFileType){
    mpGraphList=&aGraphList;
  }
  void GetGraph(int aIndex, GraphClass& oG){
    if ((unsigned)aIndex>=mpGraphList->size()) 
      throw range_error("ERROR: Cannot access element with index:"+stream_cast<string>(aIndex)+" in container of size:"+stream_cast<string>(mpGraphList->size()));
    oG=(*mpGraphList)[aIndex];
  }
  void Reset(){}
  unsigned Size(){return mpGraphList->size();}
};

//------------------------------------------------------------------------------------------------------------------------
class MemoryDataAccessorClass: public DataAccessorClass{
protected:
public:  
  vector<GraphClass> mGraphList;
public:
  MemoryDataAccessorClass(string aFileName,string aFileType=""):DataAccessorClass(aFileType){
    igzstream fin;
    fin.open(aFileName.c_str());
    if (!fin) throw range_error("ERROR: Cannot open file: "+aFileName);
    ProgressBar pb;
    cout<<"Loading file:"<<aFileName<<endl;
    while (!fin.eof()){
      GraphClass g;
      SetGraphFromFile(fin,g);
      if (!g.IsEmpty()) {
	mGraphList.push_back(g); 
	pb.Count();
      }
    }
    
  }
  void GetGraph(int aIndex, GraphClass& oG){
    if ((unsigned)aIndex>=mGraphList.size()) throw range_error("ERROR: Cannot access element with index:"+stream_cast<string>(aIndex)+" in container of size:"+stream_cast<string>(mGraphList.size()));
    oG=mGraphList[aIndex];
  }
  void Reset(){}
  unsigned Size(){return mGraphList.size();}
};

class FileDataAccessorClass: public DataAccessorClass{
protected:
  igzstream* mFin;
  string mFileName;
  int mCounter;
public:
  FileDataAccessorClass(string aFileName,string aFileType=""):DataAccessorClass(aFileType),mFin(0),mFileName(aFileName),mCounter(-1){  }

  ~FileDataAccessorClass(){mFin->close();}

  void AlignIndexToFilePointer(int aIndex){//NOTE: here there are 2 types: either the marker is in the header (gspan) or in the footer (sdf)
    if (mCounter>aIndex) throw range_error("ERROR: Cannot access preceding element. Current counter:"+stream_cast<string>(mCounter)+" requested index:"+stream_cast<string>(aIndex));
    while (!mFin->eof() && mCounter<aIndex){
      Skip(*mFin);
      mCounter++;
    }
    if (mFin->eof() && mCounter<aIndex) throw range_error("ERROR: Reached file end but still could not access element of index:"+stream_cast<string>(aIndex));
  }


  void GetGraph(int aIndex, GraphClass& oG){
    AlignIndexToFilePointer(aIndex);
    SetGraphFromFile(*mFin,oG);
    mCounter++;
  }

  void Reset(){
    string line;
    mCounter=0;
    if (mFin!=0) mFin->close();
    mFin=new igzstream();
    mFin->open(mFileName.c_str());
    if (mOffset==1) Skip(*mFin); //NOTE: if mOffset==1 skip one header 
  }

  unsigned Size(){
    Reset();
    unsigned counter=0;
    while (!mFin->eof()) {
      Skip(*mFin);
      if (!mFin->eof()) counter++;
    }
    return counter+1;//add 1 for last element
  }
};

//------------------------------------------------------------------------------------------------------------------------
/**
   Prints several informative measures given a list of predictions and a list of true targets 
*/
void OutputPerformanceMeasures(ostream& out, const vector<double>& aPredictionList, const vector<double>& aTargetList){
  assert(aPredictionList.size()==aTargetList.size());
  unsigned size=aPredictionList.size();
  unsigned error=0;
  unsigned correct=0;
  unsigned abstain=0;
  unsigned tp,tn,fp,fn;tp=tn=fp=fn=0;
  for (unsigned i=0;i<aPredictionList.size();++i){
    double prediction=aPredictionList[i];
    double target=aTargetList[i];
    if (prediction==0) abstain++;
    else {
      if (prediction!=target) error++; 
      if (prediction==target) correct++;
      if (prediction>0 && target>0) tp++;
      if (prediction>0 && target<0) fp++;
      if (prediction<0 && target>0) fn++;
      if (prediction<0 && target<0) tn++;
    }
  }	

  double pprecision=(double)tp/(tp+fp);
  double precall=(double)tp/(tp+fn);
  double pfmeasure=2*pprecision*precall/(pprecision+precall);

  double nprecision=(double)tn/(tn+fn);
  double nrecall=(double)tn/(tn+fp);
  double nfmeasure=2*nprecision*nrecall/(nprecision+nrecall);

  double balanced_fmeasure=(pfmeasure+nfmeasure)/2;

  out<<TAB<<"Size: "<<size<<endl;
  out<<TAB<<"Predicted: "<<(size-abstain)<<" ( "<<(size-abstain)*100/(double)size<<" %) Abstained: "<<abstain<<" ( "<<abstain*100/(double)size<<" %)"<<endl;
  out<<TAB<<"Correct: "<< correct <<" ( "<<correct*100/(double)(correct+error)<<" %)"<<endl;
  out<<TAB<<"Error: "<<error<<" ( "<<error*100/(double)(correct+error)<<" %)"<<endl;
  out<<TAB<<"Confusion table:"<<endl;
  out<<TAB<<"TP:"<<tp<<" FP:"<<fp<<endl;
  out<<TAB<<"FN:"<<fn<<" TN:"<<tn<<endl;
  out<<TAB<<"+Precision:"<<pprecision<<" +Recall:"<<precall<<" +F-measure:"<<pfmeasure<<endl;
  out<<TAB<<"-Precision:"<<nprecision<<" -Recall:"<<nrecall<<" -F-measure:"<<nfmeasure<<endl;
  out<<TAB<<"Balanced F-measure:"<<balanced_fmeasure<<endl;
}



//------------------------------------------------------------------------------------------------------------------------
  /**
     Data structure to: 1) facilitate rebalancing of dataset by copying
     multiple times a reference to the instance; and 2) retrieve
     prediction efficiently by overwriting a reference to the margin
     list cell element.
  */
  struct TrainItem{
    int mInstanceID;
    double mTarget;
    double* mpMargin;
    SVector* mpInstance;
  };

  /**
     Encapsulates a linear SVM model trainable with stochastic gradient
     descent over graph instances explicitely mapped by the NSPDK kernel
  */
class ModelClass{
public:
  NSPDK_FeatureGenerator* mpFeatureGenerator;
protected:
  double  mLambda;
  unsigned mEpochs;
  unsigned mBitSize;
  double mInstancePartsThreshold;
  unsigned mRadius;
  unsigned mDistance;
  bool mAbstainingTypeFlag;
  double mPositiveThreshold;
  double mNegativeThreshold;
  SVector mW;
  double  mWScale;
  double  mBias;
  map<unsigned,map<unsigned,SVector> > mGraphVertexFeatureList;
  vector<multimap<unsigned,unsigned> > mBinDataStructure;
public:
  ModelClass(NSPDK_FeatureGenerator* apFeatureGenerator):mpFeatureGenerator(apFeatureGenerator){}
  /**
     Constructor that takes in input a lambda step, number of epoch, radius and distance for the NSPDK and a bool flag if to abstain
  */
  ModelClass(NSPDK_FeatureGenerator* apFeatureGenerator, double aLambda, unsigned aEpochs, unsigned aBitSize, double aInstancePartsThreshold, unsigned aRadius, unsigned aDistance, bool aAbstainingTypeFlag ):
    mpFeatureGenerator(apFeatureGenerator),
    mLambda(aLambda),
    mEpochs(aEpochs),
    mBitSize(aBitSize),
    mInstancePartsThreshold(aInstancePartsThreshold),
    mRadius(aRadius),
    mDistance(aDistance),
    mAbstainingTypeFlag(aAbstainingTypeFlag)
  {
    //Setup variables
    mPositiveThreshold=0;
    mNegativeThreshold=0;
    mWScale=1;
    mBias=0;
  }

  void SetupFeatureGenerator(){
    mpFeatureGenerator->set_flag("radius",stream_cast<string>(mRadius));
    mpFeatureGenerator->set_flag("distance",stream_cast<string>(mDistance));
    if (PARAM.mThicknessOpt) {
      mpFeatureGenerator->set_flag("thickness_opt",stream_cast<string>("true"));
      mpFeatureGenerator->set_flag("thickness",stream_cast<string>(0));
      mpFeatureGenerator->set_flag("thickness_distance",stream_cast<string>(mDistance));
    }
    unsigned bitmask = (2 << mBitSize)-1;
    mpFeatureGenerator->set_flag("hash_bit_mask",stream_cast<string>(bitmask));
  }


  void Save(ostream& out){
    out<<"bitsize "<<mBitSize<<endl;
    out<<"instance_parts_threshold "<<mInstancePartsThreshold<<endl;
    out<<"radius "<<mRadius<<endl;
    out<<"distance "<<mDistance<<endl;
    out<<"positive_threshold "<<mPositiveThreshold<<endl;
    out<<"negative_threshold "<<mNegativeThreshold<<endl;
    out<<"bias "<<mBias<<endl;
    out<<"wscale "<<mWScale<<endl;
    out<<"w "<<mW;
  }

  void Load(istream& in){
    string attribute="";
    string expected="";
    in>>attribute>>mBitSize; expected="bitsize";if (attribute!=expected) throw range_error("Format error: expecting ["+expected+"] but found ["+attribute+"]");
    in>>attribute>>mInstancePartsThreshold; expected="instance_parts_threshold";if (attribute!=expected) throw range_error("Format error: expecting ["+expected+"] but found ["+attribute+"]");
    in>>attribute>>mRadius; expected="radius";if (attribute!=expected) throw range_error("Format error: expecting ["+expected+"] but found ["+attribute+"]");
    in>>attribute>>mDistance; expected="distance";if (attribute!=expected) throw range_error("Format error: expecting ["+expected+"] but found ["+attribute+"]");
    in>>attribute>>mPositiveThreshold; expected="positive_threshold";if (attribute!=expected) throw range_error("Format error: expecting ["+expected+"] but found ["+attribute+"]");
    in>>attribute>>mNegativeThreshold; expected="negative_threshold";if (attribute!=expected) throw range_error("Format error: expecting ["+expected+"] but found ["+attribute+"]");
    in>>attribute>>mBias; expected="bias";if (attribute!=expected) throw range_error("Format error: expecting ["+expected+"] but found ["+attribute+"]");
    in>>attribute>>mWScale; expected="wscale";if (attribute!=expected) throw range_error("Format error: expecting ["+expected+"] but found ["+attribute+"]");
    in>>attribute>>mW; assert(attribute=="w");
  }

  SVector GenerateRandomFeaturePermutedInstance (vector<SVector>& aSVDataList){
    //assign iteratively features 
    unsigned i=0;
    //start from an empty instance
    SVector x;
    //select the size of the feature vector from an intance at random
    unsigned s=randomUnsigned(aSVDataList.size());
    unsigned size=(unsigned)aSVDataList[s].sparse_size();
    unsigned safeguard_counter=0;
    do {
      //select one instance uniformly at random
      unsigned j=randomUnsigned(aSVDataList.size());
      SVector& z=aSVDataList[j];
      //if the instance has fewer features than i, then skip
      if ((unsigned)z.sparse_size()<=i){safeguard_counter++;}
      else {
	//select the i-th <feature,value> pair
	pair<unsigned, double> feature_value=z.extract_component(i);
	//only if the feature has not yet been assigned to instance we are currently building
	if (safeguard_counter>1000 || x.get((int)feature_value.first)==0){//NOTE: use safeguard counter to manage situations where it becomes extremely unlikely to add a new feature
	  //...assign it 
	  x.set((int)feature_value.first,feature_value.second);
	  i++;
	  safeguard_counter=0;
	} else {safeguard_counter++;}
      }
    } while (i<size);
    return x;
  }

  void Train(vector<double>& aTargetList, vector<unsigned>& aTrainsetIDList,vector<vector<unsigned> >& aViewPointList, DataAccessorClass* apDataset, vector<double>& oMarginList, vector<double>& oPredictionList){
    if (aTrainsetIDList.size()!=aTargetList.size()) 
      throw range_error("ERROR: Data list and Target list have not the same size: #data:"+stream_cast<string>(aTrainsetIDList.size())+" #targets:"+stream_cast<string>(aTargetList.size()));
    if (aTrainsetIDList.size()!=aViewPointList.size()) 
      throw range_error("ERROR: Data list and viewpoint list have not the same size: #data:"+stream_cast<string>(aTrainsetIDList.size())+" #targets:"+stream_cast<string>(aViewPointList.size()));
    //Create vector representation from graphs 
    vector<SVector> sv_data_list;
    GraphToSparseVector(aTrainsetIDList,aViewPointList,apDataset,sv_data_list);
    vector<SVector> sv_data_list_orig;

    //allocate local target list and compute positive/negative target counts
    vector<double> target_list;
    unsigned p,n; p=n=0;
    for (unsigned i=0;i<aTargetList.size();++i) {
      if (aTargetList[i]>0) p++; else n++;
      target_list.push_back(aTargetList[i]);
    }

    if (n==0){//if no instance has negative class then generate mOneClassSizeRatio*aTrainsetIDList.size() negative instances with randomly permuted features 
      sv_data_list_orig=sv_data_list;//copy original dataset
      unsigned neg_size=PARAM.mOneClassSizeRatio*aTrainsetIDList.size();
      cout<<"No negative instances: proceeding to generate "<<neg_size<<" negative instances by composing randomly selected (feature, values) pairs extracted from randomly selected positive instances"<<endl;
#ifdef DEBUGON
      VectorClass size_vec;
#endif
      ProgressBar pb;
      for (unsigned i=0;i<neg_size;++i){
	SVector x=GenerateRandomFeaturePermutedInstance(sv_data_list);
#ifdef DEBUGON
	size_vec.PushBack(x.sparse_size());
#endif
	sv_data_list.push_back(x);
	target_list.push_back(-1);
	pb.Count();
      }
#ifdef DEBUGON
      cout<<endl<<"Feature size statistics for generated negatives:";
      size_vec.OutputStatistics(cout);
      cout<<endl;
#endif
    }

    //clear margin list and allocate memory
    oMarginList.clear();
    for (unsigned i=0;i<target_list.size();++i) oMarginList.push_back(0);
    //...rebalance classes use pointers array to scramble and oversmple
    vector<TrainItem> balanced_dataset;
    BalanceDataset(aTrainsetIDList,target_list,oMarginList,sv_data_list, balanced_dataset);

    //train on balanced train data
    TrainModule(balanced_dataset);

    if (mAbstainingTypeFlag){   
      //test on balanced train data
      vector<double> abst_margin_list;
      for (unsigned i=0;i<balanced_dataset.size();i++) abst_margin_list.push_back((*(balanced_dataset[i].mpMargin)));
      vector<double> abst_target_list;
      for (unsigned i=0;i<balanced_dataset.size();i++) abst_target_list.push_back(balanced_dataset[i].mTarget);
      //induce thresholds
      ThresholdInduction(abst_margin_list,abst_target_list);	
    }

    //output statistics on original train data (no class balance)
    OutputInfo();
    cout<<"Performance on train set:"<<endl;
    oPredictionList=ConvertMarginToPrediction(oMarginList);
    OutputPerformanceMeasures(cout,oPredictionList,target_list);
#ifdef DEBUGON
    SaveTrainPredictions(aTrainsetIDList,oMarginList,oPredictionList,target_list);
#endif

    if (n==0){//if no instance has negative class then re-test model only on real training data to extract margins and predictions
      oMarginList.clear();
      oPredictionList.clear();
      oMarginList=Test(sv_data_list_orig);
      oPredictionList=ConvertMarginToPrediction(oMarginList);
    }
  }
      
  void BalanceDataset(vector<unsigned>& aDatasetIDList,vector<double>& aTargetList,vector<double>& oMarginList, vector<SVector>& aSVDataList, vector<TrainItem>& oDataset){
    //compute class distribution
    unsigned p,n; p=n=0;
    for (unsigned i=0;i<aTargetList.size();++i) if (aTargetList[i]>0) p++; else n++;
    cout<<"Class distribution: "<< p+n <<" (+:"<<p<<" -:"<<n<<") "<<"[+:"<<(double)p/(p+n)<<" -:"<<(double)n/(p+n)<<"]"<<endl;

    //separate positive from negative instances
    vector<TrainItem> positive_data_list;
    vector<TrainItem> negative_data_list;
    if(aTargetList.size()!=aSVDataList.size()) throw range_error("ERROR: number of target values: "+stream_cast<string>(aTargetList.size())+" is different from dataset size:"+stream_cast<string>(aSVDataList.size()));
    for (unsigned i=0;i<aTargetList.size();++i) {
      TrainItem ti;
      ti.mTarget=aTargetList[i];
      ti.mpInstance=&aSVDataList[i];
      ti.mpMargin=&oMarginList[i];
      if (i<aDatasetIDList.size()){//syntesized instaces are appended after the real instances, so the size information of the original id_list marks the start of the syntesized instances
	ti.mInstanceID=aDatasetIDList[i];
      } else ti.mInstanceID=-1;//if the instance has been syntesized then it has no correspndent original graph 
      if (aTargetList[i]==1) positive_data_list.push_back(ti);
      else if (aTargetList[i]==-1) negative_data_list.push_back(ti);
      else throw range_error("ERROR: target has to be 1 or -1; cannot be: "+stream_cast<string>(aTargetList[i]));
    }
    //randomly shuffle data
    for (unsigned i=0;i<positive_data_list.size();++i){
      unsigned j=rand()*positive_data_list.size()/RAND_MAX;
      swap(positive_data_list[i],positive_data_list[j]);
    }
    for (unsigned i=0;i<negative_data_list.size();++i){
      unsigned j=rand()*negative_data_list.size()/RAND_MAX;
      swap(negative_data_list[i],negative_data_list[j]);
    }

    //oversample minority class only if there is an imbalance higher than MIN_KFOLD_IMBALANCE and if there is at least one instance for the minority class
    vector<TrainItem> balanced_positive_data_list;
    vector<TrainItem> balanced_negative_data_list;
    const double MIN_KFOLD_IMBALANCE=1;
    if (p!=0 && p<n/MIN_KFOLD_IMBALANCE) {
      cout<<"Oversampling positive factor: "<<n/(double)p<<endl;
      unsigned ratio=n/p;
      unsigned reminder=n%p;
      //duplicate a number of times equal to ratio the datastaset itself 
      for (unsigned i=0;i<ratio;i++)
	balanced_positive_data_list.insert(balanced_positive_data_list.end(),positive_data_list.begin(),positive_data_list.end());
      //add the remainder instances
      for (unsigned i=0;i<reminder;i++)
	balanced_positive_data_list.push_back(positive_data_list[i]);
      balanced_negative_data_list=negative_data_list;
    } else if (n!=0 && n<p/MIN_KFOLD_IMBALANCE) {
      cout<<"Oversampling negative factor: "<<p/(double)n<<endl;
      unsigned ratio=p/n;
      unsigned reminder=p%n;
      for (unsigned i=0;i<ratio;i++)
	balanced_negative_data_list.insert(balanced_negative_data_list.end(),negative_data_list.begin(),negative_data_list.end());
      for (unsigned i=0;i<reminder;i++)
	balanced_negative_data_list.push_back(negative_data_list[i]);
      balanced_positive_data_list=positive_data_list;
    } else {
      balanced_positive_data_list=positive_data_list;
      balanced_negative_data_list=negative_data_list;
    }

    //compose dataset by alternating positive and negative examples
    unsigned i;
    for (i=0;i<balanced_positive_data_list.size();i++){
      oDataset.push_back(balanced_positive_data_list[i]);
      if (i<balanced_negative_data_list.size()) 
	oDataset.push_back(balanced_negative_data_list[i]);
    }
    for (unsigned j=i;j<balanced_negative_data_list.size();j++)
      oDataset.push_back(balanced_negative_data_list[i]);

    //compute new class ratio 
    unsigned bp=0,bn=0;
    for (unsigned i=0;i<oDataset.size();i++) if (oDataset[i].mTarget>0) bp++; else bn++;
    cout<<"Rebalanced dataset: "<<bp+bn<<" (+:"<<bp<<" -:"<<bn<<")"<<endl;
  }

  void ThresholdInduction(const vector<double>& aMarginList, const vector<double>& aTargetList){
    if (aMarginList.size()!=aTargetList.size()) throw range_error("ERROR: Prediction list and target list cannot have different size;  prediction size:"+stream_cast<string>(aMarginList.size())+" target size:"+stream_cast<string>(aTargetList.size()));    
    
    ClearThresholds();
    
    //compute class distribution
    unsigned p=0,n=0;
    for (unsigned i=0;i<aMarginList.size();++i) 
      if (aTargetList[i]>0) p++; 
      else n++;
    //sort margins 
    vector<pair<double,unsigned> > margin_id_list;
    for (unsigned i=0;i<aMarginList.size();++i)
      margin_id_list.push_back(make_pair(aMarginList[i],i));
    sort(margin_id_list.begin(),margin_id_list.end());

    //treat special cases of lack of instances for one class
    if (p==0){
      unsigned last=margin_id_list.size()-1;
      mPositiveThreshold=margin_id_list[last].first;
      mNegativeThreshold=mPositiveThreshold;
    } else if (n==0){
      mPositiveThreshold=margin_id_list[0].first;
      mNegativeThreshold=mPositiveThreshold;
    } else { //general case
      //find discriminant such that the balanced_accuracy is above threshold
      unsigned pn_discriminant;
      for (pn_discriminant=0; pn_discriminant < margin_id_list.size(); pn_discriminant++){
	if (margin_id_list[pn_discriminant].first>=0) break;
      }

      if (PARAM.mThresholdMode=="PERFORMANCE"){
	//proceed setting the thresholds in small steps in a symmetric way around the bias value
	double STEP=.001;
	for (double th=1;th>0;th=th-STEP){
	  unsigned positive_threshold_index=margin_id_list.size() - (margin_id_list.size() - pn_discriminant)*th;
	  unsigned negative_threshold_index=(pn_discriminant)*th;
	  unsigned tp,fp,tn,fn; tp=fp=tn=fn=0;
	  for (unsigned i=positive_threshold_index;i<margin_id_list.size();++i){
	    unsigned id=margin_id_list[i].second;
	    double target=aTargetList[id];
	    if (target>0) tp++;
	    else fp++;
	  }
	  for (unsigned i=0;i<negative_threshold_index;++i){
	    unsigned id=margin_id_list[i].second;
	    double target=aTargetList[id];
	    if (target<0) tn++;
	    else fn++;
	  }

	  double pprecision=(double)tp/(tp+fp);
	  double precall=(double)tp/(tp+fn);
	  double pfmeasure=2*pprecision*precall/(pprecision+precall);

	  double nprecision=(double)tn/(tn+fn);
	  double nrecall=(double)tn/(tn+fp);
	  double nfmeasure=2*nprecision*nrecall/(nprecision+nrecall);

	  double balanced_fmeasure=(pfmeasure+nfmeasure)/2;
	  mPositiveThreshold=margin_id_list[positive_threshold_index].first;
	  mNegativeThreshold=margin_id_list[negative_threshold_index].first;
     
	  if (balanced_fmeasure>PARAM.mDataSizeThreshold) break;
	}
      } else if (PARAM.mThresholdMode=="PROPORTIONAL"){
	unsigned positive_threshold_index=margin_id_list.size() - (margin_id_list.size() - pn_discriminant)*PARAM.mDataSizeThreshold;
	unsigned negative_threshold_index=(pn_discriminant)*PARAM.mDataSizeThreshold;
	mPositiveThreshold=margin_id_list[positive_threshold_index].first;
	mNegativeThreshold=margin_id_list[negative_threshold_index].first;
      } else throw range_error("ERROR: Unknown threshold mode :"+PARAM.mThresholdMode);
    }
  }

  void SaveTrainPredictions(vector<unsigned>& aDatasetIDList,vector<double>& aMarginList, vector<double>& aPredictionList, vector<double>& aTargetList){
    string ofs_name=PARAM.mPrefix+"train_R"+stream_cast<string>(mRadius)+"_D"+stream_cast<string>(mDistance)+".predictions";
    ofstream ofs(ofs_name.c_str());
    ofs<<"# positive_threshold:"<<mPositiveThreshold<<" negative_threshold:"<<mNegativeThreshold<<endl;
    for (unsigned i=0;i<aDatasetIDList.size();++i)
      ofs<<"id: "<<aDatasetIDList[i]<<" margin: "<<aMarginList[i]<<" prediction: "<<aPredictionList[i]<<" target: "<<aTargetList[i]<<endl;
  }

  void TrainModule(vector<TrainItem>& aDataset){
#ifdef DEBUGON
    VectorClass s;
#endif
    // Shift t in order to have a reasonable initial learning rate. This assumes |x| \approx 1.
    double maxw = 1.0 / sqrt(mLambda);
    double typw = sqrt(maxw);
    double eta0 = typw / max(1.0,dloss(-typw));
    double t = 1 / (eta0 * mLambda);
  
    //Iterate epochs times in gradient descent
    ProgressBar pb(1);
    cout<<"Training for "<<mEpochs<<" epochs."<<endl;
    for (unsigned e=0; e<mEpochs; e++ ){
      pb.Count();
      //iterate over all train instances
      for (unsigned i=0;i<aDataset.size();++i){
	double eta = 1.0 / (mLambda * t);
	double s = 1 - eta * mLambda;
	mWScale *= s;
	if (mWScale < 1e-9)
	  {
	    mW.scale(mWScale);
	    mWScale = 1;
	  }
	const SVector &x = (*aDataset[i].mpInstance);
	double y = aDataset[i].mTarget;
	double wx = dot(mW,x) * mWScale;
	double margin=(wx + mBias);
	(*(aDataset[i].mpMargin))=margin;
	double z = y * margin;
#if LOSS < LOGLOSS
	if (z < 1)
#endif
	  {
	    double etd = eta * dloss(z);
	    mW.add(x, etd * y / mWScale);
#if BIAS
	    // Slower rate on the bias because it learns at each iteration.
	    mBias += etd * y * 0.01;
#endif
	  }
	t += 1;
	if (e%2==1)//every other iteration
	  if (PARAM.mIntraInstanceRegularization==true)
	    if (aDataset[i].mInstanceID!=-1)//if the instance has a correspondent graph then it has a correspondent set of feature vectors associated to each vertex and it is legal to perform intra instance regularization
	      IntraInstanceRegularization(aDataset[i].mInstanceID,eta);
      }    
#ifdef DEBUGON
      s.PushBack(mW.sparse_size());
#endif
    }
#ifdef DEBUGON
    cout<<endl<<"W size statistics: ";s.OutputStatistics(cout);cout<<endl;
#endif
  }

  void IntraInstanceRegularization(unsigned aID, double aEta){
    double ni=PARAM.mIntraInstanceRegularizationLambda/aEta;
    //for each vertex in instance with id=aID
    map<unsigned,SVector>::iterator single_graph_vertex_map_begin=mGraphVertexFeatureList[aID].begin();
    map<unsigned,SVector>::iterator single_graph_vertex_map_end=mGraphVertexFeatureList[aID].end();
    assert(single_graph_vertex_map_begin!=single_graph_vertex_map_end);
    for (map<unsigned,SVector>::iterator ij=single_graph_vertex_map_begin;ij!=single_graph_vertex_map_end;++ij){
      ////unsigned vid=ij->first;
      //extract the feature vector associated to the vertex       
      SVector& xi=ij->second;
      //build a mask vector of xi: i.e. a vector that has all features of xi set to 1
      SVector wi;
      for (unsigned j=0;j<(unsigned)xi.sparse_size();j++){
	pair<unsigned, double> feature_value=xi.extract_component(j);	  
	wi.set((int)feature_value.first,1);
      }
      //compute vector (W.xi)/|xi| wi
      double score=dot(mW,xi);
      double norm=sqrt(dot(xi,xi));
      if (norm<=0) throw range_error("Error: Instance with ID:"+stream_cast<string>(aID)+" induces a masking feature vector with norm="+stream_cast<string>(norm));
      //update W=W -ni (W- (W.xi)/|xi| wi)
      //equivalent to: W=(1-ni)W + ni (W.xi)/|xi| wi
      mW.scale(1-ni);
      wi.scale(ni*score/norm);
      mW.add(wi);
    }
  }
 
  vector<double> Test(vector<pair<double,SVector*> >& aDataset){
    vector<double> margin_list;
    //iterate over all test instances
    for (unsigned i=0;i<aDataset.size();++i){
      const SVector &x = (*aDataset[i].second);
      double y = Predict(x);
      margin_list.push_back(y);
    }
    return margin_list;      
  }

  vector<double> Test(vector<unsigned>& aTestSetIDList, vector<vector<unsigned> >& aViewPointList, DataAccessorClass* apDataset){
    //Create vector representation from graphs 
    vector<SVector> sv_data_list;
    GraphToSparseVector(aTestSetIDList,aViewPointList,apDataset,sv_data_list);
    return Test(sv_data_list);
  }

  vector<double> Test(vector<SVector>& aSVDataList){
    vector<double> margin_list;
    //iterate over all test instances
    for (unsigned i=0;i<aSVDataList.size();++i){
      const SVector &x = aSVDataList[i];
      double y = Predict(x);
      margin_list.push_back(y);
    }
    return margin_list;
  }

  inline double Predict(const SVector& x){
    return dot(mW,x) * mWScale + mBias;
  }

  vector<double> ConvertMarginToPrediction(vector<double>& aMarginList){
    vector<double> prediction_list;
    for (unsigned i=0;i<aMarginList.size();++i){
      double y = aMarginList[i];
      double prediction=0;
      if (y>=mPositiveThreshold) prediction = 1;
      else if (y<mNegativeThreshold) prediction = -1;
      else prediction=0;
      prediction_list.push_back(prediction);
    }	
    return prediction_list;
  }

  void GraphToSparseVector(vector<unsigned>& aDatasetIDList, vector<vector<unsigned> >& aViewPointList, DataAccessorClass* apDataset, vector<SVector>& oSVectorDataset){
    SetupFeatureGenerator();
    oSVectorDataset.clear();
    mGraphVertexFeatureList.clear();
    apDataset->Reset();
    vector<unsigned> graph_viewpoint_list; graph_viewpoint_list.push_back(0);
#ifdef DEBUGON
    cout<<"Converting graph to sparse feature representation (R="<<mRadius<<",D="<<mDistance<<")."<<endl;
    ProgressBar progress_bar;
    VectorClass s,ss,sparses;
#endif
    for (unsigned i=0;i<aDatasetIDList.size();++i){
      SVector x;
      oSVectorDataset.push_back(x);
      unsigned gid=aDatasetIDList[i];
      GraphClass g;
      apDataset->GetGraph(gid,g);
      g.SetGraphID(stream_cast<string>(gid));
      g.ComputePairwiseDistanceInformation(mDistance,mRadius);
      if (!PARAM.mFeatureCache) mpFeatureGenerator->Clear();
      mpFeatureGenerator->generate_feature_vector(g,oSVectorDataset[i],aViewPointList[i]);

      if (PARAM.mInstancePartsThreshold!=1 || PARAM.mIntraInstanceRegularization==true){//if viewpoint selection or intra instance regularization is required 
	//for each selected vertex extract feature vector (to be used later to score the vertex importance)
	mpFeatureGenerator->set_flag("normalization",stream_cast<string>("false"));//NOTE: do not normalize feature vector when considering parts
	//NOTE: if viewpointlist is empty then cycle for all vertices
	unsigned vsize=0;
	if (aViewPointList[i].size()>0) vsize=aViewPointList[i].size();
	else vsize=g.VertexSize();
	for (unsigned j=0;j<vsize;++j){
	  unsigned vid;
	  if (aViewPointList[i].size()>0) vid=aViewPointList[i][j];
	  else vid=j;
	  SVector xv;
	  graph_viewpoint_list[0]=vid;
	  mpFeatureGenerator->generate_feature_vector(g,xv,graph_viewpoint_list);
	  mGraphVertexFeatureList[gid][vid]=xv;
	}
	mpFeatureGenerator->set_flag("normalization",stream_cast<string>("true"));
      }
#ifdef DEBUGON
      sparses.PushBack(oSVectorDataset[i].sparse_size());
      s.PushBack(g.VertexSize());
      ss.PushBack(g.EdgeSize());
      progress_bar.Count();
#endif
    }
#ifdef DEBUGON
    cout<<endl<<"Vertex set size statistics: ";s.OutputStatistics(cout);cout<<endl;
    cout<<endl<<"Edge set size statistics: ";ss.OutputStatistics(cout);cout<<endl;
    cout<<endl<<"Sparse vector size statistics: ";sparses.OutputStatistics(cout);cout<<endl;
#endif
  }

  /**
     Ranks the vertices of each instance graph by their score and
     returns the id list of the fraction of the top and lowest scoring
     vertices for each instance.
  */
  vector<vector<unsigned> > ExtractViewpointList(vector<unsigned>& aDatasetIDList, DataAccessorClass* apDataset){
    apDataset->Reset();
    vector<vector<unsigned> > viewpoint_list;
    //for each instance extract the graph representation
    cout<<"Extracting viewpoints."<<endl;
    ProgressBar pb;
    for (unsigned i=0;i<aDatasetIDList.size();i++){
      pb.Count();
      unsigned gid=aDatasetIDList[i];
      vector<pair<double,unsigned> > score_list;
      //for each selected vertex, extract score
      map<unsigned,SVector>::iterator single_graph_vertex_map_begin=mGraphVertexFeatureList[gid].begin();
      map<unsigned,SVector>::iterator single_graph_vertex_map_end=mGraphVertexFeatureList[gid].end();
      assert(single_graph_vertex_map_begin!=single_graph_vertex_map_end);
      for (map<unsigned,SVector>::iterator ij=single_graph_vertex_map_begin;ij!=single_graph_vertex_map_end;++ij){
	unsigned vid=ij->first;
	SVector& x=ij->second;
	double score=dot(mW,x);//NOTE: no need to add bias as only the relative order is of interest here
	score_list.push_back(make_pair(score,vid));
      }
      //sort score
      sort(score_list.begin(),score_list.end());
      //extract high and low scoring vertices as a fixed fraction of all vertices
      unsigned block_size=max ((double)score_list.size()*(mInstancePartsThreshold)/2,(double)1);
      vector<unsigned> viewpoint;
      for (unsigned k=0;k<block_size;++k)
	viewpoint.push_back(score_list[k].second);
      for (unsigned k=score_list.size()-1;k>score_list.size()-1-block_size;k--)
	viewpoint.push_back(score_list[k].second);
      viewpoint_list.push_back(viewpoint);      
    }
#ifdef DEBUGON
    cout<<endl<<"Viewpoint set statistics: ";
    VectorClass s;
    for (unsigned t=0;t<viewpoint_list.size();++t) 
      s.PushBack(viewpoint_list[t].size());
    s.OutputStatistics(cout);
    cout<<endl;
#endif
    return viewpoint_list;
  }

  map<unsigned,vector<double> > ExtractVertexMargin(vector<unsigned>& aDatasetIDList, DataAccessorClass* apDataset){
    apDataset->Reset();
    mpFeatureGenerator->set_flag("normalization",stream_cast<string>("false"));//NOTE: do not normalize feature vector when considering parts
    map<unsigned,vector<double> > vertex_margin_list;
    vector<unsigned> graph_viewpoint_list; 
    graph_viewpoint_list.push_back(0);
    if (PARAM.mGraphType=="DIRECTED") 
      graph_viewpoint_list.push_back(0);//add a second viewpoint for the reverse direction graph
    //for each instance extract the graph representation
    cout<<"Extracting scores for each vertex for each instance."<<endl;
    ProgressBar pb;
    for (unsigned i=0;i<aDatasetIDList.size();i++){
      pb.Count();
      unsigned id=aDatasetIDList[i];
      GraphClass g;
      apDataset->GetGraph(id,g);      
      g.SetGraphID(stream_cast<string>(id));
      g.ComputePairwiseDistanceInformation(mDistance,mRadius);
      if (!PARAM.mFeatureCache) mpFeatureGenerator->Clear(); //NOTE: clear cache for each instance, but retain cache within vertices
      vector<double> score_list;
      //for each vertex, extract score
      unsigned size=PARAM.mGraphType=="DIRECTED"?g.VertexSize()/2:g.VertexSize();
      for (unsigned j=0;j<size;++j){      
	SVector x;
	graph_viewpoint_list[0]=j;
	if (PARAM.mGraphType=="DIRECTED") 
	  graph_viewpoint_list[1]=j+size;
	mpFeatureGenerator->generate_feature_vector(g,x,graph_viewpoint_list);
	double score=Predict(x);
	score_list.push_back(score);
      }
      vertex_margin_list[id]=score_list;
    }
    mpFeatureGenerator->set_flag("normalization",stream_cast<string>("true"));
    return vertex_margin_list;
  }


  void ClearThresholds(){mPositiveThreshold=0;mNegativeThreshold=0;}

  bool Good(){
    if (dot(mW,mW)==0) return false;
    else return true;
  }

  void OutputInfo(){
    cout<<SEP<<endl;
    cout<<"Model information"<<endl;
    cout<<SEP<<endl;
    cout<<"Bit size:"<<mBitSize<<endl;
    cout<<"Instance Parts Threshold:"<<mInstancePartsThreshold<<endl;
    cout<<"Radius:"<<mRadius<<" Distance:"<<mDistance<<endl;
    cout<<"W Norm:"<<sqrt(dot(mW,mW))*mWScale<<endl;
    cout<<"Bias:"<<mBias<<endl;
    cout<<"Positive threshold:"<<mPositiveThreshold<<endl;
    cout<<"Negative threshold:"<<mNegativeThreshold<<endl;
    cout<<SEP<<endl;
  }


  vector<unsigned> FilterData(vector<unsigned>& aDatasetIDList, DataAccessorClass* apDataset){
    unsigned bitmask = (2 << 29)-1;
    mpFeatureGenerator->set_flag("hash_bit_mask",stream_cast<string>(bitmask));//use high bit rate for more accurate hashing
    vector<vector<unsigned> > null_viewpoint_list;
    for (unsigned j=0;j<aDatasetIDList.size();++j)
      null_viewpoint_list.push_back(vector<unsigned>());
    
    vector<SVector> svector_list;
    GraphToSparseVector(aDatasetIDList, null_viewpoint_list, apDataset, svector_list);
    ComputeBinDataStructure(svector_list);
    vector<unsigned> filtered_list=ComputeMinimallyOverlappingHighDensityCenterList(svector_list);
    vector<unsigned> absolute_filtered_list;
    for (unsigned i=0;i<filtered_list.size();i++){//retrieve original id
      unsigned id=filtered_list[i];
      absolute_filtered_list.push_back(aDatasetIDList[id]);
    }
    sort(absolute_filtered_list.begin(),absolute_filtered_list.end());//sort id in increasing order
    bitmask = (2 << mBitSize)-1;
    mpFeatureGenerator->set_flag("hash_bit_mask",stream_cast<string>(bitmask));//revert to lower bit rate for efficient learning
    mBinDataStructure.clear();    //remove hash data structure from memory
    return absolute_filtered_list;
  }

  void ComputeBinDataStructure(vector<SVector>& aSVectorList){
    cout<<"Computing hash data structure..."<<endl;
    ProgressBar progress_bar;

    //init structure
    mBinDataStructure.clear();
    for (unsigned k=0;k<PARAM.mNumHashFunctions;++k)
      mBinDataStructure.push_back(multimap<unsigned,unsigned>());
    
    //fill structure
    for (unsigned i=0;i<aSVectorList.size();++i){
      vector<unsigned> min_list=ComputeHashSignature(aSVectorList[i]);

      for (unsigned k=0;k<PARAM.mNumHashFunctions;++k)
	mBinDataStructure[k].insert(make_pair(min_list[k],i));
      progress_bar.Count();
    }
 
#ifdef DEBUGON
    //output statistics 
    cout<<endl<<"Hash bin size statistics: ";
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
#endif
  }

  inline  vector<unsigned> ComputeHashSignature(SVector& aX){
    const unsigned MAXUNSIGNED=2<<30;
    vector<unsigned> signature;
    for (unsigned k=0;k<PARAM.mNumHashFunctions*PARAM.mNumCompactHashFunctions;++k)
      signature.push_back(MAXUNSIGNED);
    unsigned size=(unsigned)aX.sparse_size();
    for (unsigned f=0;f<size;++f){//scan each component    
      unsigned hash_id=aX.extract_component(f).first;
      if (hash_id==0) throw range_error("Error: Feature ID = 0. Feature ID  has to be strictly > 0");
      for (unsigned k=0;k<PARAM.mNumHashFunctions*PARAM.mNumCompactHashFunctions;++k){//for each hash function
	unsigned new_hash=IntHash(hash_id,MAXUNSIGNED,k);
	if (new_hash<signature[k]) signature[k]=new_hash;//find min hash
      }
    }
    //combine signatures to increase sparsity
    vector<unsigned> compact_signature;
    for (unsigned i=0;i<signature.size();i=i+PARAM.mNumCompactHashFunctions){
      unsigned new_hash=0;
      for(unsigned j=0;j<PARAM.mNumCompactHashFunctions;j++)
	new_hash+=signature[i+j];
      compact_signature.push_back(new_hash);
    }
    return compact_signature;
  }

  vector<vector<unsigned> > ComputeApproximateNeighborhood(vector<SVector>& aSVectorList){
    ProgressBar pb;
    cout<<"Computing approximate neighborhoods..."<<endl;
    vector<vector<unsigned> > neighborhood_list;
    for (unsigned i=0;i<aSVectorList.size();++i){
      vector<unsigned> signature=ComputeHashSignature(aSVectorList[i]);
      neighborhood_list.push_back(ComputeApproximateNeighborhood(signature));
      pb.Count();
    }
    return neighborhood_list;
  }

  vector<unsigned> ComputeApproximateNeighborhood(const vector<unsigned>& aInstanceSignature){
    map<unsigned,int> neighborhood;
    vector<pair<unsigned, double> > vec;
    for (unsigned k=0;k<PARAM.mNumHashFunctions;++k){
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

    //trim row neighborhood, i.e. given a list of neighbours with an associated occurences count, return only a fraction of the highest count ones
    const int MIN_BINS_IN_COMMON=2;//Minimum number of bins that two instances have to have in common in order to be considered similar
    vector<unsigned> neighborhood_list;
    //sort by num occurences
    vector<pair<int,unsigned> > count_list;
    for (map<unsigned,int>::const_iterator it=neighborhood.begin();it!=neighborhood.end();++it){
      unsigned id=it->first;
      int count=it->second;
      if (count>=MIN_BINS_IN_COMMON)//NOTE: consider instances that have at least MIN_BINS_IN_COMMON 
	count_list.push_back(make_pair(-count,id));//NOTE:-count to sort from highest to lowest
    }
    unsigned effective_size=min((unsigned)count_list.size(),(unsigned)(PARAM.mNeighborhoodSize));
    sort(count_list.begin(),count_list.end());
    for (unsigned i=0;i<effective_size;++i)
      neighborhood_list.push_back(count_list[i].second);
    return neighborhood_list;
  }

  double ComputeApproximateDensity(vector<unsigned>& aApproximateNeighborhood,vector<SVector>& aSVectorList){
    double density=0;
    //compute kernel pairs between all elements in aApproximateNeighborhood
    for (unsigned i=0;i<aApproximateNeighborhood.size();i++){
      unsigned ii=aApproximateNeighborhood[i];
      for (unsigned j=0;j<aApproximateNeighborhood.size();j++){
	unsigned jj=aApproximateNeighborhood[j];
	if (ii!=jj)
	  density+=dot(aSVectorList[ii],aSVectorList[jj]);
      }
    }
    unsigned size=aApproximateNeighborhood.size();
    density=density/(size*size-size);
    return density;
  }

  vector<unsigned> ComputeMinimallyOverlappingHighDensityCenterList(vector<SVector>& aSVectorList){
    vector<vector<unsigned> > neighborhood_list=ComputeApproximateNeighborhood(aSVectorList);
    //compute density for each instance in dataset
    vector<pair<double,unsigned> > density_list;
    {
      cout<<"Computing approximate density information for each instance."<<endl;
      ProgressBar progress_bar;
      for (unsigned i=0;i<aSVectorList.size();++i){
	double density=ComputeApproximateDensity(neighborhood_list[i],aSVectorList);
	density_list.push_back(make_pair(-density,i));
	progress_bar.Count();
	
      }
    }
    sort(density_list.begin(),density_list.end());
    vector<unsigned> result;
    set<unsigned> active_neighborhood;
    unsigned effective_size=density_list.size()/2;//NOTE: consider only instances with density above median value
    {
      cout<<"Selecting "<< effective_size<<" instances with nearly non-overlapping neighborhoods [max overlap allowed "<<PARAM.mMaxIntersectionSize<<" over "<<PARAM.mNeighborhoodSize<<"]."<<endl;
      ProgressBar progress_bar;
      for (unsigned i=0;i<effective_size;i++){
	unsigned id=density_list[i].second;
	set<unsigned> neighborhood_set;
	neighborhood_set.insert(neighborhood_list[id].begin(),neighborhood_list[id].end());
	set<unsigned> intersection;
	set_intersection( active_neighborhood.begin(), active_neighborhood.end(), neighborhood_set.begin(), neighborhood_set.end(), inserter(intersection,intersection.begin()) );
	//if the intersection between the neighborhood of the current center and the union of all active neighborhoods is less than a defined constant (eg. 0) then accept the new center in the active set
	if (i==0 || intersection.size()<=PARAM.mMaxIntersectionSize){
	  active_neighborhood.insert(neighborhood_list[id].begin(),neighborhood_list[id].end());
	  result.push_back(id);
	}
	progress_bar.Count(); 
      }
    }
    return result;
  }


};

//------------------------------------------------------------------------------------------------------------------------
class ModelHierarchyClass{
public:
  ModelHierarchyClass(NSPDK_FeatureGenerator* apFeatureGenerator, string aFileName):mpFeatureGenerator(apFeatureGenerator){
    Load(aFileName);
  }

  ModelHierarchyClass(NSPDK_FeatureGenerator* apFeatureGenerator, vector<pair<unsigned,unsigned> >& aLevelList,double aLambda, unsigned aEpochs, double aInstancePartsThreshold):mpFeatureGenerator(apFeatureGenerator){
    //construct all classifiers in hierarchy
    for (unsigned i=0;i<aLevelList.size();++i){
      unsigned radius=aLevelList[i].first;
      unsigned distance=aLevelList[i].second;
      ////double instance_parts_threshold=pow(aInstancePartsThreshold,(double)(i+1));//NOTE: exponentially reduce the number of parts used for more complex instance representations
      double instance_parts_threshold=aInstancePartsThreshold;
      bool abstaining_flag;
      if (i<aLevelList.size()-1) abstaining_flag=true;
      else abstaining_flag=false;
      ModelClass* mc=new ModelClass(mpFeatureGenerator,aLambda,aEpochs,PARAM.mBitSize,instance_parts_threshold,radius,distance,abstaining_flag); 
      mpModelList.push_back(mc);
    }
  }

  void Save(string aFileName){
    ofstream ofs;
    ofs.open(aFileName.c_str());
    if (!ofs) throw range_error("Cannot open file:"+aFileName);
    unsigned counter=0;
    for (unsigned i=0;i<mpModelList.size();++i)
      if (mpModelList[i]->Good()) counter++;
    ofs<<"num_levels "<<counter<<endl;
    for (unsigned i=0;i<counter;i++)
      mpModelList[i]->Save(ofs);
  }

  void Load(string aFileName){
    ifstream ifs;
    ifs.open(aFileName.c_str());
    if (!ifs) throw range_error("Cannot open file:"+aFileName);
    cout<<endl<<"Loading model file: "<<aFileName<<endl;
    unsigned num_levels=0;
    string attribute="";
    string expected="";
    ifs>>attribute>>num_levels; expected="num_levels";if (attribute!=expected) throw range_error("Format error: expecting ["+expected+"] but found ["+attribute+"]");
    for (unsigned i=0;i<num_levels;i++){
      ModelClass* mc=new ModelClass(mpFeatureGenerator);
      mc->Load(ifs);
      mc->OutputInfo();cout<<endl;
      mpModelList.push_back(mc);
    }
  }

  void Train(vector<double> aTargetList, vector<unsigned> aTrainSetIDList, DataAccessorClass* apDataset){
    assert(aTargetList.size()==aTrainSetIDList.size());
    //wrapper for semisupervised case: self-training
    //assume unsupervised material receives 0 target
    //filter the unsupervised material and put it in separate lists
    //iterate:
    //train on supervised and test on unsupervised
    //replace 0 target with prediction
    vector<double> target_list(aTargetList);
    vector<unsigned> train_supervised_id_list;
    vector<double> train_supervised_target_list;
    vector<unsigned> train_unsupervised_id_list;
    ////map<unsigned,unsigned> id_pos_map;
    for (unsigned i=0;i<target_list.size();++i){
      unsigned id=aTrainSetIDList[i];
      double target=target_list[i];
      ////id_pos_map[id]=i;
      if (target!=0) {
	train_supervised_id_list.push_back(id);
	train_supervised_target_list.push_back(target);
      } else {
	train_unsupervised_id_list.push_back(id);
      }
    }
    if (train_unsupervised_id_list.size()>0) { 
      //if unsupervised material is present then
      //train on supervised material
      cout<<endl<<"Semisupervised training on "<<target_list.size()<<" instances"<<endl;
      cout<<TAB<<"supervised instances: "<<train_supervised_id_list.size()<<" ("<<100*train_supervised_id_list.size()/(double)target_list.size() <<"%)"<<endl;
      cout<<TAB<<"unsupervised instances: "<<train_unsupervised_id_list.size()<<" ("<<100*train_unsupervised_id_list.size()/(double)target_list.size() <<"%)"<<endl;
      TrainCore(train_supervised_target_list,train_supervised_id_list,apDataset);
      
      //repeat for a predefined number of iteration
      for (unsigned iteration=0;iteration<PARAM.mSemiSupervisedIteration;iteration++){
	//test on unsupervised material
	cout<<endl<<TAB<<"Iteration "<<iteration+1<<"/"<<PARAM.mSemiSupervisedIteration<<endl;
	cout<<"Testing on unsupervised instances: "<<train_unsupervised_id_list.size()<<endl;
	vector<double> prediction_list;
	vector<double> margin_list;
	Test(train_unsupervised_id_list,apDataset,prediction_list, margin_list);
	//find high and low threshold for margin (i.e. high confidence predictions)
	vector<double> sorted_positive_margin_list;
	vector<double> sorted_negative_margin_list;
	for (unsigned i=0;i<margin_list.size();++i) 
	  if (margin_list[i]>0) sorted_positive_margin_list.push_back(margin_list[i]);
	  else sorted_negative_margin_list.push_back(margin_list[i]);
	cout<<"Predicted class distribution:  +:"<<sorted_positive_margin_list.size()<<" ("<< 100*sorted_positive_margin_list.size()/(double)train_unsupervised_id_list.size()<<" %)"<<" -:"<<sorted_negative_margin_list.size()<<" ("<< 100*sorted_negative_margin_list.size()/(double)train_unsupervised_id_list.size()<<" %)"<<endl;
	sort(sorted_positive_margin_list.begin(),sorted_positive_margin_list.end());
	sort(sorted_negative_margin_list.begin(),sorted_negative_margin_list.end());
	unsigned high_threshold_id=(sorted_positive_margin_list.size()-1)*(1-PARAM.mSemiSupervisedThreshold);
	unsigned low_threshold_id=(sorted_negative_margin_list.size()-1)*(PARAM.mSemiSupervisedThreshold);
	double high_threshold=sorted_positive_margin_list.size()>0?sorted_positive_margin_list[high_threshold_id]:sorted_negative_margin_list[sorted_negative_margin_list.size()-1];
	double low_threshold=sorted_negative_margin_list.size()>0?sorted_negative_margin_list[low_threshold_id]: sorted_positive_margin_list[0];
	cout<<"Low score threshold:"<<low_threshold<<" High score threshold:"<<high_threshold<<endl;
	//replace 0 target with predicted target only for high confidence predictions
	map<unsigned,double> semi_supervise_augmented_target_map;
	for (unsigned i=0;i<target_list.size();++i)//copy target for supervised instances
	  if(target_list[i]!=0){ 
	    unsigned id=aTrainSetIDList[i];
	    semi_supervise_augmented_target_map[id]=target_list[i];
	  }
	unsigned counter_p=0;
	unsigned counter_n=0;
	for (unsigned i=0;i<train_unsupervised_id_list.size();++i){//copy prediction for unsupervised instances
	  double margin=margin_list[i];
	  bool margin_test=false;
	  if (PARAM.mSemiSupervisedInduceOnlyPositive)
	    margin_test=(margin >= high_threshold);
	  else if (PARAM.mSemiSupervisedInduceOnlyNegative)
	    margin_test=(margin <=low_threshold);
	  else 
	    margin_test=(margin <=low_threshold  || margin >= high_threshold);
	  if (margin_test==true) {
	    double prediction=prediction_list[i];
	    unsigned id=train_unsupervised_id_list[i];
	    assert(semi_supervise_augmented_target_map.count(id)==0);
	    semi_supervise_augmented_target_map[id]=prediction;
	    if (prediction>0) counter_p++;
	    else counter_n++;
	  } 
	}
	if (PARAM.mSemiSupervisedInduceOnlyPositive) cout<<"Adding only predicted positives"<<endl;
	if (PARAM.mSemiSupervisedInduceOnlyNegative) cout<<"Adding only predicted negatives"<<endl;
	cout<<"Added +:"<<counter_p<<" and -:"<< counter_n<<" instances from unsupervised set to training set of size "<<train_supervised_id_list.size()<<endl;
	//compose indices vectors for training instances and target
	vector<unsigned> train_semi_supervise_augmented_id_list;
	vector<double> train_semi_supervise_augmented_target_list;
	for (map<unsigned,double>::iterator it=semi_supervise_augmented_target_map.begin();it!=semi_supervise_augmented_target_map.end();++it){
	  unsigned id=it->first;
	  double target=it->second;
	  train_semi_supervise_augmented_id_list.push_back(id);
	  train_semi_supervise_augmented_target_list.push_back(target);
	}
	//retrain 
	TrainCore(train_semi_supervise_augmented_target_list,train_semi_supervise_augmented_id_list,apDataset);
      }

    } else { //if no unsupervised material is present then train directly
      TrainCore(target_list,aTrainSetIDList,apDataset);
    }
  }
  void TrainCore(vector<double> aTargetList, vector<unsigned> aTrainSetIDList, DataAccessorClass* apDataset){
    unsigned model_level=0;

    vector<unsigned> trainset_id_list;
    vector<double> target_list;

    //perform initial data filtering
    if (PARAM.mDataFilter==true){
      trainset_id_list=mpModelList[model_level]->FilterData(aTrainSetIDList, apDataset);
      //sync target data
      for (unsigned i=0;i<trainset_id_list.size();++i){
	unsigned id=trainset_id_list[i];
	target_list.push_back(aTargetList[id]);
      }
    } else {
      trainset_id_list=aTrainSetIDList;
      target_list=aTargetList;
    }

    vector<vector<unsigned> > viewpoint_list;
    for (unsigned j=0;j<trainset_id_list.size();++j)
      viewpoint_list.push_back(vector<unsigned>());
    
    //perform initial training
    cout<<endl<<"Training on "<<trainset_id_list.size()<<" ("<<100*(double)trainset_id_list.size()/(double)aTrainSetIDList.size()<<"%) at level "<<model_level+1<<"/"<<mpModelList.size()<<endl;
    vector<double> margin_list;
    vector<double> prediction_list;
    mpModelList[model_level]->Train(target_list,trainset_id_list,viewpoint_list,apDataset,margin_list,prediction_list);

    vector<unsigned> tmp_trainset_id_list;
    vector<double> tmp_target_list;
    for(unsigned i=0;i<prediction_list.size();i++){
      if (prediction_list[i]==0) {
	tmp_trainset_id_list.push_back(trainset_id_list[i]);
	tmp_target_list.push_back(target_list[i]);
      } else {
	string graph_id=stream_cast<string>(trainset_id_list[i]);
	mpModelList[model_level]->mpFeatureGenerator->ClearCache(graph_id);
      }
    }
    trainset_id_list=tmp_trainset_id_list;
    target_list=tmp_target_list;

    if (trainset_id_list.size()>0){ //extract viewpoints for next level
      if (PARAM.mInstancePartsThreshold==1){//if no viewpoint selection is required then instantiate empty data structures
	viewpoint_list.clear();
	for (unsigned j=0;j<trainset_id_list.size();++j)
	  viewpoint_list.push_back(vector<unsigned>());
      } else 
	viewpoint_list=mpModelList[model_level]->ExtractViewpointList(trainset_id_list,apDataset);
    }
    model_level++;

    

    //iterate if there are instances that are ambigous
    while (trainset_id_list.size()>0 && model_level<mpModelList.size()) {
      cout<<endl<<endl;cout<<"Training on "<<trainset_id_list.size()<<" ("<<100*(double)trainset_id_list.size()/(double)aTrainSetIDList.size()<<"%) at level "<<model_level+1<<"/"<<mpModelList.size()<<endl;

      vector<double> margin_list;
      vector<double> prediction_list;
      mpModelList[model_level]->Train(target_list,trainset_id_list,viewpoint_list,apDataset,margin_list,prediction_list);

      vector<unsigned> tmp_trainset_id_list;
      vector<double> tmp_target_list;
      for(unsigned i=0;i<prediction_list.size();i++){
	if (prediction_list[i]==0) {
	  tmp_trainset_id_list.push_back(trainset_id_list[i]);
	  tmp_target_list.push_back(target_list[i]);
	} else {
	  string graph_id=stream_cast<string>(trainset_id_list[i]);
	  mpModelList[model_level]->mpFeatureGenerator->ClearCache(graph_id);
	}
      }
      trainset_id_list=tmp_trainset_id_list;
      target_list=tmp_target_list;
	
      //extract viewpoints for next level
      if (trainset_id_list.size()>0){ //extract viewpoints for next level
	if (PARAM.mInstancePartsThreshold==1){//if no viewpoint selection is required then instantiate empty data structures
	  viewpoint_list.clear();
	  for (unsigned j=0;j<trainset_id_list.size();++j)
	    viewpoint_list.push_back(vector<unsigned>());
	} else 
	  viewpoint_list=mpModelList[model_level]->ExtractViewpointList(trainset_id_list,apDataset);
      }
    
      model_level++;      
    }
  }

  void Test(vector<unsigned> aTestSetIDList,DataAccessorClass* apDataset, vector<double>& oPredictionList, vector<double>& oMarginList){
    unsigned model_level=0;
    vector<unsigned> testset_id_list=aTestSetIDList;
    vector<vector<unsigned> > viewpoint_list;
    for (unsigned j=0;j<aTestSetIDList.size();++j)
      viewpoint_list.push_back(vector<unsigned>());

    map<unsigned,double> prediction_map;
    map<unsigned,double> margin_map;

    //perform initial testing
    vector<double> margin_list=mpModelList[model_level]->Test(testset_id_list,viewpoint_list,apDataset);
    vector<double> prediction_list=mpModelList[model_level]->ConvertMarginToPrediction(margin_list);

    unsigned counter=0;
    vector<unsigned> tmp_testset_id_list;
    for(unsigned i=0;i<prediction_list.size();i++){
      if (prediction_list[i]==0) {
	tmp_testset_id_list.push_back(testset_id_list[i]);
      } else {
	prediction_map[testset_id_list[i]]=prediction_list[i];
	margin_map[testset_id_list[i]]=margin_list[i];
	counter++;
      
	string graph_id=stream_cast<string>(testset_id_list[i]);
	mpModelList[model_level]->mpFeatureGenerator->ClearCache(graph_id);
      }
    }
    cout<<endl<<"At level "<<model_level+1<<"/"<<mpModelList.size()<<endl;
    cout<<"Tested: "<<testset_id_list.size()<<endl;
    cout<<"Predicted: "<<counter<<" ( "<< 100*(double)counter/(double)testset_id_list.size()<<" % ) "<<endl;
    cout<<"Abstained: "<<testset_id_list.size()-counter<<" ( "<< 100*(double)(testset_id_list.size()-counter)/(double)testset_id_list.size()<<" % ) "<<endl;

    testset_id_list=tmp_testset_id_list;

    if (testset_id_list.size()>0){ //extract viewpoints for next level
      if (PARAM.mInstancePartsThreshold==1){//if no viewpoint selection is required then instantiate empty data structures
	viewpoint_list.clear();
	for (unsigned j=0;j<testset_id_list.size();++j)
	  viewpoint_list.push_back(vector<unsigned>());
      } else 
	viewpoint_list=mpModelList[model_level]->ExtractViewpointList(testset_id_list,apDataset);
    }
    model_level++;



    //iterate if there are instances that are ambigous
    while (testset_id_list.size()>0 && model_level<mpModelList.size()) {
      vector<double> margin_list=mpModelList[model_level]->Test(testset_id_list,viewpoint_list,apDataset);
      vector<double> prediction_list=mpModelList[model_level]->ConvertMarginToPrediction(margin_list);

      unsigned counter=0;
      vector<unsigned> tmp_testset_id_list;
      for(unsigned i=0;i<prediction_list.size();i++){
	if (prediction_list[i]==0) {
	  tmp_testset_id_list.push_back(testset_id_list[i]);
	} else {
	  prediction_map[testset_id_list[i]]=prediction_list[i];
	  margin_map[testset_id_list[i]]=margin_list[i];
	  counter++;
	
	  string graph_id=stream_cast<string>(testset_id_list[i]);
	  mpModelList[model_level]->mpFeatureGenerator->ClearCache(graph_id);
	}
      }
      cout<<endl<<"At level "<<model_level+1<<"/"<<mpModelList.size()<<endl;
      cout<<"Tested: "<<testset_id_list.size()<<endl;
      cout<<"Predicted: "<<counter<<" ( "<< 100*(double)counter/(double)testset_id_list.size()<<" % ) "<<endl;
      cout<<"Abstained: "<<testset_id_list.size()-counter<<" ( "<< 100*(double)(testset_id_list.size()-counter)/(double)testset_id_list.size()<<" % ) "<<endl;

      testset_id_list=tmp_testset_id_list;
      
      //extract viewpoints for next level
      if (testset_id_list.size()>0){ //extract viewpoints for next level
	if (PARAM.mInstancePartsThreshold==1){//if no viewpoint selection is required then instantiate empty data structures
	  viewpoint_list.clear();
	  for (unsigned j=0;j<testset_id_list.size();++j)
	    viewpoint_list.push_back(vector<unsigned>());
	} else 
	  viewpoint_list=mpModelList[model_level]->ExtractViewpointList(testset_id_list,apDataset);
      }
      
      model_level++;      
    }

    //copy predictions from map to vector
    for (map<unsigned,double>::iterator it=prediction_map.begin();it!=prediction_map.end();++it){
      //unsigned id=it->first;
      double prediction=it->second;
      oPredictionList.push_back(prediction);
    }      
    for (map<unsigned,double>::iterator it=margin_map.begin();it!=margin_map.end();++it){
      //unsigned id=it->first;
      double margin=it->second;
      oMarginList.push_back(margin);
    }      
  }

  void TestPart(vector<unsigned> aTestSetIDList,DataAccessorClass* apDataset, map<unsigned,vector<double> >& oPartMarginList){
    unsigned model_level=0;
    vector<unsigned> testset_id_list=aTestSetIDList;
    vector<vector<unsigned> > viewpoint_list;
    for (unsigned j=0;j<aTestSetIDList.size();++j)
      viewpoint_list.push_back(vector<unsigned>());

    //perform initial testing
    vector<double> margin_list=mpModelList[model_level]->Test(testset_id_list,viewpoint_list,apDataset);
    vector<double> prediction_list=mpModelList[model_level]->ConvertMarginToPrediction(margin_list);

    unsigned counter=0;
    vector<unsigned> tmp_testset_id_list;
    vector<unsigned> curr_testset_id_list;
    for(unsigned i=0;i<prediction_list.size();i++){
      if (prediction_list[i]==0) {
	tmp_testset_id_list.push_back(testset_id_list[i]);
      } else {
	curr_testset_id_list.push_back(testset_id_list[i]);
	counter++;
      }
    }
    cout<<endl<<"At level "<<model_level+1<<"/"<<mpModelList.size()<<endl;
    cout<<"Tested: "<<testset_id_list.size()<<endl;
    cout<<"Predicted: "<<counter<<" ( "<< 100*(double)counter/(double)testset_id_list.size()<<" % ) "<<endl;
    cout<<"Abstained: "<<testset_id_list.size()-counter<<" ( "<< 100*(double)(testset_id_list.size()-counter)/(double)testset_id_list.size()<<" % ) "<<endl;

    //compute vertex margin for accepted instances
    map<unsigned,vector<double> > part_margin_list;
    part_margin_list=mpModelList[model_level]->ExtractVertexMargin(curr_testset_id_list,apDataset);
    oPartMarginList.insert(part_margin_list.begin(),part_margin_list.end());

    testset_id_list=tmp_testset_id_list;

    //extract viewpoints for next level
    if (testset_id_list.size()>0){ 
      viewpoint_list=mpModelList[model_level]->ExtractViewpointList(testset_id_list,apDataset);
    }
    model_level++;

    //iterate if there are instances that are ambigous
    while (testset_id_list.size()>0 && model_level<mpModelList.size()) {
      unsigned counter=0;
      vector<unsigned> tmp_testset_id_list;
      vector<unsigned> curr_testset_id_list;
      for(unsigned i=0;i<prediction_list.size();i++){
	if (prediction_list[i]==0) {
	  tmp_testset_id_list.push_back(testset_id_list[i]);
	} else {
	  curr_testset_id_list.push_back(testset_id_list[i]);
	  counter++;
	}
      }
      cout<<endl<<"At level "<<model_level+1<<"/"<<mpModelList.size()<<endl;
      cout<<"Tested: "<<testset_id_list.size()<<endl;
      cout<<"Predicted: "<<counter<<" ( "<< 100*(double)counter/(double)testset_id_list.size()<<" % ) "<<endl;
      cout<<"Abstained: "<<testset_id_list.size()-counter<<" ( "<< 100*(double)(testset_id_list.size()-counter)/(double)testset_id_list.size()<<" % ) "<<endl;

      map<unsigned,vector<double> > part_margin_list;
      part_margin_list=mpModelList[model_level]->ExtractVertexMargin(testset_id_list,apDataset);
      oPartMarginList.insert(part_margin_list.begin(),part_margin_list.end());

      testset_id_list=tmp_testset_id_list;
      
      //extract viewpoints for next level
      if (testset_id_list.size()>0){ 
	viewpoint_list=mpModelList[model_level]->ExtractViewpointList(testset_id_list,apDataset);
      }
      model_level++;      
    }
  }
protected:
  NSPDK_FeatureGenerator* mpFeatureGenerator;
  vector<ModelClass*> mpModelList;
};


//------------------------------------------------------------------------------------------------------------------------
void LoadTargetList(const string& aFileName, vector<double>& oList){
  oList.clear();
  cout<<"Reading target file: "<<aFileName<<" ..";
  ifstream fin;
  fin.open(aFileName.c_str());
  if (!fin) throw range_error("Cannot open file:"+aFileName);
  while (!fin.eof()){
    string line;
    getline(fin,line);
    stringstream ss;
    ss<<line<<endl;
    while (!ss.eof()){
      double target(0);
      ss>>target;
      if (ss.good()){
	assert(target==1 || target==-1 || target==0);
	oList.push_back(target);
      }
    }
  }
  fin.close();
  cout<<".. read: "<<oList.size()<<" targets."<<endl;
}

//------------------------------------------------------------------------------------------------------------------------
void Train(){
  // load data
  vector<double> train_target_list;
  LoadTargetList(PARAM.mTargetFileName,train_target_list);
  DataAccessorClass* pda;
  if (PARAM.mMode=="FILE"){
    pda= new FileDataAccessorClass (PARAM.mDataFileName,PARAM.mFileType);
  } else if (PARAM.mMode=="MEMORY"){
    pda= new MemoryDataAccessorClass (PARAM.mDataFileName,PARAM.mFileType);
  } else
    throw range_error("ERROR: Unknown mode parameter: "+PARAM.mMode);

  NSPDK_FeatureGenerator* p_feature_generator=new MNSPDK_FeatureGenerator();
  if (PARAM.mKernelType=="HARD"){}
  else if (PARAM.mKernelType=="SOFT")
    p_feature_generator->set_flag("match_type","soft");
  else 
    throw range_error("ERROR: Unknown kernel type: "+PARAM.mKernelType);
#ifdef DEBUGON
  p_feature_generator->set_flag("verbosity",stream_cast<string>(1));
#endif
  ModelHierarchyClass myc(p_feature_generator,PARAM.mLevelList, PARAM.mLambda, PARAM.mEpochs, PARAM.mInstancePartsThreshold);

  cout<<endl<<SEP<<endl<<"Train phase"<<endl<<SEP<<endl;
  ProgressBar pb;pb.Count();

  vector<unsigned> train_id_list;
  for (unsigned i=0;i<train_target_list.size();++i) train_id_list.push_back(i);
  myc.Train(train_target_list,train_id_list,pda);

  string model_file_name(PARAM.mModelFileName);
  cout<<endl<<"Saving model to file: "<<model_file_name<<endl;
  myc.Save(model_file_name);
#ifdef DEBUGON
  string feature_map_file=PARAM.mPrefix+"feature.map";
  cout<<endl<<"Saving feaure map to file: "<<feature_map_file<<endl;
  ofstream ofs(feature_map_file.c_str());
  p_feature_generator->OutputFeatureMap(ofs);
  ofs.close();
#endif
  cout<<endl<<"Train phase completed:";
}

//------------------------------------------------------------------------------------------------------------------------
void Test(){
  // load data
  DataAccessorClass* pda;
  if (PARAM.mMode=="FILE"){
    pda= new FileDataAccessorClass (PARAM.mDataFileName,PARAM.mFileType);
  } else if (PARAM.mMode=="MEMORY"){
    pda= new MemoryDataAccessorClass (PARAM.mDataFileName,PARAM.mFileType);
  } else
    throw range_error("ERROR: Unknown mode parameter: "+PARAM.mMode);

  NSPDK_FeatureGenerator* p_feature_generator=new MNSPDK_FeatureGenerator();
  if (PARAM.mKernelType=="HARD"){}
  else if (PARAM.mKernelType=="SOFT")
    p_feature_generator->set_flag("match_type","soft");
  else 
    throw range_error("ERROR: Unknown kernel type: "+PARAM.mKernelType);
  ModelHierarchyClass myc(p_feature_generator,PARAM.mModelFileName);
#ifdef DEBUGON
  p_feature_generator->set_flag("verbosity",stream_cast<string>(1));
#endif

  cout<<SEP<<endl<<"Test phase"<<endl<<SEP<<endl;
  ProgressBar pb;pb.Count();
  vector<double> prediction_list;
  vector<double> margin_list;
  vector<unsigned> test_id_list;
  unsigned size=pda->Size();
  cout<<"Testing on "<<size<<" instances."<<endl;
  for (unsigned i=0;i<size;++i) test_id_list.push_back(i);
  myc.Test(test_id_list,pda,prediction_list, margin_list);
    
  string ofs_name=PARAM.mPrefix+"output.predictions";
  ofstream ofs(ofs_name.c_str());
  for (unsigned i=0;i<prediction_list.size();++i)
    ofs<<prediction_list[i]<<" "<<margin_list[i]<<endl;
  cout<<endl<<"Predictions and margins saved in file: "<<ofs_name<<endl;
  cout<<endl<<"Test phase completed:"<<endl;
}
  
//------------------------------------------------------------------------------------------------------------------------
void CrossValidation(){
  ProgressBar pbt;pbt.Count();
  // load data
  vector<double> data_target_list;
  LoadTargetList(PARAM.mTargetFileName,data_target_list);
  DataAccessorClass* pda;
  if (PARAM.mMode=="FILE"){
    pda= new FileDataAccessorClass (PARAM.mDataFileName,PARAM.mFileType);
  } else if (PARAM.mMode=="MEMORY"){
    pda= new MemoryDataAccessorClass (PARAM.mDataFileName,PARAM.mFileType);
  } else
    throw range_error("ERROR: Unknown mode parameter: "+PARAM.mMode);

  NSPDK_FeatureGenerator* p_feature_generator=new MNSPDK_FeatureGenerator();
  if (PARAM.mKernelType=="HARD"){}
  else if (PARAM.mKernelType=="SOFT")
    p_feature_generator->set_flag("match_type","soft");
  else 
    throw range_error("ERROR: Unknown kernel type: "+PARAM.mKernelType);

  //main
  vector<pair<double,double> > prediction_list;
  vector<pair<double,double> > margin_list;

  //randomly shuffle indices
  unsigned size=pda->Size();
  vector<unsigned> data_id_list;
  for (unsigned i=0;i<size;++i) data_id_list.push_back(i);
  for (unsigned i=0;i<size;++i){
    unsigned j=rand()*size/RAND_MAX;
    swap(data_id_list[i],data_id_list[j]);
  }

  //loop to build train-test split in the crossvalidation way 
  map<unsigned, vector<double> > test_result_map;
  for (unsigned f=0;f<PARAM.mCrossValidationNumFolds;f++){
    ProgressBar pbf; pbf.Count();
    cout<<SEP<<endl;
    cout<<TAB<<TAB<<"Fold: "<<f+1<<" of "<<PARAM.mCrossValidationNumFolds<<endl;
    vector<unsigned> train_id_list;
    vector<unsigned> test_id_list;
    map<unsigned,unsigned> class_counter_map;
    for (unsigned i=0;i<size;++i){
      unsigned id=data_id_list[i];
      double target=data_target_list[id];
      if (class_counter_map.count(target)==0) class_counter_map[target]=1;
      else class_counter_map[target]++;
      if (target!=0 && class_counter_map[target]%PARAM.mCrossValidationNumFolds==f) //NOTE: exclude unsupervised material from test element list 
	test_id_list.push_back(id);
      else 
	train_id_list.push_back(id);
    }
    //sort the indices in order to guarantee sequential file access
    sort(train_id_list.begin(),train_id_list.end());
    sort(test_id_list.begin(),test_id_list.end());
    //extract target list for training
    vector<double> train_target_list;
    for (unsigned i=0;i<train_id_list.size();i++){
      unsigned id=train_id_list[i];
      train_target_list.push_back(data_target_list[id]);
    }
    //extract target list for testing
    vector<double> test_target_list;
    for (unsigned i=0;i<test_id_list.size();i++){
      unsigned id=test_id_list[i];
      test_target_list.push_back(data_target_list[id]);
    }
    ModelHierarchyClass myc(p_feature_generator,PARAM.mLevelList, PARAM.mLambda, PARAM.mEpochs, PARAM.mInstancePartsThreshold);
    //perform training
    myc.Train(train_target_list,train_id_list,pda);
//perform testing
    vector<double> fold_prediction_list;
    vector<double> fold_margin_list;
    myc.Test(test_id_list,pda,fold_prediction_list,fold_margin_list);
    //add to test_result_map
    assert(fold_prediction_list.size()==fold_margin_list.size());
    assert(fold_prediction_list.size()==test_id_list.size());
    for (unsigned i=0;i<fold_prediction_list.size();i++){
      unsigned id=test_id_list[i];
      //pack all the result fields sequentially in a vector
      double target=test_target_list[i];
      double prediction=fold_prediction_list[i];
      double margin=fold_margin_list[i];
      vector<double> res;
      res.push_back(target);
      res.push_back(prediction);
      res.push_back(margin);
      //memoize the result vector with the test id
      test_result_map[id]=res;      
    }
    cout<<"Fold phase concluded in:"<<endl;
  }

  vector<double> cv_prediction_list;
  vector<double> cv_target_list;
  string ofs_name=PARAM.mPrefix+"output.cv_predictions";
  ofstream ofs(ofs_name.c_str());
  //for all test ids read in order
  for (map<unsigned,vector<double> >::iterator it=test_result_map.begin();it!=test_result_map.end();++it){
    unsigned id=it->first;
    //unpack the result fields from the result vector memoized with the test id
    double target=it->second[0];
    double prediction=it->second[1];
    double margin=it->second[2];
    ofs<<id<<" "<<target<<" "<<prediction<<" "<<margin<<endl;
    cv_prediction_list.push_back(prediction);
    cv_target_list.push_back(target);
  }
  cout<<SEP<<endl<<"Performance on data set in cross validation:"<<endl;
  OutputPerformanceMeasures(cout,cv_prediction_list,cv_target_list);

  cout<<endl<<"Instance id, true target, prediction and margin saved in file: "<<ofs_name<<endl;
  cout<<"Crossvalidation concluded in:"<<endl;
}

//------------------------------------------------------------------------------------------------------------------------
void LearningCurve(){
  cout<<SEP<<endl;
  cout<<"Computing learning curve with "<<PARAM.mLearningCurveNumPoints<<" folds."<<endl;
  ProgressBar pbt;pbt.Count();
  // load data
  vector<double> data_target_list;
  LoadTargetList(PARAM.mTargetFileName,data_target_list);
  DataAccessorClass* pda;
  if (PARAM.mMode=="FILE"){
    pda= new FileDataAccessorClass (PARAM.mDataFileName,PARAM.mFileType);
  } else if (PARAM.mMode=="MEMORY"){
    pda= new MemoryDataAccessorClass (PARAM.mDataFileName,PARAM.mFileType);
  } else
    throw range_error("ERROR: Unknown mode parameter: "+PARAM.mMode);

  NSPDK_FeatureGenerator* p_feature_generator=new MNSPDK_FeatureGenerator();
  if (PARAM.mKernelType=="HARD"){}
  else if (PARAM.mKernelType=="SOFT")
    p_feature_generator->set_flag("match_type","soft");
  else 
    throw range_error("ERROR: Unknown kernel type: "+PARAM.mKernelType);

  //main
  vector<pair<double,double> > prediction_list;
  vector<pair<double,double> > margin_list;

  //randomly shuffle indices
  unsigned size=pda->Size();
  vector<unsigned> data_id_list;
  for (unsigned i=0;i<size;++i) data_id_list.push_back(i);
  for (unsigned i=0;i<size;++i){
    unsigned j=rand()*size/RAND_MAX;
    swap(data_id_list[i],data_id_list[j]);
  }

  //build train-test split in the learning curve way 
  vector<unsigned> test_id_list;
  //test data is the first fold
  for (unsigned i=0;i<size/PARAM.mLearningCurveNumPoints;++i){
    unsigned id=data_id_list[i];
    test_id_list.push_back(id);
  }
  //sort the indices in order to guarantee sequential file access
  sort(test_id_list.begin(),test_id_list.end());
  //extract target list for testing
  vector<double> test_target_list;
  for (unsigned i=0;i<test_id_list.size();i++){
    unsigned id=test_id_list[i];
    test_target_list.push_back(data_target_list[id]);
  }
  //training data is built incrementally adding 1/PARAM.mLearningCurveNumPoints * size instances     
  for (unsigned f=1;f<PARAM.mLearningCurveNumPoints;f++){//NOTE: start from 1 as the first fold is used for the test data
    ProgressBar pbf; pbf.Count();
    cout<<SEP<<endl;
    cout<<TAB<<TAB<<"Fold: "<<f<<" of "<<PARAM.mLearningCurveNumPoints-1<<endl;
    //generate the training set
    vector<unsigned> train_id_list;
    for (unsigned i=size/PARAM.mLearningCurveNumPoints;i<size*(f+1)/PARAM.mLearningCurveNumPoints;++i){
      unsigned id=data_id_list[i];
      train_id_list.push_back(id);
    }
    //sort the indices in order to guarantee sequential file access
    sort(train_id_list.begin(),train_id_list.end());
    //extract target list for training
    vector<double> train_target_list;
    for (unsigned i=0;i<train_id_list.size();i++){
      unsigned id=train_id_list[i];
      train_target_list.push_back(data_target_list[id]);
    }
    
    ModelHierarchyClass myc(p_feature_generator,PARAM.mLevelList, PARAM.mLambda, PARAM.mEpochs, PARAM.mInstancePartsThreshold);
    //perform training
    myc.Train(train_target_list,train_id_list,pda);
    //perform testing
    vector<double> fold_prediction_list;
    vector<double> fold_margin_list;
    myc.Test(test_id_list,pda,fold_prediction_list,fold_margin_list);
    assert(fold_prediction_list.size()==fold_margin_list.size());
    assert(fold_prediction_list.size()==test_id_list.size());
    cout<<SEP<<endl<<"Performance on test set:"<<endl;
    OutputPerformanceMeasures(cout,fold_prediction_list,test_target_list);

    //save results to file
    string ofs_name=PARAM.mPrefix+"output.lc_predictions_fold"+stream_cast<string>(f);
    ofstream ofs(ofs_name.c_str());
    for (unsigned i=0;i<fold_prediction_list.size();i++){
      unsigned id=test_id_list[i];
      double target=test_target_list[i];
      double prediction=fold_prediction_list[i];
      double margin=fold_margin_list[i];
      ofs<<id<<" "<<target<<" "<<prediction<<" "<<margin<<endl;
    }
    ofs.close();
    cout<<endl<<"Instance id, true target, prediction and margin saved in file: "<<ofs_name<<endl;
    cout<<"Fold phase concluded in:"<<endl;
  }
  cout<<"Lerning curve concluded in:"<<endl;
}

//------------------------------------------------------------------------------------------------------------------------
void TestPart(){
  // load data
  DataAccessorClass* pda;
  if (PARAM.mMode=="FILE"){
    pda= new FileDataAccessorClass (PARAM.mDataFileName,PARAM.mFileType);
  } else if (PARAM.mMode=="MEMORY"){
    pda= new MemoryDataAccessorClass (PARAM.mDataFileName,PARAM.mFileType);
  } else
    throw range_error("ERROR: Unknown mode parameter: "+PARAM.mMode);
  NSPDK_FeatureGenerator* p_feature_generator=new MNSPDK_FeatureGenerator();
  if (PARAM.mKernelType=="HARD"){}
  else if (PARAM.mKernelType=="SOFT")
    p_feature_generator->set_flag("match_type","soft");
  else 
    throw range_error("ERROR: Unknown kernel type: "+PARAM.mKernelType);
  ModelHierarchyClass myc(p_feature_generator,PARAM.mModelFileName);
#ifdef DEBUGON
  p_feature_generator->set_flag("verbosity",stream_cast<string>(1));
#endif

  cout<<SEP<<endl<<"Test Part phase"<<endl<<SEP<<endl;
  ProgressBar pb;pb.Count();
  map<unsigned,vector<double> > part_margin_list;
  vector<unsigned> test_id_list;
  unsigned size=pda->Size();
  cout<<"Testing on "<<size<<" instances."<<endl;
  for (unsigned i=0;i<size;++i) test_id_list.push_back(i);
  myc.TestPart(test_id_list,pda,part_margin_list);
    
  string ofs_name=PARAM.mPrefix+"output.vertex_margins";
  ofstream ofs(ofs_name.c_str());
  for (map<unsigned,vector<double> >::iterator it=part_margin_list.begin();it!=part_margin_list.end();++it){
    unsigned instance_id=it->first;
    vector<double>& part_margin=it->second;
    for (unsigned j=0;j<part_margin.size();++j)
      ofs<<instance_id<<" "<<j<<" "<<part_margin[j]<<endl;
  }
  cout<<endl<<"Margins for each vertex for each instance saved in file: "<<ofs_name<<endl;
  cout<<endl<<"Test phase completed:"<<endl;
}

//------------------------------------------------------------------------------------------------------------------------
int main(int argc, const char **argv)
{
  cout<<SEP<<endl<<PROG_CREDIT<<endl<<SEP<<endl;
  PARAM.Init(argc, argv);
  srand(PARAM.mRandSeed);
  if (PARAM.mAction=="TRAIN") Train();
  else if (PARAM.mAction=="TEST") Test();
  else if (PARAM.mAction=="CROSSVALIDATION") CrossValidation();
  else if (PARAM.mAction=="LEARNINGCURVE") LearningCurve();
  else if (PARAM.mAction=="TEST-PART") TestPart();
  else
    throw range_error("ERROR: Unknown action parameter: "+PARAM.mAction);
  return 0;
}
