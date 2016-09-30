#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "DataFormats/Common/interface/Handle.h"

#include "DataFormats/Math/interface/deltaR.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/JetReco/interface/JetCollection.h"

#include <map>
#include <string>
#include <fstream>
#include <iomanip>
#include "TTree.h"
#include "Python.h"



// python script
static const std::string __script = "\
import os, numpy as np\n\
outArray = None\n\
def fillArray(*values):\n\
  print \"!!!fillArray!!!\"\n\
  global outArray\n\
  outArray = np.array(values).astype(np.float32)\n\
  return\n\
def saveFile(fout):\n\
  print \"!!!saveFile!!! \"+fout\n\
  np.savetxt(fout,outArray)\n\
  return\n\
def myprint(text):\n\
  print 'text passed: '+text+'. it works!'\n\
  return\n\
";


class MyMiniAODAnalyzer : public edm::EDAnalyzer {
  
public:
  MyMiniAODAnalyzer(const edm::ParameterSet& cfg);
  virtual ~MyMiniAODAnalyzer();
  
  virtual void analyze (const edm::Event& event, const edm::EventSetup & eventSetup);
  virtual void beginJob();
  virtual void endJob();
  virtual void beginRun(const edm::Run & run,    const edm::EventSetup & eventSetup) {};
  virtual void endRun  (const edm::Run & run,    const edm::EventSetup & eventSetup) {};
  virtual void beginLuminosityBlock  (const edm::LuminosityBlock& lumi, const edm::EventSetup& eventSetup);
  virtual void endLuminosityBlock  (const edm::LuminosityBlock& lumi, const edm::EventSetup& eventSetup);
  
private:
  
  template<typename jetCollection>
  void fillJets(const edm::Handle<jetCollection> &, std::string );
  void initialize();
  template<typename T>
  void computeQuantiles(std::vector<T>*, std::vector<T>*, std::vector<double>);
  std::map<int, std::vector<std::pair<int, int> > > readJSONFile(const std::string&);
  bool AcceptEventByRunAndLumiSection(const int& , const int& ,std::map<int, std::vector<std::pair<int, int> > >&);

  
  // // python objects
  // PyObject* _pyContext;
  // // python functions
  // PyObject* _pyPrint;

  // PyObject* _pyFillArray;      //function
  // PyObject* _pyPt;             //function arguments
  // PyObject* _pyEta;            //function arguments
  // PyObject* _pyPhi;  //function arguments
  // PyObject* _pyArray;          //array

  // PyObject* _pySaveFile;


  /// file service and tree
  edm::Service<TFileService> outfile_;

  TTree* outTree_;
  int    runId_;
  int    lumiId_;
  float  lumi_;
  int    isSig_;

  std::vector<float>* PFJetPt_;
  std::vector<float>* PFJetEta_;
  std::vector<float>* PFJetPhi_;
  std::vector<int>*   nVtx_;

  std::vector<float>* qPFJetPt_;
  std::vector<float>* qPFJetEta_;
  std::vector<float>* qPFJetPhi_;
  std::vector<int>*   qNVtx_;

  edm::EDGetTokenT<pat::JetCollection>    PFJetToken_;
  edm::EDGetTokenT<reco::VertexCollection>   vtxToken_;

  double maxJetEta_;
  double minJetPt_;

  std::string lumiFile_;
  std::map<int,std::map<int,float> > lumiMap;

  std::vector<double> quantiles_;

  std::vector<std::string> subsystemNames_;
  std::vector<bool>* subsystemQuality_;

  std::vector<std::string> qualityFiles_;

  typedef std::map<int, std::vector<std::pair<int, int> > > jsonMapType;
  jsonMapType myJsonMap;
  std::map<std::string,jsonMapType> qualityMaps;
};



MyMiniAODAnalyzer::MyMiniAODAnalyzer(const edm::ParameterSet& cfg): 
  PFJetToken_               (consumes<pat::JetCollection>(cfg.getUntrackedParameter<edm::InputTag>("PFJetTag"))),
  vtxToken_                 (consumes<reco::VertexCollection>(cfg.getUntrackedParameter<edm::InputTag>("vtx"))),
  //params for wide jet calculation
  maxJetEta_                (cfg.getUntrackedParameter<double>("maxJetEta")),
  minJetPt_                 (cfg.getUntrackedParameter<double>("minJetPt")),
  lumiFile_                 (cfg.getUntrackedParameter<std::string>("lumiFile")),
  quantiles_                (cfg.getUntrackedParameter<std::vector<double> >("quantiles")),
  subsystemNames_           (cfg.getUntrackedParameter<std::vector<std::string> >("subsystems")),
  qualityFiles_             (cfg.getUntrackedParameter<std::vector<std::string> >("qualityFiles"))
{
  // // initialize python
  // Py_Initialize();
  // PyEval_InitThreads();

  // // run the __script
  // PyObject* pyMainModule = PyImport_AddModule("__main__");
  // PyObject *pyMainDict = PyModule_GetDict(pyMainModule);
  // _pyContext = PyDict_Copy(pyMainDict);
  // PyRun_String(__script.c_str(), Py_file_input, _pyContext, _pyContext);

  // // functions
  // _pyPrint     = PyDict_GetItemString(_pyContext, "myprint");
  // _pyFillArray = PyDict_GetItemString(_pyContext, "fillArray");
  // _pySaveFile  = PyDict_GetItemString(_pyContext, "saveFile");

  
  // _pyFillArrayArgs = PyTuple_New(3);


  // //test print
  // if (PyCallable_Check(_pyPrint))
  //   {
  //     PyObject* pValue=Py_BuildValue("(z)",(char*)"something");
  //     PyErr_Print();
  //     printf("Let's give this a shot!\n");
  //     PyObject* presult=PyObject_CallObject(_pyPrint,pValue);
  //     PyErr_Print();
  //   }
  // else
  //   {
  //     std::cout << "function not callable" << std::endl;
  //     PyErr_Print();
  //   }
}

void MyMiniAODAnalyzer::initialize()
{
  lumiId_ = -1;
  lumi_  = -1;
  runId_ = -1;
  isSig_ = -1;

  PFJetPt_->clear();
  PFJetEta_->clear();
  PFJetPhi_->clear();
  nVtx_->clear();

  qPFJetPt_->clear();
  qPFJetEta_->clear();
  qPFJetPhi_->clear();
  qNVtx_->clear();

  subsystemQuality_->clear();
}

template<typename jetCollection>
void MyMiniAODAnalyzer::fillJets(const edm::Handle<jetCollection> & jets, std::string type)
{
  // Selected jets
  reco::CaloJetCollection recojets;
  typename jetCollection::const_iterator i = jets->begin();
  for(;i != jets->end(); i++){
    if(std::abs(i->eta()) < maxJetEta_ && i->pt() >= minJetPt_)
      {
	
	if(type.compare(std::string("PF")) == 0)
	  {
	    PFJetPt_->push_back(i->pt());
	    PFJetEta_->push_back(i->eta());
	    PFJetPhi_->push_back(i->phi());
	  }
	
      }
  }
  return;
}
 

template<typename T>
void MyMiniAODAnalyzer::computeQuantiles(std::vector<T>* myDistr, std::vector<T>* myQuan, std::vector<double> qq)
{
  //need to sort the distr to compute quantiles
  std::vector<T> dummyDistr = *myDistr;
  std::sort(dummyDistr.begin(), dummyDistr.end());

  //scan the vector and find quantiles
  unsigned int qItr = 0;
  for (unsigned int itr = 0; itr < dummyDistr.size(); ++itr)
    {
      if(qItr >= qq.size())
	return;

      float prob = ((float)itr+1.)/(float)dummyDistr.size();
      if(prob >= qq[qItr])
	{
	  //fill pyTuple
	  //PyTuple_SetItem(_pyFillArrayArgs, pos, PyFloat_FromDouble(dummyDistr.at(itr)));
	  //PyObject* dummyVal = Py_BuildValue("(f)",PyFloat_FromDouble(dummyDistr.at(itr)));

	  //fill root tree
	  myQuan->push_back(dummyDistr.at(itr));
	  ++qItr;
	}
    }
  return;
}






//===================== beginJob and Analyze =========================
void MyMiniAODAnalyzer::beginJob() {

  TH1::SetDefaultSumw2() ;
  outTree_ = outfile_-> make<TTree>("MyTree","MyTree");

  outTree_->Branch("runId",     &runId_,       "runId_/I");
  outTree_->Branch("lumiId",    &lumiId_,      "lumiId_/I");
  outTree_->Branch("lumi",      &lumi_,        "lumi_/F");
  outTree_->Branch("isSig",     &isSig_,       "isSig_/I");

  PFJetPt_ = new std::vector<float>;
  PFJetEta_ = new std::vector<float>;
  PFJetPhi_ = new std::vector<float>;
  nVtx_     = new std::vector<int>;
  // outTree_->Branch("PFJetPt",     "std::vector<std::float>",     &PFJetPt_);
  // outTree_->Branch("PFJetEta",    "std::vector<std::float>",     &PFJetEta_);
  // outTree_->Branch("PFJetPhi",    "std::vector<std::float>",     &PFJetPhi_);
  outTree_->Branch("nVtx",           "std::vector<std::int>",       &nVtx_);

  qPFJetPt_ = new std::vector<float>;
  qPFJetEta_ = new std::vector<float>;
  qPFJetPhi_ = new std::vector<float>;
  qNVtx_     = new std::vector<int>;
  outTree_->Branch("qPFJetPt",     "std::vector<std::float>",     &qPFJetPt_);
  outTree_->Branch("qPFJetEta",    "std::vector<std::float>",     &qPFJetEta_);
  outTree_->Branch("qPFJetPhi",    "std::vector<std::float>",     &qPFJetPhi_);
  outTree_->Branch("qNVtx",        "std::vector<std::int>",       &qNVtx_);

  subsystemQuality_ = new std::vector<bool>;
  outTree_->Branch("subsystemQuality", "std::vector<bool>",        &subsystemQuality_);
  outTree_->Branch("subsystemNames",   "std::vector<std::string>", &subsystemNames_);


  //open and load luminosity file
  std::ifstream infile;
  infile.open(lumiFile_);
  while(!infile.eof())
    {
      int run, ls;
      float lumi;
      infile >> run >> ls >> lumi;
      lumiMap[run][ls] = lumi;
    }
  infile.close();

  // Load the json files
  for(unsigned int itr = 0; itr < subsystemNames_.size(); ++itr)
    qualityMaps[subsystemNames_.at(itr)] = readJSONFile(qualityFiles_[itr]);

}

void MyMiniAODAnalyzer::endJob() 
{

  // if (PyCallable_Check(_pySaveFile))
  //   {
  //     PyObject* oName=Py_BuildValue("(z)",(char*)"outTest.txt");
  //     PyErr_Print();
  //     std::cout << "saving file..." << std::endl;
  //     PyObject* presult=PyObject_CallObject(_pySaveFile,oName);
  //     PyErr_Print();
  //   }
  // else
  //   {
  //     std::cout << "_pySaveFile function not callable" << std::endl;
  //     PyErr_Print();
  //   }


  delete PFJetPt_;
  delete PFJetEta_;
  delete PFJetPhi_;

  delete qPFJetPt_;
  delete qPFJetEta_;
  delete qPFJetPhi_;
  delete qNVtx_;

  delete subsystemQuality_;

}

MyMiniAODAnalyzer::~MyMiniAODAnalyzer()
{
  // cleanup python objects
  // Py_DECREF(_pyFillArrayArgs);
  // Py_DECREF(_pyFillArray);
  // Py_DECREF(_pySaveFile);
  // Py_DECREF(_pyPrint);
  // Py_DECREF(_pyContext);
  // Py_Finalize();
}

void MyMiniAODAnalyzer::beginLuminosityBlock (const edm::LuminosityBlock & lumi, const edm::EventSetup &eventSetup) 
{
  initialize();

  lumiId_            = lumi.luminosityBlock();
  runId_             = lumi.run();
  lumi_              = lumiMap[runId_][lumiId_];

  //store quality info
  unsigned int totQuality = 0;
  for(unsigned int itr = 0; itr < subsystemNames_.size(); ++itr)
    {
      bool qq = AcceptEventByRunAndLumiSection(runId_, lumiId_, qualityMaps[subsystemNames_.at(itr)]);
      subsystemQuality_->push_back(qq);
      if (qq == true)
	++totQuality;
    }
  if (totQuality < subsystemNames_.size())
    isSig_ = 0;
  else
    isSig_ = 1;

}

void MyMiniAODAnalyzer::endLuminosityBlock (const edm::LuminosityBlock & lumi, const edm::EventSetup &eventSetup) 
{
  //compute and store quantiles
  computeQuantiles(PFJetPt_, qPFJetPt_, quantiles_);
  computeQuantiles(PFJetEta_,qPFJetEta_,quantiles_);
  computeQuantiles(PFJetPhi_,qPFJetPhi_,quantiles_);
  computeQuantiles(nVtx_,    qNVtx_,    quantiles_);

  //fill tree one event per LS
  outTree_->Fill();

  // //fill np array one event per LS
  // if (PyCallable_Check(_pyFillArray))
  //   {
  //     _pyArray = PyObject_CallObject(_pyFillArray,_pyFillArrayArgs);
  //     PyErr_Print();
  //   }
  // else
  //   {
  //     std::cout << "_pyFillArray function not callable" << std::endl;
  //     PyErr_Print();
  //   }
}


void MyMiniAODAnalyzer::analyze (const edm::Event &event, const edm::EventSetup &eventSetup) 
{
  //fill Jets
  edm::Handle<pat::JetCollection> PFJets;
  event.getByToken(PFJetToken_,PFJets);

  if(PFJets.isValid())
    fillJets(PFJets, std::string("PF"));  

  edm::Handle<reco::VertexCollection> recVtxs;
  event.getByToken(vtxToken_,recVtxs);
  if(recVtxs.isValid())
    nVtx_->push_back(recVtxs->size());
}


//====================================================
//=== functions to parse json files the manual way ===
//====================================================
std::map<int, std::vector<std::pair<int, int> > >
MyMiniAODAnalyzer::readJSONFile(const std::string& inFileName)
{
  std::ifstream inFile(inFileName.c_str(), std::ios::in);
  
  std::string line;
  while(!inFile.eof())
    {
      std::string buffer;
      inFile >> buffer;
      line += buffer;
    }
  
  
  
  // define map with result
  std::map<int, std::vector<std::pair<int, int> > > jsonMap;
    
  
  
  // loop on JSON file
  for(std::string::const_iterator it = line.begin(); it < line.end(); ++it)
    {
      // find run number
      if( (*(it) == '"') && (*(it+7) == '"') )   
	{
	  std::string run(it+1, it+7);
	  //std::cout << "found run " << run << std::endl;
      
      
      
	  // find lumi sections
	  std::vector<std::pair<int, int> > lumisections;
	  for(std::string::const_iterator it2 = it+10; it2 < line.end(); ++it2)
	    {
	      if( (*(it2) == ']') && (*(it2-1) == ']') ) break;
	      if( *(it2) != '[' ) continue;
        
	      std::string::const_iterator it_beg = it2;
	      std::string::const_iterator it_mid;
	      std::string::const_iterator it_end;
        
	      for(std::string::const_iterator it3 = it_beg; it3 < line.end(); ++it3)
		{
		  if( *(it3) == ',' ) it_mid = it3;
		  if( *(it3) == ']' )
		    {
		      it_end = it3;
		      break;
		    }
		}
            
            
            
	      std::string lumi_beg(it_beg+1, it_mid);
	      std::string lumi_end(it_mid+1, it_end);
	      //std::cout << "[" << lumi_beg;
	      //std::cout << ",";
	      //std::cout << lumi_end << "]" << std::endl;
        
	      std::pair<int, int> tempLS(atoi(lumi_beg.c_str()), atoi(lumi_end.c_str()));
	      lumisections.push_back(tempLS);
        
	      it2 = it_end;
	    }
      
      
	  jsonMap[atoi(run.c_str())] = lumisections;
	} // find run number
    
    } // loop on JSON file
  
  
  
  return jsonMap;
}






bool MyMiniAODAnalyzer::AcceptEventByRunAndLumiSection(const int& runId, const int& lumiId,
						       std::map<int, std::vector<std::pair<int, int> > >& jsonMap)
{
  // select by runId
  if( jsonMap.find(runId) == jsonMap.end() ) return false;
  
  
  
  // select by lumiId
  std::vector<std::pair<int, int> > lumisections = jsonMap[runId];
  
  int skipEvent = true;
  for(unsigned int i = 0; i < lumisections.size(); ++i)
    if( (lumiId >= lumisections.at(i).first) &&
        (lumiId <= lumisections.at(i).second) )
      skipEvent = false;
  
  if( skipEvent == true ) return false;
  
  
  return true;
}









DEFINE_FWK_MODULE(MyMiniAODAnalyzer);
