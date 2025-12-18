#ifndef DiHadronCorrelationEvt_
#define DiHadronCorrelationEvt_
#include <vector>
#include <TVector3.h>
#include <TLorentzVector.h>

class DiHadronCorrelationEvt
{

 public:

  std::vector< std:: vector<TVector3> > pVect;
  std::vector< std:: vector<TVector3> > pVect_trg;
  std::vector< std:: vector<TVector3> > pVect_ass;

  std::vector< std:: vector<double> > chgVect;
  std::vector< std:: vector<double> > chgVect_trg;
  std::vector< std:: vector<double> > chgVect_ass;

  std::vector< std:: vector<double> > weightVect;
  std::vector< std:: vector<double> > wgtallVect;
  //std::vector< std:: vector<double> > weightVect_trg;
  //std::vector< std:: vector<double> > weightVect_ass;

  
  int run;
  int event;
  double zvtx;
  float cent;
  
  bool operator<(const DiHadronCorrelationEvt & b) const{
    if (zvtx != b.zvtx) return zvtx < b.zvtx;
    else return false;
  }



  DiHadronCorrelationEvt()
    {
      pVect.resize(10);
      pVect_trg.resize(10);
      pVect_ass.resize(10);

      chgVect.resize(10);
      chgVect_trg.resize(10);
      chgVect_ass.resize(10);

      weightVect.resize(10);
      wgtallVect.resize(10);
      //weightVect_trg.resize(10);
      //weightVect_ass.resize(10);
      
    }

  DiHadronCorrelationEvt(unsigned int size, unsigned int size_trg, unsigned int size_ass)
    {
      pVect.resize(size);
      chgVect.resize(size);
      weightVect.resize(size);
      wgtallVect.resize(size);
            
      pVect_trg.resize(size_trg);
      chgVect_trg.resize(size_trg);
      //weightVect_trg.resize(size_trg);
      
      pVect_ass.resize(size_ass);
      chgVect_ass.resize(size_ass);
      //weightVect_ass.resize(size_ass);
    }


  ~DiHadronCorrelationEvt()
    {
      pVect.clear();
      chgVect.clear();
      weightVect.clear();
      wgtallVect.clear();
      std::vector< std:: vector<TVector3> >().swap(pVect);
      std::vector< std:: vector<double> >().swap(chgVect);
      std::vector< std:: vector<double> >().swap(weightVect);
      std::vector< std:: vector<double> >().swap(wgtallVect);
      
      pVect_trg.clear();
      chgVect_trg.clear();
      //weightVect_trg.clear();
      std::vector< std:: vector<TVector3> >().swap(pVect_trg);
      std::vector< std:: vector<double> >().swap(chgVect_trg);
      //std::vector< std:: vector<double> >().swap(weightVect_trg);
      
      pVect_ass.clear();
      chgVect_ass.clear();
      //weightVect_ass.clear();
      std::vector< std:: vector<TVector3> >().swap(pVect_ass);
      std::vector< std:: vector<double> >().swap(chgVect_ass);
      //std::vector< std:: vector<double> >().swap(weightVect_ass);
    }


  void reset()
  {
    unsigned int size = 0;

    size = pVect.size();
    pVect.clear();
    std::vector< std:: vector<TVector3> >().swap(pVect);
    pVect.resize(size);

    size = chgVect.size();
    chgVect.clear();
    std::vector< std:: vector<double> >().swap(chgVect);
    chgVect.resize(size);

    size = weightVect.size();
    weightVect.clear();
    std::vector< std:: vector<double> >().swap(weightVect);
    weightVect.resize(size);

    size = wgtallVect.size();
    wgtallVect.clear();
    std::vector< std:: vector<double> >().swap(wgtallVect);
    wgtallVect.resize(size);

    size = pVect_trg.size();
    pVect_trg.clear();
    std::vector< std:: vector<TVector3> >().swap(pVect_trg);
    pVect_trg.resize(size);

    size = chgVect_trg.size();
    chgVect_trg.clear();
    std::vector< std:: vector<double> >().swap(chgVect_trg);
    chgVect_trg.resize(size);

    /*  size = weightVect_trg.size();
    weightVect_trg.clear();
    std::vector< std:: vector<double> >().swap(weightVect_trg);
    weightVect_trg.resize(size);
    */
    
    size = pVect_ass.size();
    pVect_ass.clear();
    std::vector< std:: vector<TVector3> >().swap(pVect_ass);
    pVect_ass.resize(size);

    size = chgVect_ass.size();
    chgVect_ass.clear();
    std::vector< std:: vector<double> >().swap(chgVect_ass);
    chgVect_ass.resize(size);

    /*
    size = weightVect_ass.size();
    weightVect_ass.clear();
    std::vector< std:: vector<double> >().swap(weightVect_ass);
    weightVect_ass.resize(size);
    */
    
    run = -999;
    event = -999;
    zvtx = -999.;
    cent = -999.;
  }

};
#endif  //DiHadronCorrelationEvt_
