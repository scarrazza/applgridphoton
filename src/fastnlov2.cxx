/**
   Interface to fastNLO grids, v2.x with fixed additive contributions
   
   Author: Giovanni Stagnitto

   We first read the whole fastNLO table, which usually contains contributions
   for different kinematic bins in rapidity. This implies that we contruct a
   vector of appl::grid object:

   >> m_grid = std::vector<appl::grid*>(Nrapidity);

   with lenght equal to the number of rapidity bins.

 **/


#include "appl_grid/fastnlo.h"
#include "appl_igrid.h"

#include "fastNLOTools.h"
#include "fastNLOConstants.h"
#include "fastNLOCoeffAddBase.h"
#include "fastNLOCoeffBase.h"

#include <iostream>
#include <fstream>
#include <istream>
#include <sstream>
#include <vector>
#include <cmath>
#include <sys/stat.h>

template<class A>
void out( A& a ) {
  std::ofstream sout( "/dev/null" );
  // std::ostream& sout = std::cout;
  sout << a << std::endl; }

template<class A>
void out( const std::string& s, A& a ) {
  std::ofstream sout( "/dev/null" );
  // std::ostream& sout = std::cout;
  sout << s << ":\t " << a << std::endl; }


typedef std::vector< std::pair<std::string, std::vector< std::pair<double, double>>> > DictVecd;

DictVecd FlattenAllButOne(std::vector < std::vector <std::pair<double,double> > >& vec) {
  
  DictVecd res;

  std::string old = "";
  std::vector<std::pair<double, double>> pts;
  
  for (unsigned int i=0; i<vec.size(); i++){
    const std::vector <std::pair<double,double>>& vi = vec[i];
    std::stringstream ss;
    for (unsigned int j=0; j<vi.size()-1; j++){
      ss << "Dim" << j+1 << ": " << vi[j].first << " - " << vi[j].second << " ";
    }
    std::string binlabel = ss.str();
    if (old == "" || binlabel == old) {
      pts.push_back(std::make_pair(vi[vi.size()-1].first, vi[vi.size()-1].second));
    } else {
      res.push_back(std::make_pair(old, pts));
      pts.clear();
      // pts.push_back( vi[vi.size()-1].first + (vi[vi.size()-1].second - vi[vi.size()-1].first)/2 );
      pts.push_back(std::make_pair(vi[vi.size()-1].first, vi[vi.size()-1].second));
    }
    old = binlabel;
  }
  
  res.push_back(std::make_pair(old, pts));

  return res;

}

// It seems that in the new version of fastNLO they have changed the order of the subprocess
int sevenpoints(int& i){
  switch (i) {
  case 0 : return 0;
  case 2 : return 6;
  case 1 : return 5;
  case 6 : return 4;
  case 4 : return 2;
  case 5 : return 3;
  case 3 : return 1;
  default : std::cerr << "Check subprocess index" << std::endl; exit(1);
  }
}


void fastnlo::constructv2( const std::string& filename ) {

  // can we open the fastnlo grid file
  struct stat stfileinfo;
  if ( stat(filename.c_str(),&stfileinfo) )   {
    std::cerr << "fastnlo::fastnlo() cannot open fastnlo grid file " << filename << std::endl;
    return;
  }

  // read the fastnlo grid
  std::ifstream table( filename.c_str() );

  std::string dummy;

  // -----------------------------------------------------------------------
  // HEADER
 
  fastNLOTools::ReadMagicNo(table);
  long long tabversion=0;   table >> tabversion; out( "tabversion", tabversion);
  std::string ScenName; table >> ScenName;   out( "ScenName", ScenName );

  int Ncontrib,Ndata,Nmult;
  table >> Ncontrib;
  table >> Nmult ; // not used any longer
  table >> Ndata;
  fastNLOTools::ReadUnused(table); // NuserString
  fastNLOTools::ReadUnused(table); // NuserInt
  fastNLOTools::ReadUnused(table); // NuserFloat
  fastNLOTools::ReadUnused(table); // Imachine
  fastNLOTools::ReadMagicNo(table);
  fastNLOTools::PutBackMagicNo(table);

  int nCoeff = Ncontrib+Ndata;
  // -----------------------------------------------------------------------

  // -----------------------------------------------------------------------
  // SCENARIO
  
  fastNLOTools::ReadMagicNo(table);

  int Ipublunits; table >> Ipublunits;   out( "Ipublunits", Ipublunits );   // # Publ. x section (10^-Ipublunits b)

  std::vector <std::string> ScDescript;
  fastNLOTools::ReadFlexibleVector(ScDescript,table);
  for ( unsigned int i=0 ; i<ScDescript.size() ; i++ ) { out( "label", ScDescript[i] );} 
  
  double Ecms; table >> Ecms; out( "Ecms", Ecms );  // # Centre-of-mass energy (Ecms/GeV)  
  int ILOord; table >> ILOord; out( "ILOord", ILOord ); // # Power in a_s of LO process (ILOord)
  unsigned int NObsBin; table >> NObsBin; out( "NObsBin", NObsBin ); // # No. of observable bins (NObsBin)
  unsigned int NDim; table >> NDim; out( "NDim", NDim ); // # Dim. of observable binning (NDim)

  std::vector <std::string> DimLabel;  
  char buffer[257];
  DimLabel.resize(NDim);
  table.getline(buffer,256);
  for(int i=NDim-1;i>=0;i--){
    table.getline(buffer,256);
    DimLabel[i] = buffer;
  }
  for ( unsigned int i=0 ; i<DimLabel.size() ; i++ ) { out( "label", DimLabel[i] );} 

  std::vector <int> IDiffBin;
  IDiffBin.resize(NDim);
  for(int i=NDim-1;i>=0;i--){  table >>  IDiffBin[i];  }
  // for ( unsigned int i=0 ; i<DimLabel.size() ; i++ )
  //   { std::cout << "IDiffBin " << i << " " << IDiffBin[i] << std::endl; }

  // Every bin has a lower and upper bin boundary and belongs to a 'dimension'.
  // In a point-wise differential measurement, the upper bin boundary is equal to the lower one.
  // Bin[i][j] contains (first, second) where first-second are limits of dim_j in bin_i
  // where i is the index of the bin and j is the index of dimension i.e. rap1, rap2, ..., pt
  std::vector < std::vector <std::pair<double,double> > > Bin;
  Bin.resize(NObsBin);
  for(unsigned int i=0;i<NObsBin;i++){
    Bin[i].resize(NDim);
    for(int j=NDim-1;j>=0;j--){
      table >> Bin[i][j].first;
      if (IDiffBin[j]==0 || IDiffBin[j]==2) {
	table >> Bin[i][j].second;
      } else {
	// For point-wise differential, IDiffBin = 1, set UpBin equal to LoBin
	Bin[i][j].second = Bin[i][j].first;
      }
      //std::cout << "Bin " << i << " " << j << " first " << Bin[i][j].first << std::endl;
      //std::cout << "Bin " << i << " " << j << " second " << Bin[i][j].second << std::endl;
    }
  }


  DictVecd RapBins = FlattenAllButOne(Bin);
  
  // for (unsigned int irap = 0; irap < RapBins.size(); irap++){
  //   out( RapBins[irap].first );
  //   for (unsigned int ipt = 0; ipt < RapBins[irap].second.size(); ipt++){
  //     std::cout << RapBins[irap].second[ipt].first << " " << RapBins[irap].second[ipt].second << std::endl;
  //   }
  // }

  // BinSize: width of each bin i.e. pthigh_i - ptlow_i
  std::vector <double> BinSize;
  fastNLOTools::ReadFlexibleVector(BinSize,table,NObsBin);
  //for (int i=0; i<BinSize.size(); i++)
  //  { std::cout << "BinSize " << i << " " << BinSize[i] << std::endl;}

  int INormFlag; table >> INormFlag; out ( "INormFlag", INormFlag );
  if (INormFlag !=0) std::cerr << "ERROR, INormFlag != 0, can have unexpected behavior" << INormFlag << std::endl;

  fastNLOTools::ReadMagicNo(table);
  fastNLOTools::PutBackMagicNo(table);


  //*****************************************************************************************
  // In these vectors we collect all the (additive) contributions we will find.
  // They will have all the same dimensions, since they will be filled at the same time.
  // If this is a number, this means that we are supposing that property to be the same for all contributions.
  
  // fnlosigmatilde[0] is the first SigmaTilde found in the file
  fastNLO::v6d fnlosigmatilde;
  // this is the vector of powers
  std::vector<int> powers;
  // number of PDFs in the initial state;
  int NPDF;
  // number of subprocess to consider;
  std::vector<int> NSubs;
  // this contains informations about the scale nodes;
  std::vector<fastNLO::v4d> SNODES;
  std::vector<fastNLO::v2d> SFACT;
  // this contains informations about the scale nodes;
  fastNLO::v2d XNODES;
  // total number of xnodes (each element of the inner vector is the nmax for each bin)
  std::vector< std::vector<int> > nxmaxs;  
  // xsections units
  std::vector< int > xsunits;  
  //weights
  fastNLO::WgtStat fWgt;
  //***********************************************************************

  
  // -----------------------------------------------------------------------
  // COEFF TABLE
  //!< read nCoeff Coefficient tables (additive, multiplicative and data)
  for (int i=0; i<nCoeff; i++) {
    
    int fVersionRead; table >> fVersionRead;
    if (fVersionRead != fastNLO::tablemagicno)
      std::cerr << "ERROR, version > 22000, can have unexpected behavior" << std::endl;
    
    int IXsectUnits;   table >> IXsectUnits;   out( "IXsectUnits", IXsectUnits );
    int IDataFlag;     table >> IDataFlag;     out( "IDataFlag", IDataFlag );
    int IAddMultFlag;  table >> IAddMultFlag;  out( "IAddMultFlag", IAddMultFlag );
    int IContrFlag1;   table >> IContrFlag1;   out( "IContrFlag1", IContrFlag1 ); // if 1, it's LO, NLO or NNLO
    int IContrFlag2;   table >> IContrFlag2;   out( "IContrFlag2", IContrFlag2 ); // 1: LO, 2: NLO, 3: NNLO
    int NScaleDep;     table >> NScaleDep;     out( "NScaleDep", NScaleDep );

    std::vector < std::string > CtrbDescript;
    std::vector < std::string > CodeDescript;
    fastNLOTools::ReadFlexibleVector(CtrbDescript,table);
    fastNLOTools::ReadFlexibleVector(CodeDescript,table);

    fastNLOTools::PrintVector(CtrbDescript, "CtrbDescript");
    fastNLOTools::PrintVector(CodeDescript, "CodeDescript");
    
    bool add = false;
    bool mult = false;
    
    if ( IDataFlag == 0 && IAddMultFlag == 0 ) {
      std::cout << "Found additive fixed order contribution (v2.0). Now reading in." << std::endl;
      add = true;
    } else if ( IDataFlag == 0 && IAddMultFlag == 1 ) {
      std::cout << "Found multiplicative contribution. Now reading in." << std::endl;
      mult = true;
    } else {
      std::cerr << "Could not identify coefficient table! (or not implemented) ... " << std::endl;
      exit(1);
    }

    if (mult) {
      // legacy reasons
      unsigned int fNObsBins = NObsBin;
      
      char buffer[5257];

      int Nuncorrel; table >> Nuncorrel;
      std::vector < std::string > UncDescr;
      UncDescr.resize(Nuncorrel);
      table.getline(buffer,5256);
      for(int i=0;i<Nuncorrel;i++){
	table.getline(buffer,5256);
	UncDescr[i] = buffer;
	//std::cout << "UncDescr " << i << " = " << UncDescr[i] << std::endl;
      }
      
      int Ncorrel; table >> Ncorrel;
      std::vector < std::string > CorDescr;
      CorDescr.resize(Ncorrel);
      table.getline(buffer,5256);
      for(int i=0;i<Ncorrel;i++){
	table.getline(buffer,5256);
	CorDescr[i] = buffer;
	//std::cout << "CorDescr " << i << " = " << CorDescr[i] << std::endl;
      }
      
      fastNLO::v2d UncorLo;
      fastNLO::v2d UncorHi;
      fastNLO::v2d CorrLo;
      fastNLO::v2d CorrHi;
      fastNLO::v1d fact;
      fact.resize(fNObsBins);
      UncorLo.resize(fNObsBins);
      UncorHi.resize(fNObsBins);
      CorrLo.resize(fNObsBins);
      CorrHi.resize(fNObsBins);
      for(unsigned int i=0;i<fNObsBins;i++){
	table >> fact[i];
	// std::cout << "fact " << i << " = " << fact[i] << std::endl;
	UncorLo[i].resize(Nuncorrel);
	UncorHi[i].resize(Nuncorrel);
	for(int j=0;j<Nuncorrel;j++){
	  table >> UncorLo[i][j];
	  table >> UncorHi[i][j];
	  // std::cout << "UncorLo " << i << " " << j << " = " << UncorLo[i][j] << std::endl;
	  // std::cout << "UncorHi " << i << " " << j << " = " << UncorHi[i][j] << std::endl;
	}
	CorrLo[i].resize(Ncorrel);
	CorrHi[i].resize(Ncorrel);
	for(int j=0;j<Ncorrel;j++){
	  table >> CorrLo[i][j];
	  table >> CorrHi[i][j];
	  // std::cout << "CorrLo " << i << " " << j << " = " << CorrLo[i][j] << std::endl;
	  // std::cout << "CorrHi " << i << " " << j << " = " << CorrHi[i][j] << std::endl;
	}
      }

    }

    if (add) {

      int stest; if ( tabversion>=24000 ) table >> stest; //"fastNLO_CoeffAddBase"
      if ( tabversion>=24000 ) fastNLOTools::ReadUnused(table);

      int IRef; table >> IRef; out( "IRef", IRef );
      int IScaleDep; table >> IScaleDep; out( "IScaleDep", IScaleDep );

      if ( tabversion >= 24000 ) {
	int Nevt; table >> Nevt;
	table >> fWgt.WgtNevt;
	table >> fWgt.NumTable;
	table >> fWgt.WgtNumEv;
	table >> fWgt.WgtSumW2;
	table >> fWgt.SigSumW2;
	table >> fWgt.SigSum;
	fastNLOTools::ReadFlexibleVector ( fWgt.WgtObsSumW2, table );
	fastNLOTools::ReadFlexibleVector ( fWgt.SigObsSumW2, table );
	fastNLOTools::ReadFlexibleVector ( fWgt.SigObsSum, table );
	fastNLOTools::ReadFlexibleVector ( fWgt.WgtObsNumEv, table );
      }

      else{
	double Nevt; table >> Nevt; out( "Nevt", Nevt );
	double readNevt = Nevt;
	if ( Nevt <= 0 ) {//v2300
	  table >> Nevt;
	  table >> fWgt.WgtNevt;

	  if ( readNevt<=-2 ) table >> fWgt.NumTable;
	  table >> fWgt.WgtNumEv;
	  table >> fWgt.WgtSumW2;
	  table >> fWgt.SigSumW2;
	  table >> fWgt.SigSum;

	  fastNLOTools::ReadFlexibleVector ( fWgt.WgtObsSumW2, table );
	  fastNLOTools::ReadFlexibleVector ( fWgt.SigObsSumW2, table );
	  fastNLOTools::ReadFlexibleVector ( fWgt.SigObsSum, table );
	  fastNLOTools::ReadFlexibleVector ( fWgt.WgtObsNumEv, table );
	}
      }
      
      int Npow; table >> Npow; out( "Npow", Npow );
      //int NPDF; initialised out of the loop
      table >> NPDF; out( "NPDF", NPDF );
				       
      std::vector < int > NPDFPDG;
      if(NPDF>0){
	NPDFPDG.resize(NPDF);
	for(int i=0;i<NPDF;i++){
	  table >>  NPDFPDG[i]; out( "NPDFPDG", NPDFPDG[i] );
	}
      }
      int NPDFDim = 0; table >> NPDFDim; out( "NPDFDim", NPDFDim );

      
      if (NPDFDim != 1){
	std::cerr << "NPDFDim != 1 (x-interpolation storage different from Half-Matrix)" << std::endl;
	std::cerr << "NPDF " << NPDF << std::endl;
	exit(1);
      } 
         
      
      std::vector < int > NFFPDG;
      int NFragFunc; table >> NFragFunc; out( "NFragFunc", NFragFunc );

      if (NFragFunc > 0){
	std::cerr << "NFragFun > 0 not yet implemented" << std::endl;
	exit(1);
      }      

      if(NFragFunc>0){
	NFFPDG.resize(NFragFunc);
	for(int i=0;i<NFragFunc;i++){
	  table >>  NFFPDG[i]; out( "NFFPDG", NFFPDG[i] );
	}
      }
      int NFFDim = 0;  table >> NFFDim; out( "NFFDim", NFFDim );
      int NSubproc = 0; table >> NSubproc; out( "NSubproc", NSubproc );

      // from fastNLOGeneratorConstants.h
      //!< Flag 1 to define PDF linear combinations of partonic subprocesses (e.g. hh --> jets: 3)
      int IPDFdef1 = 0; table >> IPDFdef1; out( "IPDFdef1", IPDFdef1 );
      //!< Flag 2 to define PDF linear combinations (dep. on IPDFdef1; for 3 e.g. 1 for jet specific LCs, 121 for generic 11x11 matrix)
      int IPDFdef2 = 0; table >> IPDFdef2; out( "IPDFdef2", IPDFdef2 );
      //!< Flag 3 to define PDF LCs at   LO (dep. on IPDFdef1, IPDFdef2; for 3, 1 e.g. 6 subprocesses, ignored for IPDFdef2==121)
      //!< Flag 3 to define PDF LCs at  NLO (dep. on IPDFdef1, IPDFdef2; for 3, 1 e.g. 7 subprocesses, ignored for IPDFdef2==121)
      //!< Flag 3 to define PDF LCs at NNLO (dep. on IPDFdef1, IPDFdef2; for 3, 1 e.g. 7 subprocesses, ignored for IPDFdef2==121)
      int IPDFdef3 = 0; table >> IPDFdef3; out( "IPDFdef3", IPDFdef3 );      

      // IF THIS CONDITION IS TRUE
      // if ( c->GetIPDFdef2()==1 &&  (c->GetIPDFdef3()==1 || c->GetIPDFdef3()==2))
      // THEN
      //  CalcPDFLinearComb is used to calculate the
      //  linear combinations of different parton flavors
      //  according to the used subprocesses of your
      //  calculation/table.
      
      if (IPDFdef1 == 0 || IPDFdef2 == 0) {
	std::cerr << "IPDFdef1 == 0 || IPDFdef2 == 0 not yet implemented" << std::endl;
	exit(1);
      }

      fastNLO::v2d XNode1, XNode2, ZNode;
      
      XNode1.resize(NObsBin);
      for(unsigned int i=0;i<NObsBin;i++){
	int xtot;
	table >> xtot;
	XNode1[i].resize(xtot);
	for(int j=0;j<xtot;j++){
	  table >> XNode1[i][j];
	  //std::cout << "XNode1 " << i << " " << j << " " << XNode1[i][j] << std::endl;
	}
      }
       
      if(NPDFDim==2){
	XNode2.resize(NObsBin);
	for(unsigned int i=0;i<NObsBin;i++){
	  int xtot;
	  table >> xtot;
	  XNode2[i].resize(xtot);
	  for(int j=0;j<xtot;j++){
            table >> XNode2[i][j];
	  }
	}
      }

      if(NFragFunc>0){
	unsigned int fNObsBins = NObsBin;
	std::vector <int> Nztot;
	Nztot.resize(fNObsBins);
	ZNode.resize(fNObsBins);
	for(int i=0;i<fNObsBins;i++){
         table >> Nztot[i];
         ZNode[i].resize(Nztot[i]);
         for(int j=0;j<Nztot[i];j++){
	   table >> ZNode[i][j];
         }
	}
      }
      
      int NScales = 0; table >> NScales; out( "NScales", NScales );
      int NScaleDim = 0; table >> NScaleDim; out( "NScaleDim", NScaleDim );
      std::vector < int > Iscale;
      Iscale.resize(NScales);
      for(int i=0;i<NScales;i++){
	table >> Iscale[i];
      }
      std::vector < std::vector < std::string > > ScaleDescript;
      int NscaleDescript;
      ScaleDescript.resize(NScaleDim);
      for(int i=0;i<NScaleDim;i++){
	table >> NscaleDescript; out( "NscaleDescript", NscaleDescript);
	ScaleDescript[i].resize(NscaleDescript);
	table.getline(buffer,256);
	for(int j=0;j<NscaleDescript;j++){
	  table.getline(buffer,256);
	  ScaleDescript[i][j] = buffer;
	  //std::cout << "ScaleDescript " << i << " " << j << " = " << ScaleDescript[i][j] << std::endl;
	}
      }

      if ( tabversion>=24000 ) fastNLOTools::ReadUnused(table);
      if ( tabversion>=24000 ) fastNLOTools::ReadUnused(table);

      std::vector < int > Nscalevar; Nscalevar.resize(NScaleDim);
      std::vector<int> Nscalenode(NScaleDim);
      for(int i=0;i<NScaleDim;i++){
	table >> Nscalevar[i];
	//std::cout << "Nscalevar " << i << " = " << Nscalevar[i] << std::endl;
	table >> Nscalenode[i];
	//std::cout << "Nscalenode " << i << " = " << Nscalenode[i] << std::endl;
      }

      fastNLO::v2d ScaleFac; ScaleFac.resize(NScaleDim);
      for(int i=0;i<NScaleDim;i++){
	ScaleFac[i].resize(Nscalevar[i]);
	for(int j=0;j<Nscalevar[i];j++){
	  table >> ScaleFac[i][j];
	  //std::cout << "ScaleFac " << i << " " << j << " = " << ScaleFac[i][j] << std::endl;
	}
      }

      fastNLO::v4d ScaleNode;
      fastNLOTools::ResizeVector( ScaleNode , NObsBin, 1 , Nscalevar[0] , Nscalenode[0] );
      int nsn = fastNLOTools::ReadVector( ScaleNode , table );
      //std::cout << "Read "<<nsn<<" lines of ScaleNode." << std::endl;
#if 0
      for (unsigned int i1 = 0; i1<ScaleNode.size(); i1++){
	for (unsigned int i2 = 0; i2<ScaleNode[i1].size(); i2++){
	  for (unsigned int i3 = 0; i3<ScaleNode[i1][i2].size(); i3++){
	    for (unsigned int i4 = 0; i4<ScaleNode[i1][i2][i3].size(); i4++){
	      std::cout << "ScaleNode " << i1 << " " << i2 << " "
			<< i3 << " "  << i4 << " " << ScaleNode[i1][i2][i3][i4] << std::endl;
	    }
	  }
	}
      }
#endif

      
      int TotalScalevars=1;
      for(int scaledim=0;scaledim<NScaleDim;scaledim++){
	TotalScalevars *= Nscalevar[scaledim];
      }

      int TotalScalenodes;
      TotalScalenodes = !ScaleNode.empty() ? ScaleNode[0][0][0].size() : 0;
      
      fastNLO::v5d SigmaTilde; // units are (p)barn * Nevt / BinSize
      SigmaTilde.resize(NObsBin);
      out("NObsBin : ", NObsBin);
      out("TotalScalevars : ", TotalScalevars);
      out("TotalScalenodes : ", TotalScalenodes);
      out("NSubproc : ", NSubproc);

      std::vector<int> nxmaxbins;
      
      for( unsigned int i=0 ; i<NObsBin ; i++ ){

	int nxmax = 0;
	switch (NPDFDim) {
	case 0: nxmax = (int)XNode1[i].size();
	  break;
	case 1: nxmax = ((int)pow((double)XNode1[i].size(),2)+XNode1[i].size())/2;
	  break;
	case 2: nxmax = XNode1[i].size()*XNode2[i].size();
	  break;
	default: ;
	}
	nxmaxbins.push_back(nxmax);
	//std::cout << "nxmax : " << nxmax << std::endl;
	
	SigmaTilde[i].resize(TotalScalevars);
	for( int k=0 ; k<TotalScalevars ; k++ ){

	  SigmaTilde[i][k].resize(TotalScalenodes);
	  for( int l=0 ; l<TotalScalenodes; l++ ){

	    SigmaTilde[i][k][l].resize(nxmax);
	    for( int m=0 ; m<nxmax ; m++ ){

	      SigmaTilde[i][k][l][m].resize(NSubproc);
	      for( int n=0 ; n<NSubproc ; n++ ){

		SigmaTilde[i][k][l][m][n] = 0.;
		//std::cout << "SigmaTilde " << i << " " << k << " " << l << " " << m << std::endl;
	      }
	    }
	  }
	}
      }
      
      int nst = fastNLOTools::ReadVector( SigmaTilde , table );
      std::cout << "Read "<<nsn+nst<<" lines of fastNLO v2 tables."<<std::endl;

      fnlosigmatilde.push_back(SigmaTilde);
      powers.push_back(Npow);
      NSubs.push_back(NSubproc);
      SFACT.push_back(ScaleFac);
      SNODES.push_back(ScaleNode);
      XNODES = XNode1;
      nxmaxs.push_back(nxmaxbins);

      xsunits.push_back(IXsectUnits);
    }
  }

  //*****************************************************************************************
  // We're now ready to create the grids

  // for legacy reason with the previous code
  int Nrapidity = RapBins.size();
  
  // This is the main object: a vector of grids, each for every rapidity bin
  m_grid = std::vector<appl::grid*>(Nrapidity);

  //*****************************************************************************************
  // Before there was a variable called 'ischeme', now not present
  // ischeme == 1 <-> 'dis' case <-> transform = f4
  // ischeme == 2 <-> 'pp and ppbar' case <-> transform = f3
  // I didn't manage to find where is this information now
  // At the moment decucing transform form NPDF
  std::string transform;
  switch ( NPDF ) { 
  case 1:  
    // dis
    transform = "f4";
    //   fy = &appl::igrid::_fy4;
    //   fx = &appl::igrid::_fx4;
    break;
  case 2: 
    // pp and ppbar  
    transform = "f3";
    //   fy = appl::igrid::_fy3;
    //   fx = appl::igrid::_fx3;
    break;
  default:
    std::cerr << "fastnlo::fastnlo() transform not defined V2" << std::endl;
    return;
  }
  out( "tranform: ", transform );

  // There was another variable before called 'ireaction'
  // with values:
  // - 1: DIS
  // - 2: ppx
  // - 3: ppbar
  // This was used to set the number of subprocess
  // Now we have directly that information in NSubProc

  // This variable is used in grid construction
  // Again, deducing this from NPDF
  // This is the name of the PDF combination to be used!
  std::string pdfname; 
  if ( NPDF == 1 ) pdfname = "dis";
  if ( NPDF == 2 ) pdfname = "nlojet";
  //if ( NPDF == 2 ) pdfname = "nlojetpp"; // ppbar

  // This variable was also used to set this boolean used in grid construction
  // Again, deducing this from NPDF
  bool DISgrid = false;
  if ( NPDF == 1 )  DISgrid = true;
  out( "DISgrid? ", DISgrid );
  //*****************************************************************************************
  //*****************************************************************************************

  // create dummy igrid, purely so we can call the fx and fy functions!
  // needed for checking purpose
  int dummysubproc = 7;
  appl::igrid grid_dummy(   5, 10, 100, 0,   5, 0.01, 0.1, 0,  transform,
			    dummysubproc, DISgrid );
  
  // this variable is useful if we want to read from vector with size = NObsBin
  int ibin = 0;
  // this is the main loop over the rapidity bins
  for (unsigned int irap=0; irap<RapBins.size(); irap++){

    std::vector<double> ptlims;
    for (unsigned int i=0; i< RapBins[irap].second.size(); i++){
      ptlims.push_back(RapBins[irap].second[i].first);
    }
    ptlims.push_back(RapBins[irap].second.back().second);

    // for (int i=0; i< ptlims.size(); i++) { out( ptlims[i] ); }
    
    m_grid[irap] = new appl::grid(ptlims, pdfname, powers[0], powers.back() - powers[0],
				  transform );

    // needed to guarantee that the x-nodes are the same as in fastnlo table
    m_grid[irap]->symmetrise(true);

    m_grid[irap]->setCMSScale(Ecms);
    m_grid[irap]->setDocumentation(RapBins[irap].first);
    //m_grid[irap]->setNormalised();
    
    for (unsigned int ipt=0; ipt<ptlims.size()-1; ipt++) {

      // need this width, since the grid is already divided by it, 
      // so when the applgrid normalises to this, it will be wrong.
      // Before it was: ptbin[irap][ipt+1]-ptbin[irap][ipt];
      // Now it's already in the table
      double width = BinSize[ibin];

      // In the previous version, the informations given were:
      // the total number of xnodes + the lower limit of the xnodes array used for each bin
      // Now we have directly access to the xnodes.
      // However, only need to make the limits, using appropriate fx/fy transforms 
      // It should be checked that by passing to appl_grid Nxtot, xl, xu
      // we should recover the values listed in the table.
           
      int Nxtot = XNODES[ibin].size();
      double xl = XNODES[ibin].front();
      double xu = XNODES[ibin].back();      


#if 0      
      // check if the upper x limit in the grid is the same as the one
      // calculated using transform, starting from xl and Nxtot
      double hylim = -grid_dummy.fy( xl );      
      double hyl = hylim;
      double hyu = hylim/Nxtot;
      double Txl = grid_dummy.fx(-hyl);
      double Txu = grid_dummy.fx(-hyu);
      // these should give the same result
      std::cout << "upper in the table  : " << xu << std::endl;
      std::cout << "upper via transform : " << Txu << std::endl;
      // and it is indeed the case
#endif
      
      for ( unsigned int iord=0 ; iord<powers.size() ; iord++ ) {

	const int& Nsubproc = NSubs[iord];

	int icentral = 0;
	for ( unsigned int imu=0; imu<SFACT[iord].size(); imu++ ) { 
	  if ( std::fabs(SFACT[iord][0][imu]-1)<1e-5 ) icentral = imu;
	} 
	
	// The same said for x values apply also for Q values	
	double Q2min = SNODES[iord][ibin][0][icentral].front();
	// std::cout << iord << " " << ibin << " " << Q2min << std::endl;
	Q2min *= Q2min;
	double Q2max = SNODES[iord][ibin][0][icentral].back();
	// std::cout << iord << " " << ibin << " " << Q2max << std::endl;
	Q2max *= Q2max;

	const unsigned int& Nscalebin = SNODES[iord][ibin][0][icentral].size();

#if 0
	std::cout << "Q2 " << Q2min
		  << "\t"  << Q2max
		  << "tau " << grid_dummy.ftau(Q2min)
		  << "\t"  << grid_dummy.ftau(Q2max)
		  << "\tNscale "  << Nscalebin 
		  << "\tDISgrid " << DISgrid << std::endl;  
#endif
	
	appl::igrid* g = new appl::igrid( Nscalebin, Q2min, Q2max, 0, Nxtot, xl, xu, 0,
					  transform, Nsubproc, DISgrid );

	g->symmetrise(true);
	
	// Where is this information???
	// However it seems it's not necessary to reweight
	//g->reweight(true);
	
	m_grid[irap]->add_igrid(ipt, iord, g);

	for ( int isub=0 ; isub<Nsubproc ; isub++ ) {

	  SparseMatrix3d& sgrid = *g->weightgrid()[isub];

	  int ix1 = 0;
	  int ix2 = 0;
	  
	  int Nx1 = g->weightgrid()[isub]->Ny();
	  int Nx2 = g->weightgrid()[isub]->Nz();
	  if ( DISgrid ) Nx2 = 1;

	  for (int ix = 0; ix < nxmaxs[iord][ibin]; ix++) {

#if 0
	    std::cout << "iord  = " << iord << ", ix = " << ix
		      << ", ix1 = " << ix1 << " (x1 = " << XNODES[ibin][Nx1-1-ix1] << ")"
		      << ", ix2 = " << ix2 << " (x2 = " << XNODES[ibin][Nx2-1-ix2] << ")"
		      << ", or x1 = " << XNODES[ibin][ix1]
		      << " and x2 = " << XNODES[ibin][ix2]
		      << std::endl;
#endif
	    
	    //for (int ix = 0; i < nxmaxs[iord]; ix++) {  
	    for ( unsigned int iQbin=0 ; iQbin<Nscalebin ; iQbin++ ) {       
			
	      // Before it was: // double d = array[ibin][ix][isub][iposition[iord]][iQbin];
	      // the second index is 0 since we want the central scale
 
	      // different indices btw applgrid and fastnlo?
	      int iappl = sevenpoints(isub);
	      //int iappl = isub;
		
	      double d = fnlosigmatilde[iord][ibin][icentral][iQbin][ix][iappl];

	      // std::cout << "\tarray[" << iord << "][" << iQbin << "][" << ix << "][" << isub << "]" << std::endl;
	      //std::cout << "  d = " << d << std::endl;
	      //std::cout << "  d*width = " << d*width << std::endl;

	      // NB!! It seems that in fastNLO 2.x the weights are referred to x*f(x),
	      // while in applgrid simpy f(x) is used. Then we multiply here by the product x1*x2.
	      double xprod = XNODES[ibin][ix1] * XNODES[ibin][ix2];
	      sgrid(iQbin, Nx2-1-ix2, Nx1-1-ix1) = d*width*xprod;
	    }

	    ix1++;
	    if (ix1 > ix2){
	      ix1 = 0;
	      ix2++;
	    }

	  }  // end of loops over x1, x2 and Q2 bins
	  
	} // end of isub loop
	
      } // end of iord loop

      ibin++;

    } // end of ipt loop

  } // end of irap loop


  // for (unsigned int idx=0; idx<m_grid.size(); idx++){
  //   std::cout << m_grid[idx]->getDocumentation() << std::endl;
  // }
  
  
}
