#ifndef __fnlotools__
#define __fnlotools__

#include <string>
#include <vector>
#include "fastNLOConstants.h"
#include "speaker.h"

namespace fastNLOTools {


   //! - Reading std::vectors from disk
   template<typename T> int ReadVector( std::vector<T>& v, std::istream& table , double nevts = 1);
   int ReadVector( std::vector<double>& v, std::istream& table , double nevts = 1);

   template<typename T>  int ReadFlexibleVector(std::vector<T>& v, std::istream& table, int nProcLast=0 , double nevts = 1);
   int ReadFlexibleVector( std::vector<std::string >& v, std::istream& table , int nProcLast = 0 , double nevts = 1 );
   int ReadFlexibleVector( std::vector<double >& v, std::istream& table , int nProcLast = 0 , double nevts = 1 );
   int ReadFlexibleVector( std::vector<int >& v, std::istream& table , int nProcLast = 0 , double nevts = 1 );
   int ReadFlexibleVector( std::vector<unsigned long long >& v, std::istream& table , int nProcLast = 0 , double nevts = 1 );
   int ReadUnused( std::istream& table );

   //! - Resizing tools
   template<typename T> void ResizeFlexibleVector(std::vector<T>& v, const std::vector<T>& nom);
   void ResizeFlexibleVector( std::vector<double >& v, const std::vector<double >& nom);
   void ResizeFlexibleVector( std::vector<unsigned long long >& v, const std::vector<double >& nom);

   //! - Clearing tools
   template<typename T> void ClearVector(std::vector<std::vector<T > >& v);
   template<typename T> void ClearVector(std::vector<T>& v);

   // there are nicer  options in c++11
   void ResizeVector( fastNLO::v1d& v, int dim0 );
   void ResizeVector( fastNLO::v2d& v, int dim0 , int dim1 );
   void ResizeVector( fastNLO::v3d& v, int dim0 , int dim1, int dim2 );
   void ResizeVector( fastNLO::v4d& v, int dim0 , int dim1, int dim2, int dim3 );
   void ResizeVector( fastNLO::v5d& v, int dim0 , int dim1, int dim2, int dim3, int dim4 );
   void ResizeVector( fastNLO::v6d& v, int dim0 , int dim1, int dim2, int dim3, int dim4, int dim5 );
   void ResizeVector( fastNLO::v7d& v, int dim0 , int dim1, int dim2, int dim3, int dim4, int dim5, int dim6 );

   //! - Printout of std::vectors
   template<typename T> void PrintVector( const std::vector<T>& v, std::string name, std::string prefix="");

   //! - useful i/o
   bool CheckVersion(int version); //!< check version and exit if failed.
   //bool CheckVersion(const std::string& version) {return CheckVersion((int)std::stoi(version));} ; //!< check version and exit if failed.
   bool ReadMagicNo(std::istream& table);         //!< Read and check magic number from table.
   void PutBackMagicNo(std::istream& table);      //!< Reset magic number, such that it can be recognized by other reading routines

};



//________________________________________________________________________________________________________________
// Reading functions
template<typename T>
int fastNLOTools::ReadVector( std::vector<T>& v, std::istream& table , double nevts){
   //! Read values according to the size() of the given std::vector
   //! from table (v2.0 format).
   int nn = 0;
   for( unsigned int i=0 ; i<v.size() ; i++ ){
      nn += ReadVector(v[i],table, nevts);
   }
   return nn;
};


template<typename T>
int fastNLOTools::ReadFlexibleVector(std::vector<T>& v, std::istream& table, int nProcLast, double nevts ){
   int nn = 0;
   int size = 0;
   table >> size; nn++;
   v.resize(size);
   for(unsigned int i0=0;i0<v.size();i0++){
      nn += ReadFlexibleVector(v[i0],table,nProcLast,nevts);
   }
   return nn;
};

//________________________________________________________________________________________________________________
// Resizing functions
template<typename T>
void fastNLOTools::ResizeFlexibleVector(std::vector<T>& v, const std::vector<T>& nom) {
   v.resize(nom.size());
   for (unsigned int i = 0 ; i<v.size() ; i++) {
      ResizeFlexibleVector(v[i],nom[i]);
   }
};


//________________________________________________________________________________________________________________
// Clearing
template<typename T>
void fastNLOTools::ClearVector(std::vector<std::vector<T > >& v) {
   for (unsigned int i = 0 ; i<v.size() ; i++) {
      ClearVector(v[i]);
   }
};

template<typename T>
void fastNLOTools::ClearVector(std::vector<T >& v) {
   for (unsigned int i = 0 ; i<v.size() ; i++) {
      v[i]=0;
   }
};

template<typename T>
void fastNLOTools::PrintVector( const std::vector<T>& v, std::string name, std::string prefix){
   std::cout<<" "<<prefix<<" "<<name<<std::endl;
   for(unsigned int i=0;i<v.size();i++){
      std::cout<<" "<<prefix<<"   "<<i<<"\t"<<v[i]<<std::endl;
   }
}


#endif
